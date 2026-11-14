#include "netcdf.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "common.hpp"
#include "csv.hpp"
#include "json.hpp"

namespace climate {
namespace {

constexpr std::uint32_t NC_DIMENSION = 10;
constexpr std::uint32_t NC_VARIABLE = 11;
constexpr std::uint32_t NC_ATTRIBUTE = 12;
constexpr std::uint32_t STREAMING = 0xFFFFFFFFu;

std::size_t type_size(int type) {
    switch (type) {
        case 1:
        case 2: return 1;
        case 3: return 2;
        case 4:
        case 5: return 4;
        case 6: return 8;
        default: throw std::runtime_error("unknown NetCDF type code " + std::to_string(type));
    }
}

std::int16_t to_i16(const unsigned char* b) {
    return static_cast<std::int16_t>((static_cast<std::uint16_t>(b[0]) << 8) |
                                     static_cast<std::uint16_t>(b[1]));
}

std::int32_t to_i32(const unsigned char* b) {
    return static_cast<std::int32_t>((static_cast<std::uint32_t>(b[0]) << 24) |
                                     (static_cast<std::uint32_t>(b[1]) << 16) |
                                     (static_cast<std::uint32_t>(b[2]) << 8) |
                                     static_cast<std::uint32_t>(b[3]));
}

double to_f32(const unsigned char* b) {
    const std::uint32_t bits = static_cast<std::uint32_t>(to_i32(b));
    float value = 0.0f;
    std::memcpy(&value, &bits, sizeof(value));
    return static_cast<double>(value);
}

double to_f64(const unsigned char* b) {
    std::uint64_t bits = 0;
    for (int i = 0; i < 8; ++i) bits = (bits << 8) | static_cast<std::uint64_t>(b[i]);
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

class HeaderReader {
public:
    explicit HeaderReader(std::ifstream& input) : input_(input) {}

    std::uint32_t u32() {
        unsigned char b[4];
        read(b, 4);
        return (static_cast<std::uint32_t>(b[0]) << 24) |
               (static_cast<std::uint32_t>(b[1]) << 16) |
               (static_cast<std::uint32_t>(b[2]) << 8) | static_cast<std::uint32_t>(b[3]);
    }

    std::uint64_t u64() {
        const std::uint64_t high = u32();
        const std::uint64_t low = u32();
        return (high << 32) | low;
    }

    std::string name() {
        const std::uint32_t length = u32();
        std::string text(length, '\0');
        if (length > 0) read(reinterpret_cast<unsigned char*>(&text[0]), length);
        align();
        return text;
    }

    std::vector<double> numbers(int type, std::size_t count) {
        const std::size_t size = type_size(type);
        std::vector<unsigned char> raw(count * size);
        if (!raw.empty()) read(raw.data(), raw.size());
        std::vector<double> values;
        values.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            const unsigned char* element = raw.data() + i * size;
            switch (type) {
                case 1: values.push_back(static_cast<double>(static_cast<std::int8_t>(element[0]))); break;
                case 2: values.push_back(static_cast<double>(element[0])); break;
                case 3: values.push_back(static_cast<double>(to_i16(element))); break;
                case 4: values.push_back(static_cast<double>(to_i32(element))); break;
                case 5: values.push_back(to_f32(element)); break;
                case 6: values.push_back(to_f64(element)); break;
                default: throw std::runtime_error("unsupported NetCDF type");
            }
        }
        align();
        return values;
    }

    std::string chars(std::size_t count) {
        std::string text(count, '\0');
        if (count > 0) read(reinterpret_cast<unsigned char*>(&text[0]), count);
        align();
        return text;
    }

    void align() {
        const std::streamoff position = input_.tellg();
        if (position < 0) throw std::runtime_error("cannot determine NetCDF header position");
        const std::streamoff remainder = position % 4;
        if (remainder != 0) input_.seekg(4 - remainder, std::ios::cur);
    }

private:
    void read(unsigned char* destination, std::size_t count) {
        input_.read(reinterpret_cast<char*>(destination), static_cast<std::streamsize>(count));
        if (!input_) throw std::runtime_error("unexpected end of NetCDF header");
    }

    std::ifstream& input_;
};

bool read_list_header(HeaderReader& reader, std::uint32_t expected_tag, std::uint32_t& count) {
    const std::uint32_t tag = reader.u32();
    count = reader.u32();
    if (tag == 0 && count == 0) return false;
    if (tag != expected_tag)
        throw std::runtime_error("malformed NetCDF list header (tag " + std::to_string(tag) + ")");
    return true;
}

NcAttribute read_attribute(HeaderReader& reader) {
    NcAttribute attribute;
    attribute.name = reader.name();
    attribute.type = static_cast<int>(reader.u32());
    attribute.count = static_cast<std::size_t>(reader.u32());
    if (attribute.type == 2) {
        attribute.text = reader.chars(attribute.count);
    } else {
        attribute.numbers = reader.numbers(attribute.type, attribute.count);
    }
    return attribute;
}

std::vector<NcAttribute> read_attributes(HeaderReader& reader) {
    std::vector<NcAttribute> attributes;
    std::uint32_t count = 0;
    if (!read_list_header(reader, NC_ATTRIBUTE, count)) return attributes;
    attributes.reserve(count);
    for (std::uint32_t i = 0; i < count; ++i) attributes.push_back(read_attribute(reader));
    return attributes;
}

std::size_t product_from(const std::vector<std::size_t>& shape, std::size_t start) {
    std::size_t total = 1;
    for (std::size_t i = start; i < shape.size(); ++i) total *= shape[i];
    return total;
}

std::string lower_copy(std::string text) {
    for (auto& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return text;
}

bool contains(const std::string& text, const std::string& needle) {
    return text.find(needle) != std::string::npos;
}

const NcAttribute* find_attribute(const NcVariable& variable, const std::string& name) {
    for (const auto& attribute : variable.attributes)
        if (attribute.name == name) return &attribute;
    return nullptr;
}

double numeric_attribute(const NcVariable& variable, const std::string& name, double fallback,
                         bool* found) {
    const NcAttribute* attribute = find_attribute(variable, name);
    if (attribute == nullptr || attribute->numbers.empty()) {
        if (found != nullptr) *found = false;
        return fallback;
    }
    if (found != nullptr) *found = true;
    return attribute->numbers.front();
}

// --- calendar helpers (proleptic Gregorian, days since 1970-01-01) ---

long long days_from_civil(int year, unsigned month, unsigned day) {
    year -= month <= 2;
    const long long era = (year >= 0 ? year : year - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(year - era * 400);
    const unsigned doy = (153u * (month + (month > 2 ? -3u : 9u)) + 2u) / 5u + day - 1u;
    const unsigned doe = yoe * 365u + yoe / 4u - yoe / 100u + doy;
    return era * 146097ll + static_cast<long long>(doe) - 719468ll;
}

struct CivilDate {
    int year = 1970;
    int month = 1;
    int day = 1;
};

CivilDate civil_from_days(long long z) {
    z += 719468ll;
    const long long era = (z >= 0 ? z : z - 146096ll) / 146097ll;
    const unsigned doe = static_cast<unsigned>(z - era * 146097ll);
    const unsigned yoe = (doe - doe / 1460u + doe / 36524u - doe / 146096u) / 365u;
    const long long year = static_cast<long long>(yoe) + era * 400ll;
    const unsigned doy = doe - (365u * yoe + yoe / 4u - yoe / 100u);
    const unsigned mp = (5u * doy + 2u) / 153u;
    const unsigned day = doy - (153u * mp + 2u) / 5u + 1u;
    const unsigned month = mp + (mp < 10u ? 3u : static_cast<unsigned>(-9));
    CivilDate date;
    date.year = static_cast<int>(year + (month <= 2u ? 1 : 0));
    date.month = static_cast<int>(month);
    date.day = static_cast<int>(day);
    return date;
}

bool is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

int days_in_year(int year) { return is_leap_year(year) ? 366 : 365; }

bool parse_date(const std::string& text, CivilDate& date) {
    int year = 0;
    int month = 0;
    int day = 0;
    if (std::sscanf(text.c_str(), "%d-%d-%d", &year, &month, &day) != 3) return false;
    if (month < 1 || month > 12 || day < 1 || day > 31) return false;
    date.year = year;
    date.month = month;
    date.day = day;
    return true;
}

}  // namespace

std::string netcdf_type_name(int type) {
    switch (type) {
        case 1: return "byte";
        case 2: return "char";
        case 3: return "short";
        case 4: return "int";
        case 5: return "float";
        case 6: return "double";
        default: return "type" + std::to_string(type);
    }
}

bool is_hdf5_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) return false;
    unsigned char magic[8] = {0};
    input.read(reinterpret_cast<char*>(magic), 8);
    static const unsigned char expected[8] = {0x89, 'H', 'D', 'F', '\r', '\n', 0x1A, '\n'};
    return std::memcmp(magic, expected, 8) == 0;
}

NetCDFFile NetCDFFile::open(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open NetCDF file: " + path);

    unsigned char magic[4];
    input.read(reinterpret_cast<char*>(magic), 4);
    if (!input) throw std::runtime_error("file is too short to be NetCDF: " + path);
    if (!(magic[0] == 'C' && magic[1] == 'D' && magic[2] == 'F'))
        throw std::runtime_error("not a NetCDF-3 classic file: " + path +
                                 " (NetCDF-4/HDF5 files must be converted first)");
    if (magic[3] != 1 && magic[3] != 2)
        throw std::runtime_error("unsupported NetCDF version byte " +
                                 std::to_string(static_cast<int>(magic[3])));

    NetCDFFile file;
    file.path_ = path;
    file.version_ = static_cast<int>(magic[3]);
    HeaderReader reader(input);

    const std::uint32_t numrecs = reader.u32();
    file.record_count_ = numrecs == STREAMING ? 0 : numrecs;

    std::uint32_t count = 0;
    if (read_list_header(reader, NC_DIMENSION, count)) {
        file.dimensions_.reserve(count);
        for (std::uint32_t i = 0; i < count; ++i) {
            NcDimension dimension;
            dimension.name = reader.name();
            dimension.length = reader.u32();
            dimension.unlimited = dimension.length == 0;
            file.dimensions_.push_back(dimension);
        }
    }

    file.global_attributes_ = read_attributes(reader);

    if (read_list_header(reader, NC_VARIABLE, count)) {
        file.variables_.reserve(count);
        for (std::uint32_t i = 0; i < count; ++i) {
            NcVariable variable;
            variable.name = reader.name();
            const std::uint32_t ndims = reader.u32();
            variable.dimension_ids.resize(ndims);
            for (std::uint32_t d = 0; d < ndims; ++d)
                variable.dimension_ids[d] = static_cast<int>(reader.u32());
            variable.attributes = read_attributes(reader);
            variable.type = static_cast<int>(reader.u32());
            variable.vsize = reader.u32();
            variable.data_offset = file.version_ == 1 ? reader.u32() : reader.u64();

            variable.shape.resize(ndims);
            for (std::uint32_t d = 0; d < ndims; ++d) {
                const int id = variable.dimension_ids[d];
                if (id < 0 || static_cast<std::size_t>(id) >= file.dimensions_.size())
                    throw std::runtime_error("variable '" + variable.name +
                                             "' references an invalid dimension id");
                const NcDimension& dimension = file.dimensions_[static_cast<std::size_t>(id)];
                if (dimension.unlimited && d == 0) {
                    variable.record = true;
                    variable.shape[d] = file.record_count_;
                } else {
                    variable.shape[d] = dimension.length;
                }
            }
            variable.element_count = product_from(variable.shape, 0);
            file.variables_.push_back(variable);
        }
    }

    // Record stride is the sum of the per-record sizes of all record variables.
    for (const auto& variable : file.variables_)
        if (variable.record) file.record_stride_ += variable.vsize;

    if (numrecs == STREAMING && file.record_stride_ > 0) {
        input.seekg(0, std::ios::end);
        const std::streamoff size = input.tellg();
        std::uint64_t data_start = 0;
        for (const auto& variable : file.variables_)
            if (variable.record)
                data_start = data_start == 0 ? variable.data_offset
                                             : std::min<std::uint64_t>(data_start,
                                                                       variable.data_offset);
        if (data_start > 0 && size > static_cast<std::streamoff>(data_start)) {
            file.record_count_ = static_cast<std::size_t>(
                (static_cast<std::uint64_t>(size) - data_start) / file.record_stride_);
            for (auto& variable : file.variables_)
                if (variable.record && !variable.shape.empty())
                    variable.shape[0] = file.record_count_;
            for (auto& variable : file.variables_)
                variable.element_count = product_from(variable.shape, 0);
        }
    }

    return file;
}

const NcVariable* NetCDFFile::find_variable(const std::string& name) const {
    for (const auto& variable : variables_)
        if (variable.name == name) return &variable;
    const std::string lowered = lower_copy(name);
    for (const auto& variable : variables_)
        if (lower_copy(variable.name) == lowered) return &variable;
    return nullptr;
}

double NetCDFFile::attribute_number(const NcVariable& variable, const std::string& name,
                                    double fallback) const {
    return numeric_attribute(variable, name, fallback, nullptr);
}

std::string NetCDFFile::attribute_text(const NcVariable& variable, const std::string& name) const {
    const NcAttribute* attribute = find_attribute(variable, name);
    if (attribute == nullptr) return std::string();
    return attribute->text;
}

std::vector<double> NetCDFFile::read_variable(const std::string& name) const {
    const NcVariable* variable = find_variable(name);
    if (variable == nullptr) throw std::runtime_error("NetCDF variable not found: " + name);

    std::ifstream input(path_, std::ios::binary);
    if (!input) throw std::runtime_error("cannot reopen NetCDF file: " + path_);

    const std::size_t size = type_size(variable->type);
    const std::size_t per_record =
        variable->record ? product_from(variable->shape, 1) : variable->element_count;
    const std::size_t records = variable->record ? record_count_ : 1;

    std::vector<double> values;
    values.reserve(per_record * records);
    std::vector<unsigned char> buffer(per_record * size);

    for (std::size_t record = 0; record < records; ++record) {
        const std::uint64_t offset =
            variable->data_offset +
            (variable->record ? static_cast<std::uint64_t>(record) * record_stride_ : 0);
        input.seekg(static_cast<std::streamoff>(offset));
        if (!input) throw std::runtime_error("cannot seek to NetCDF variable data");
        if (!buffer.empty())
            input.read(reinterpret_cast<char*>(buffer.data()),
                       static_cast<std::streamsize>(buffer.size()));
        if (!input) throw std::runtime_error("truncated NetCDF variable data: " + name);

        for (std::size_t i = 0; i < per_record; ++i) {
            const unsigned char* element = buffer.data() + i * size;
            switch (variable->type) {
                case 1:
                    values.push_back(static_cast<double>(static_cast<std::int8_t>(element[0])));
                    break;
                case 2: values.push_back(static_cast<double>(element[0])); break;
                case 3: values.push_back(static_cast<double>(to_i16(element))); break;
                case 4: values.push_back(static_cast<double>(to_i32(element))); break;
                case 5: values.push_back(to_f32(element)); break;
                case 6: values.push_back(to_f64(element)); break;
                default: throw std::runtime_error("unsupported NetCDF variable type");
            }
        }
    }
    return values;
}

double decimal_year_from_time(double value, const std::string& units) {
    const std::string lowered = lower_copy(units);
    const std::size_t since = lowered.find(" since ");
    if (since == std::string::npos)
        throw std::runtime_error("time units must contain 'since': " + units);

    const std::string unit = lowered.substr(0, since);
    const std::string epoch_text = lowered.substr(since + 7);

    CivilDate epoch;
    const std::string date_part = epoch_text.substr(0, epoch_text.find(' '));
    if (!parse_date(date_part, epoch))
        throw std::runtime_error("cannot parse time epoch: " + epoch_text);

    const long long epoch_days = days_from_civil(epoch.year, static_cast<unsigned>(epoch.month),
                                                 static_cast<unsigned>(epoch.day));

    if (contains(unit, "day") || contains(unit, "hour") || contains(unit, "minute") ||
        contains(unit, "second")) {
        double days = value;
        if (contains(unit, "hour")) days = value / 24.0;
        if (contains(unit, "minute")) days = value / 1440.0;
        if (contains(unit, "second")) days = value / 86400.0;

        const double absolute = static_cast<double>(epoch_days) + days;
        const long long whole_days = static_cast<long long>(std::floor(absolute));
        const double fraction = absolute - static_cast<double>(whole_days);
        const CivilDate date = civil_from_days(whole_days);
        const int doy = static_cast<int>(whole_days - days_from_civil(date.year, 1, 1)) + 1;
        return static_cast<double>(date.year) +
               (static_cast<double>(doy) - 1.0 + fraction) / static_cast<double>(days_in_year(date.year));
    }

    if (contains(unit, "month")) {
        const double months = value;
        const double index = static_cast<double>(epoch.year * 12 + (epoch.month - 1)) + months;
        const long long whole = static_cast<long long>(std::floor(index));
        const double fraction = index - static_cast<double>(whole);
        const long long year = whole >= 0 ? whole / 12 : (whole - 11) / 12;
        const int month = static_cast<int>(whole - year * 12) + 1;
        return static_cast<double>(year) + (static_cast<double>(month) - 1.0 + fraction) / 12.0;
    }

    if (contains(unit, "year") || contains(unit, "annual"))
        return static_cast<double>(epoch.year) + value;

    throw std::runtime_error("unsupported time unit: " + unit);
}

NetCDFExportResult export_netcdf_to_csv(const std::string& path, const std::string& variable_name,
                                        const std::string& mode, const std::string& output) {
    NetCDFExportResult result;
    result.output_path = output;
    try {
        const NetCDFFile file = NetCDFFile::open(path);
        const NcVariable* variable = file.find_variable(variable_name);
        if (variable == nullptr) {
            std::ostringstream message;
            message << "variable '" << variable_name << "' not found in " << path
                    << "; available variables:";
            for (const auto& candidate : file.variables()) message << " " << candidate.name;
            throw std::runtime_error(message.str());
        }

        std::vector<std::string> dimension_names;
        for (int id : variable->dimension_ids)
            dimension_names.push_back(file.dimensions()[static_cast<std::size_t>(id)].name);

        const double scale = file.attribute_number(*variable, "scale_factor", 1.0);
        const double offset = file.attribute_number(*variable, "add_offset", 0.0);
        bool has_fill = false;
        const double fill = numeric_attribute(*variable, "_FillValue",
                                              std::numeric_limits<double>::quiet_NaN(), &has_fill);

        const std::size_t count = variable->element_count;
        const std::vector<double> raw = file.read_variable(variable_name);
        if (raw.size() != count)
            throw std::runtime_error("variable element count does not match its shape");

        const auto cleaned = [&](double value) -> double {
            if (!std::isfinite(value)) return std::numeric_limits<double>::quiet_NaN();
            if (has_fill && value == fill) return std::numeric_limits<double>::quiet_NaN();
            return value * scale + offset;
        };

        // Locate the time dimension and its coordinate variable.
        int time_index = -1;
        for (std::size_t i = 0; i < dimension_names.size(); ++i) {
            if (contains(lower_copy(dimension_names[i]), "time")) {
                time_index = static_cast<int>(i);
                break;
            }
        }
        if (time_index < 0 && variable->shape.size() > 1) time_index = 0;
        if (time_index < 0) time_index = 0;

        int lat_index = -1;
        int lon_index = -1;
        for (std::size_t i = 0; i < dimension_names.size(); ++i) {
            const std::string name = lower_copy(dimension_names[i]);
            if (contains(name, "lat")) lat_index = static_cast<int>(i);
            if (contains(name, "lon")) lon_index = static_cast<int>(i);
        }

        const NcVariable* time_variable = nullptr;
        std::vector<std::string> candidates = {dimension_names[static_cast<std::size_t>(time_index)],
                                               "time", "Time", "TIME"};
        for (const auto& candidate : candidates) {
            const NcVariable* found = file.find_variable(candidate);
            if (found != nullptr && found->shape.size() == 1 &&
                found->shape[0] == variable->shape[static_cast<std::size_t>(time_index)]) {
                time_variable = found;
                break;
            }
        }

        std::vector<double> time_values;
        std::string time_units;
        if (time_variable != nullptr) {
            time_values = file.read_variable(time_variable->name);
            time_units = file.attribute_text(*time_variable, "units");
            if (time_units.empty()) time_units = "days since 1970-01-01";
        }

        std::vector<double> lat_values;
        std::vector<double> lon_values;
        if (lat_index >= 0) {
            const NcVariable* lat_variable = file.find_variable(dimension_names[static_cast<std::size_t>(lat_index)]);
            if (lat_variable == nullptr) lat_variable = file.find_variable("lat");
            if (lat_variable != nullptr) lat_values = file.read_variable(lat_variable->name);
        }
        if (lon_index >= 0) {
            const NcVariable* lon_variable = file.find_variable(dimension_names[static_cast<std::size_t>(lon_index)]);
            if (lon_variable == nullptr) lon_variable = file.find_variable("lon");
            if (lon_variable != nullptr) lon_values = file.read_variable(lon_variable->name);
        }

        std::vector<std::size_t> strides(variable->shape.size(), 1);
        for (std::size_t i = variable->shape.size(); i-- > 1;) {
            strides[i - 1] = strides[i] * variable->shape[i];
        }

        std::ostringstream csv;
        csv.precision(12);
        json::Value summary = json::Value::object();
        summary.set("source", json::Value::from_string(path));
        summary.set("variable", json::Value::from_string(variable_name));
        summary.set("type", json::Value::from_string(netcdf_type_name(variable->type)));
        summary.set("mode", json::Value::from_string(mode));
        summary.set("records", json::Value::from_number(static_cast<double>(file.record_count())));

        json::Value dims = json::Value::array();
        for (std::size_t i = 0; i < dimension_names.size(); ++i) {
            json::Value item = json::Value::object();
            item.set("name", json::Value::from_string(dimension_names[i]));
            item.set("length", json::Value::from_number(static_cast<double>(variable->shape[i])));
            dims.push(item);
        }
        summary.set("dimensions", dims);

        if (mode == "timeseries") {
            if (variable->shape.empty())
                throw std::runtime_error("timeseries export needs at least one dimension");
            const std::size_t time_count = variable->shape[static_cast<std::size_t>(time_index)];
            if (!time_values.empty() && time_values.size() != time_count)
                throw std::runtime_error("time coordinate length does not match the data");

            csv << "time,value\n";
            std::size_t written = 0;
            for (std::size_t t = 0; t < time_count; ++t) {
                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                for (std::size_t flat = 0; flat < count; ++flat) {
                    const std::size_t index_in_time =
                        (flat / strides[static_cast<std::size_t>(time_index)]) %
                        variable->shape[static_cast<std::size_t>(time_index)];
                    if (index_in_time != t) continue;

                    double weight = 1.0;
                    if (lat_index >= 0) {
                        const std::size_t lat_position =
                            (flat / strides[static_cast<std::size_t>(lat_index)]) %
                            variable->shape[static_cast<std::size_t>(lat_index)];
                        if (lat_position < lat_values.size())
                            weight = std::cos(deg_to_rad(lat_values[lat_position]));
                        else
                            weight = 1.0;
                    }
                    const double value = cleaned(raw[flat]);
                    if (!std::isfinite(value)) continue;
                    weighted_sum += value * weight;
                    weight_sum += weight;
                }
                if (weight_sum <= 0.0) continue;
                const double time = time_values.empty()
                                        ? static_cast<double>(t)
                                        : decimal_year_from_time(time_values[t], time_units);
                csv << time << "," << (weighted_sum / weight_sum) << "\n";
                ++written;
            }
            summary.set("points", json::Value::from_number(static_cast<double>(written)));
            summary.set("time_units", json::Value::from_string(time_units));
        } else if (mode == "grid") {
            if (lat_index < 0 || lon_index < 0)
                throw std::runtime_error("grid export requires latitude and longitude dimensions");
            if (lat_values.empty() || lon_values.empty())
                throw std::runtime_error("grid export requires latitude and longitude coordinates");

            const std::size_t time_count = variable->shape[static_cast<std::size_t>(time_index)];
            const std::size_t last_time = time_count == 0 ? 0 : time_count - 1;
            const std::size_t lat_count = variable->shape[static_cast<std::size_t>(lat_index)];
            const std::size_t lon_count = variable->shape[static_cast<std::size_t>(lon_index)];
            if (lat_values.size() != lat_count || lon_values.size() != lon_count)
                throw std::runtime_error("coordinate lengths do not match the data grid");

            const auto edges = [](const std::vector<double>& coordinates) {
                std::vector<double> bounds(coordinates.size() + 1, 0.0);
                if (coordinates.empty()) return bounds;
                if (coordinates.size() == 1) {
                    bounds[0] = coordinates[0] - 0.5;
                    bounds[1] = coordinates[0] + 0.5;
                    return bounds;
                }
                for (std::size_t i = 0; i + 1 < coordinates.size(); ++i) {
                    const double mid = 0.5 * (coordinates[i] + coordinates[i + 1]);
                    if (i == 0) bounds[0] = coordinates[0] - (mid - coordinates[0]);
                    bounds[i + 1] = mid;
                }
                bounds[coordinates.size()] =
                    coordinates.back() +
                    (coordinates.back() - bounds[coordinates.size() - 1]);
                return bounds;
            };
            const std::vector<double> lat_bounds = edges(lat_values);
            const std::vector<double> lon_bounds = edges(lon_values);

            csv << "lat_min,lat_max,lon_min,lon_max,concentration\n";
            std::size_t written = 0;
            for (std::size_t i = 0; i < lat_count; ++i) {
                for (std::size_t j = 0; j < lon_count; ++j) {
                    const std::size_t flat = last_time * strides[static_cast<std::size_t>(time_index)] +
                                             i * strides[static_cast<std::size_t>(lat_index)] +
                                             j * strides[static_cast<std::size_t>(lon_index)];
                    if (flat >= count) continue;
                    const double value = cleaned(raw[flat]);
                    if (!std::isfinite(value)) continue;
                    csv << lat_bounds[i] << "," << lat_bounds[i + 1] << "," << lon_bounds[j] << ","
                        << lon_bounds[j + 1] << "," << value << "\n";
                    ++written;
                }
            }
            summary.set("cells", json::Value::from_number(static_cast<double>(written)));
            summary.set("time_index", json::Value::from_number(static_cast<double>(last_time)));
        } else {
            throw std::runtime_error("unknown export mode '" + mode + "'; use timeseries or grid");
        }

        write_text_file(output, csv.str());
        summary.set("output", json::Value::from_string(output));
        result.summary_json = json::dump(summary, 2);
        result.ok = true;
        return result;
    } catch (const std::exception& exception) {
        result.error = exception.what();
        return result;
    }
}

}  // namespace climate
