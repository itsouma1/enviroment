#include "openfetch.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "fetch.hpp"
#include "json.hpp"

namespace climate {
namespace {

std::string trimmed(const std::string& text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin]))) ++begin;
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) --end;
    return text.substr(begin, end - begin);
}

std::string lowered(std::string text) {
    for (char& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return text;
}

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> cells;
    std::string current;
    for (char c : line) {
        if (c == ',') {
            cells.push_back(trimmed(current));
            current.clear();
        } else {
            current.push_back(c);
        }
    }
    cells.push_back(trimmed(current));
    return cells;
}

bool numeric_cell(const std::string& raw, double* value) {
    std::string text = trimmed(raw);
    std::replace(text.begin(), text.end(), '*', ' ');
    text = trimmed(text);
    if (text.empty()) return false;
    if (text == "-999" || text == "-999.0" || text == "-9999" || lowered(text) == "nan") return false;
    char* end = nullptr;
    const double parsed = std::strtod(text.c_str(), &end);
    if (end == text.c_str() || !std::isfinite(parsed)) return false;
    *value = parsed;
    return true;
}

double decimal_year(int year, int month) {
    const double value = static_cast<double>(year) + (static_cast<double>(month) - 0.5) / 12.0;
    return std::round(value * 10000.0) / 10000.0;
}

std::string utc_now() {
    const std::time_t now = std::time(nullptr);
    std::tm utc{};
#if defined(_WIN32)
    gmtime_s(&utc, &now);
#else
    gmtime_r(&now, &utc);
#endif
    char buffer[40];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S+00:00", &utc);
    return std::string(buffer);
}

// One HTTPS GET. curl is present on Linux, macOS and Windows 10+; PowerShell and
// Python are only fallbacks, so nothing has to be installed for the open sources.
bool http_get(const std::string& url, std::string* body, std::string* error) {
    const CommandResult curl = run_command("curl -fsSL --max-time 90 \"" + url + "\"");
    if (curl.exit_code == 0 && !curl.output.empty()) {
        *body = curl.output;
        return true;
    }
    const CommandResult powershell = run_command(
        "powershell -NoProfile -Command \"(Invoke-WebRequest -UseBasicParsing -Uri '" + url +
        "').Content\"");
    if (powershell.exit_code == 0 && !powershell.output.empty()) {
        *body = powershell.output;
        return true;
    }
    const CommandResult python = run_command(
        "python3 -c \"import urllib.request;print(urllib.request.urlopen('" + url +
        "',timeout=90).read().decode())\"");
    if (python.exit_code == 0 && !python.output.empty()) {
        *body = python.output;
        return true;
    }
    *error = "could not reach " + url +
             " — no HTTP client answered (tried curl, PowerShell, python3)";
    return false;
}

bool parse_gistemp(const std::string& text, OpenSeries* series, std::string* error) {
    std::istringstream stream(text);
    std::string line;
    bool header_seen = false;
    std::size_t rows = 0;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (trimmed(line).empty()) continue;
        const std::vector<std::string> cells = split_csv_line(line);
        if (!header_seen) {
            if (!cells.empty() && lowered(cells[0]) == "year") header_seen = true;
            continue;
        }
        if (cells.empty()) continue;
        const std::string year_text = trimmed(cells[0]);
        if (year_text.size() != 4 ||
            !std::all_of(year_text.begin(), year_text.end(), [](unsigned char c) {
                return std::isdigit(c) != 0;
            }))
            continue;
        const int year = std::atoi(year_text.c_str());
        for (int month = 1; month <= 12; ++month) {
            if (static_cast<std::size_t>(month) >= cells.size()) break;
            double value = 0.0;
            if (!numeric_cell(cells[static_cast<std::size_t>(month)], &value)) continue;
            series->time.push_back(decimal_year(year, month));
            series->value.push_back(value);
            ++rows;
        }
    }
    if (rows < 12) {
        *error = "the GISTEMP table had no usable monthly values";
        return false;
    }
    return true;
}

bool parse_power(const std::string& text, const std::string& parameter, OpenSeries* series,
                 std::string* error) {
    std::istringstream stream(text);
    std::string line;
    std::vector<std::string> header;
    const char* kMonths[] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                             "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
    std::size_t rows = 0;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string stripped = trimmed(line);
        if (stripped.empty() || stripped.substr(0, 12) == "-END HEADER-") continue;
        const std::vector<std::string> cells = split_csv_line(stripped);
        if (header.empty()) {
            if (std::find(cells.begin(), cells.end(), "YEAR") != cells.end()) header = cells;
            continue;
        }
        if (cells.size() != header.size()) continue;
        std::size_t year_index = header.size();
        std::size_t parameter_index = header.size();
        for (std::size_t i = 0; i < header.size(); ++i) {
            if (header[i] == "YEAR") year_index = i;
            if (header[i] == "PARAMETER") parameter_index = i;
        }
        if (year_index >= cells.size()) continue;
        if (parameter_index < cells.size() && cells[parameter_index] != parameter) continue;
        double year_value = 0.0;
        if (!numeric_cell(cells[year_index], &year_value)) continue;
        const int year = static_cast<int>(year_value);
        for (int month = 1; month <= 12; ++month) {
            std::size_t column = header.size();
            for (std::size_t i = 0; i < header.size(); ++i)
                if (header[i] == kMonths[month - 1]) column = i;
            if (column >= cells.size()) continue;
            double value = 0.0;
            if (!numeric_cell(cells[column], &value)) continue;
            series->time.push_back(decimal_year(year, month));
            series->value.push_back(value);
            ++rows;
        }
    }
    if (rows < 12) {
        *error = "the POWER response had no usable monthly values for " + parameter;
        return false;
    }
    return true;
}

bool parse_seaice(const std::vector<std::pair<int, std::string>>& files, OpenSeries* series,
                  std::string* error) {
    std::vector<std::pair<int, int>> keys;
    std::vector<double> values;
    std::string last_url;
    for (const auto& entry : files) {
        std::string body;
        std::string http_error;
        if (!http_get(entry.second, &body, &http_error)) {
            *error = http_error;
            return false;
        }
        last_url = entry.second;
        std::istringstream stream(body);
        std::string line;
        while (std::getline(stream, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            const std::vector<std::string> cells = split_csv_line(line);
            if (cells.size() < 5) continue;
            const std::string year_text = trimmed(cells[0]);
            if (year_text.size() != 4 ||
                !std::all_of(year_text.begin(), year_text.end(), [](unsigned char c) {
                    return std::isdigit(c) != 0;
                }))
                continue;
            double month_value = 0.0;
            double extent = 0.0;
            if (!numeric_cell(cells[1], &month_value)) continue;
            if (!numeric_cell(cells[4], &extent)) continue;
            keys.emplace_back(std::atoi(year_text.c_str()), static_cast<int>(month_value));
            values.push_back(extent);
        }
    }
    if (keys.size() < 12) {
        *error = "the Sea Ice Index files had no usable monthly extent";
        return false;
    }
    std::vector<std::size_t> order(keys.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });
    for (std::size_t index : order) {
        series->time.push_back(decimal_year(keys[index].first, keys[index].second));
        series->value.push_back(values[index]);
    }
    series->url = last_url;
    return true;
}

}  // namespace

bool fetch_open_series(const std::string& source, double latitude, double longitude,
                       const std::string& start, const std::string& end,
                       const std::string& region, OpenSeries* series, std::string* error) {
    if (series == nullptr) return false;
    OpenSeries out;
    out.source = source;
    out.latitude = latitude;
    out.longitude = longitude;
    std::string body;

    if (source == "gistemp") {
        out.key = "live_gistemp";
        out.label = "GISTEMP v4 global mean anomaly";
        out.provider = "NASA Goddard Institute for Space Studies";
        out.variable = "temperature anomaly";
        out.units = "degC vs 1951-1980";
        out.url = "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv";
        if (!http_get(out.url, &body, error)) return false;
        if (!parse_gistemp(body, &out, error)) return false;
    } else if (source == "power") {
        const int start_year =
            start.size() >= 4 ? std::atoi(start.substr(0, 4).c_str()) : 1981;
        const int end_year = end.size() >= 4 ? std::atoi(end.substr(0, 4).c_str()) : 2024;
        std::ostringstream url;
        url << "https://power.larc.nasa.gov/api/temporal/monthly/point"
            << "?parameters=T2M&community=RE"
            << "&longitude=" << longitude << "&latitude=" << latitude
            << "&start=" << start_year << "&end=" << end_year << "&format=CSV";
        out.url = url.str();
        std::ostringstream key;
        key << "live_power_t2m_" << (latitude < 0 ? "m" : "")
            << std::llround(std::fabs(latitude)) << "_" << (longitude < 0 ? "m" : "")
            << std::llround(std::fabs(longitude));
        out.key = key.str();
        out.label = "NASA POWER surface meteorology";
        out.provider = "NASA Langley Research Center POWER project";
        out.variable = "T2M";
        out.units = "degC";
        if (!http_get(out.url, &body, error)) return false;
        if (!parse_power(body, "T2M", &out, error)) return false;
    } else if (source == "seaice") {
        const bool north = lowered(region) != "s";
        const std::string hemisphere = north ? "north" : "south";
        const std::string prefix = north ? "N" : "S";
        const std::string base =
            "https://noaadata.apps.nsidc.org/NOAA/G02135/" + hemisphere + "/monthly/data/";
        std::vector<std::pair<int, std::string>> files;
        for (int month = 1; month <= 12; ++month) {
            std::ostringstream name;
            name << prefix << "_" << std::setw(2) << std::setfill('0') << month
                 << "_extent_v4.0.csv";
            files.emplace_back(month, base + name.str());
        }
        out.key = std::string("live_seaice_") + (north ? "n" : "s");
        out.label = "Sea Ice Index monthly extent";
        out.provider = "NASA NSIDC DAAC (G02135)";
        out.variable = std::string("extent_") + (north ? "N" : "S");
        out.units = "million km^2";
        if (!parse_seaice(files, &out, error)) return false;
    } else {
        *error = "unknown open source '" + source +
                 "' — this path covers gistemp, power and seaice";
        return false;
    }

    if (out.time.size() < 12) {
        *error = "the response contained fewer than twelve monthly values";
        return false;
    }
    out.retrieved_at = utc_now();
    *series = out;
    return true;
}

bool write_open_series(const OpenSeries& series, const std::string& directory,
                       std::string* csv_path, std::string* error) {
    std::error_code ignored;
    std::filesystem::create_directories(directory, ignored);
    const std::string csv = directory + "/" + series.key + ".csv";
    std::ofstream output(csv, std::ios::binary);
    if (!output) {
        if (error != nullptr) *error = "cannot write " + csv;
        return false;
    }
    output << "time,value\n";
    output.setf(std::ios::fixed);
    output << std::setprecision(4);
    for (std::size_t i = 0; i < series.time.size(); ++i) {
        output << series.time[i] << ',';
        output << std::setprecision(6) << series.value[i] << '\n';
        output << std::setprecision(4);
    }
    output.close();

    json::Value root = json::Value::object();
    root.set("url", json::Value::from_string(series.url));
    root.set("variable", json::Value::from_string(series.variable));
    root.set("units", json::Value::from_string(series.units));
    root.set("source", json::Value::from_string(series.source));
    root.set("key", json::Value::from_string(series.key));
    root.set("file", json::Value::from_string(csv));
    root.set("samples", json::Value::from_number(static_cast<double>(series.time.size())));
    root.set("first_time", json::Value::from_number(series.time.front()));
    root.set("last_time", json::Value::from_number(series.time.back()));
    root.set("first_value", json::Value::from_number(series.value.front()));
    root.set("last_value", json::Value::from_number(series.value.back()));
    root.set("retrieved_at", json::Value::from_string(
                                  series.retrieved_at.empty() ? utc_now() : series.retrieved_at));
    if (series.source == "power") {
        root.set("latitude", json::Value::from_number(series.latitude));
        root.set("longitude", json::Value::from_number(series.longitude));
    }
    std::ofstream sidecar(directory + "/" + series.key + ".json", std::ios::binary);
    if (sidecar) sidecar << json::dump(root, 2) << "\n";
    if (csv_path != nullptr) *csv_path = csv;
    return true;
}

}  // namespace climate
