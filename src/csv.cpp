#include "csv.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace climate {
namespace {

std::vector<std::string> split_line(const std::string& line, char delimiter) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (quoted) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    field.push_back('"');
                    ++i;
                } else {
                    quoted = false;
                }
            } else {
                field.push_back(c);
            }
        } else if (c == '"' && trim_copy(field).empty()) {
            field.clear();
            quoted = true;
        } else if (c == delimiter) {
            fields.push_back(field);
            field.clear();
        } else {
            field.push_back(c);
        }
    }
    fields.push_back(field);
    return fields;
}

bool is_missing_marker(const std::string& text) {
    std::string lowered;
    lowered.reserve(text.size());
    for (char c : text) lowered.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    return lowered.empty() || lowered == "na" || lowered == "n/a" || lowered == "nan" ||
           lowered == "null" || lowered == "-9999" || lowered == "-999" || lowered == "none";
}

}  // namespace

std::string trim_copy(const std::string& text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin]))) ++begin;
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) --end;
    return text.substr(begin, end - begin);
}

int Table::column(const std::string& name) const {
    const std::string target = trim_copy(name);
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (trim_copy(header[i]) == target) return static_cast<int>(i);
    }
    for (std::size_t i = 0; i < header.size(); ++i) {
        std::string a = trim_copy(header[i]);
        std::string b = target;
        for (auto& c : a) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        for (auto& c : b) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (a == b) return static_cast<int>(i);
    }
    return -1;
}

std::string Table::cell(std::size_t row, int column_index) const {
    if (row >= rows.size()) return std::string();
    if (column_index < 0 || static_cast<std::size_t>(column_index) >= rows[row].size())
        return std::string();
    return rows[row][static_cast<std::size_t>(column_index)];
}

Table read_csv(const std::string& path, char delimiter) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open file: " + path);

    Table table;
    std::string line;
    bool have_header = false;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (trim_copy(line).empty()) continue;
        std::vector<std::string> fields = split_line(line, delimiter);
        if (!have_header) {
            for (auto& field : fields) table.header.push_back(trim_copy(field));
            have_header = true;
            continue;
        }
        for (auto& field : fields) field = trim_copy(field);
        if (fields.size() < table.header.size()) fields.resize(table.header.size());
        if (fields.size() > table.header.size()) fields.resize(table.header.size());
        table.rows.push_back(std::move(fields));
    }
    if (!have_header) throw std::runtime_error("file has no header row: " + path);
    return table;
}

bool parse_number(const std::string& raw, double& out) {
    const std::string text = trim_copy(raw);
    if (is_missing_marker(text)) return false;
    char* end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (end == nullptr) return false;
    while (end != nullptr && *end != '\0' && std::isspace(static_cast<unsigned char>(*end))) ++end;
    if (end == nullptr || *end != '\0') return false;
    out = value;
    return true;
}

std::vector<double> numeric_column(const Table& table, int column_index) {
    std::vector<double> values;
    values.reserve(table.rows.size());
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        double value = 0.0;
        if (parse_number(table.cell(i, column_index), value))
            values.push_back(value);
        else
            values.push_back(std::nan(""));
    }
    return values;
}

void write_text_file(const std::string& path, const std::string& content) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("cannot write file: " + path);
    output.write(content.data(), static_cast<std::streamsize>(content.size()));
    if (!output) throw std::runtime_error("failed while writing file: " + path);
}

}  // namespace climate
