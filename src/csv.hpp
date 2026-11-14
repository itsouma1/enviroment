#pragma once

#include <string>
#include <vector>

namespace climate {

struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;

    std::size_t size() const { return rows.size(); }
    bool empty() const { return rows.empty(); }

    // Returns -1 when the column is absent.
    int column(const std::string& name) const;
    std::string cell(std::size_t row, int column_index) const;
};

// Reads a delimited text file with a header row. Quoted fields, escaped quotes
// and CRLF line endings are supported. Throws std::runtime_error on failure.
Table read_csv(const std::string& path, char delimiter = ',');

// Parses a numeric field. Recognised missing-value markers ("", "NA", "N/A",
// "NaN", "null", "-9999") return false and leave the output untouched.
bool parse_number(const std::string& text, double& out);

// Returns the numeric values of one column, using NaN for missing entries.
std::vector<double> numeric_column(const Table& table, int column_index);

void write_text_file(const std::string& path, const std::string& content);

std::string trim_copy(const std::string& text);

}  // namespace climate
