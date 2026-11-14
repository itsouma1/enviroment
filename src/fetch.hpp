#pragma once

#include <string>
#include <vector>

namespace climate {

// A NASA Earthdata product that can be downloaded and converted by the Python
// wrapper. C++ never parses HDF5/NetCDF-4 itself; the wrapper converts the
// downloaded granules to the CSV schema this tool analyses.
struct RemoteDataset {
    std::string key;          // CLI key
    std::string catalog_id;   // matching catalog variable id
    std::string short_name;   // CMR short name
    std::string version;      // CMR version
    std::string provider;     // DAAC
    std::string variable;     // preferred science variable
    std::string units;
    std::string mode;         // timeseries | extent | grid
    std::string description;
};

const std::vector<RemoteDataset>& remote_datasets();

const RemoteDataset* find_remote_dataset(const std::string& key);

// Quotes a value for safe inclusion in a POSIX shell command.
std::string shell_quote(const std::string& text);

// Resolves the Python interpreter: $CLIMATE_PYTHON, then .venv/bin/python,
// then python3.
std::string resolve_python_interpreter();

// Resolves the wrapper script: $CLIMATE_WRAPPER, then tools/nasa_fetch.py.
std::string resolve_fetch_wrapper();

struct CommandResult {
    int exit_code = 0;
    std::string output;
};

CommandResult run_command(const std::string& command);

// Builds a shell command from an argv-like vector, quoting every argument.
std::string build_command(const std::vector<std::string>& arguments);

}  // namespace climate
