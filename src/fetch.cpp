#include "fetch.hpp"

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#if !defined(_WIN32)
#include <sys/wait.h>
#endif

namespace climate {

namespace {
const std::vector<RemoteDataset> kDatasets = {
    {"sst",
     "sst",
     "MODISA_L3m_SST",
     "",
     "PO.DAAC",
     "sst,sst4",
     "degree_C",
     "timeseries",
     "MODIS Aqua level-3 mapped SST (4 km, 8-day, nighttime thermal IR)."},
    {"ndvi",
     "ndvi",
     "MOD13A3",
     "",
     "LP DAAC",
     "NDVI,ndvi",
     "",
     "timeseries",
     "MODIS/Terra vegetation indices, monthly L3 global 1 km."},
    {"co2",
     "co2",
     "OCO2_L2_Lite_FP",
     "",
     "GES DISC",
     "xco2,XCO2",
     "ppm",
     "timeseries",
     "OCO-2 level-2 lite geolocated XCO2 retrievals."},
    {"sea_ice_extent",
     "sea_ice_extent",
     "NSIDC-0051",
     "",
     "NSIDC DAAC",
     "seaice_conc,cdr_seaice_conc,seaice_conc_monthly",
     "fraction",
     "extent",
     "Sea ice concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS (CDR)."},
    {"sea_ice_concentration",
     "sea_ice_concentration",
     "NSIDC-0051",
     "",
     "NSIDC DAAC",
     "seaice_conc,cdr_seaice_conc,seaice_conc_monthly",
     "fraction",
     "timeseries",
     "Area-weighted mean sea ice concentration from the NSIDC CDR."},
    {"precipitation",
     "precipitation",
     "GPCPMON",
     "",
     "GES DISC",
     "precip,precipitation",
     "mm/day",
     "timeseries",
     "GPCP monthly precipitation analysis, 0.5 degree."},
    {"imerg",
     "precipitation",
     "GPM_3IMERGM",
     "",
     "GES DISC",
     "precipitation",
     "mm/hr",
     "timeseries",
     "IMERG monthly precipitation, 0.1 degree."},
    {"air_temperature",
     "air_temperature",
     "GLDAS_NOAH025_M",
     "",
     "GES DISC",
     "Tair_f_inst,tair",
     "K",
     "timeseries",
     "GLDAS Noah monthly land surface air temperature, 0.25 degree."},
    {"ozone",
     "ozone",
     "OMDOAO3e",
     "",
     "GES DISC",
     "ColumnAmountO3,ozone",
     "DU",
     "timeseries",
     "OMI total column ozone, monthly level-3 global grid."},
    {"aerosol",
     "aerosol",
     "OMAERO",
     "",
     "GES DISC",
     "AerosolIndex,aerosol_index",
     "unitless",
     "timeseries",
     "OMI UV aerosol index, level-2 swath retrievals."},
    {"soil_moisture",
     "soil_moisture",
     "SPL3SMP_E",
     "",
     "NSIDC DAAC",
     "soil_moisture,sm",
     "m3/m3",
     "timeseries",
     "SMAP L3 enhanced soil moisture, 9 km."},
    {"snow_cover",
     "snow_cover",
     "MOD10CM",
     "",
     "NSIDC DAAC",
     "Snow_Cover,snow_cover",
     "percent",
     "timeseries",
     "MODIS/Terra monthly snow cover extent, CMG 0.05 degree."},
    {"sea_level",
     "sea_level",
     "SWOT_L2_LR_SSH_D",
     "",
     "PO.DAAC",
     "ssha,sea_surface_height_anomaly",
     "m",
     "timeseries",
     "Multi-mission altimetry sea surface height anomaly, level-4 gridded."},
    {"ocean_color",
     "ocean_color",
     "MODISA_L3m_CHL",
     "",
     "OB.DAAC",
     "chlor_a,chlorophyll",
     "mg/m3",
     "timeseries",
     "MODIS Aqua level-3 mapped chlorophyll-a, 4 km monthly."},
};
}  // namespace

const std::vector<RemoteDataset>& remote_datasets() { return kDatasets; }

const RemoteDataset* find_remote_dataset(const std::string& key) {
    for (const auto& dataset : kDatasets)
        if (dataset.key == key) return &dataset;
    return nullptr;
}

std::string shell_quote(const std::string& text) {
    std::string quoted = "'";
    for (char c : text) {
        if (c == '\'')
            quoted += "'\\''";
        else
            quoted.push_back(c);
    }
    quoted.push_back('\'');
    return quoted;
}

namespace {

bool file_exists(const std::string& path) {
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr) return false;
    std::fclose(file);
    return true;
}

}  // namespace

std::string resolve_python_interpreter() {
    const char* configured = std::getenv("CLIMATE_PYTHON");
    if (configured != nullptr && *configured != '\0') return configured;
    if (file_exists(".venv/bin/python")) return ".venv/bin/python";
    return "python3";
}

std::string resolve_fetch_wrapper() {
    const char* configured = std::getenv("CLIMATE_WRAPPER");
    if (configured != nullptr && *configured != '\0') return configured;
    return "tools/nasa_fetch.py";
}

std::string build_command(const std::vector<std::string>& arguments) {
    std::ostringstream command;
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        if (i > 0) command << ' ';
        command << shell_quote(arguments[i]);
    }
    return command.str();
}

CommandResult run_command_with_env(const std::string& command,
                                   const std::vector<std::pair<std::string, std::string>>& env) {
    CommandResult result;
#if defined(_WIN32)
    std::string prefixed;
    for (const auto& item : env) prefixed += "set " + item.first + "=" + item.second + "&& ";
    return run_command(prefixed + command);
#else
    int fds[2];
    if (pipe(fds) != 0) {
        result.exit_code = -1;
        result.output = "failed to create a pipe";
        return result;
    }
    const pid_t pid = fork();
    if (pid == 0) {
        close(fds[0]);
        dup2(fds[1], STDOUT_FILENO);
        dup2(fds[1], STDERR_FILENO);
        close(fds[1]);
        for (const auto& item : env) setenv(item.first.c_str(), item.second.c_str(), 1);
        execl("/bin/sh", "sh", "-c", command.c_str(), static_cast<char*>(nullptr));
        _exit(127);
    }
    if (pid < 0) {
        close(fds[0]);
        close(fds[1]);
        result.exit_code = -1;
        result.output = "failed to start command: " + command;
        return result;
    }
    close(fds[1]);
    char buffer[4096];
    ssize_t got = 0;
    while ((got = read(fds[0], buffer, sizeof(buffer))) > 0)
        result.output.append(buffer, static_cast<std::size_t>(got));
    close(fds[0]);
    int status = 0;
    waitpid(pid, &status, 0);
    result.exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : status;
    return result;
#endif
}

CommandResult run_command(const std::string& command) {
    CommandResult result;
    const std::string full = command + " 2>&1";
    std::FILE* pipe = popen(full.c_str(), "r");
    if (pipe == nullptr) {
        result.exit_code = -1;
        result.output = "failed to start command: " + command;
        return result;
    }
    char buffer[4096];
    while (std::fgets(buffer, sizeof(buffer), pipe) != nullptr) result.output += buffer;
    const int status = pclose(pipe);
#if defined(_WIN32)
    result.exit_code = status;
#else
    if (status == -1)
        result.exit_code = -1;
    else if (WIFEXITED(status))
        result.exit_code = WEXITSTATUS(status);
    else
        result.exit_code = status;
#endif
    return result;
}

}  // namespace climate
