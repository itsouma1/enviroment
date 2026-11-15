#pragma once

#include <string>
#include <utility>
#include <vector>

namespace climate {

// One fetched monthly series, ready to be written as time,value.
struct OpenSeries {
    std::string source;    // gistemp | power | seaice
    std::string key;       // file stem, e.g. live_gistemp
    std::string label;     // human description
    std::string provider;  // who serves it
    std::string variable;  // variable name
    std::string units;
    std::string url;       // the URL actually requested
    std::string retrieved_at;  // UTC ISO 8601 timestamp of the fetch
    double latitude = 0.0;
    double longitude = 0.0;
    std::vector<double> time;   // decimal years
    std::vector<double> value;
};

// Fetches one of the open NASA series over HTTPS. None of these need an account:
// GISTEMP v4 (GISS), POWER monthly point (LaRC) and the NSIDC Sea Ice Index.
//
// HTTP goes through curl when it is present, which it is on Linux and macOS and
// on Windows 10 and later; PowerShell and Python are tried as fallbacks. No
// Python package, credentials or extra install is required for these three
// sources.
bool fetch_open_series(const std::string& source, double latitude, double longitude,
                       const std::string& start, const std::string& end,
                       const std::string& region, const std::string& parameter,
                       OpenSeries* series, std::string* error);

// Signs in to Earthdata Login with a username and password. The credentials are
// read from the child environment through variable references, so they never
// appear on a command line, in a file or in a log.
bool verify_earthdata_password(const std::vector<std::pair<std::string, std::string>>& env,
                               std::vector<std::string>* usernames, std::string* detail);

// Checks an Earthdata Login bearer token against the EDL profile endpoint. No
// Python and no password are involved: the token is sent as a bearer header.
// Returns true when the token was accepted; `detail` always explains what
// happened, including a rejected or unreachable check.
bool verify_earthdata_token(const std::string& token, std::string* username,
                            std::string* detail);

// Writes time,value CSV plus a JSON sidecar next to it, the same schema the
// Python helper produces. Returns the CSV path.
bool write_open_series(const OpenSeries& series, const std::string& directory,
                       std::string* csv_path, std::string* error);

}  // namespace climate
