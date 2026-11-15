#ifndef CLIMATE_MATURITY_HPP
#define CLIMATE_MATURITY_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace climate {

struct MaturityScore {
    int scientific = 1;
    int preservation = 1;
    int societal = 1;
};

bool maturity_accept(const MaturityScore& score, int scientific_required,
                     int preservation_required, int societal_required);

int maturity_imbalance(const MaturityScore& score);

// ---- GHRSST CDAF ---------------------------------------------------------

struct CdafMetrics {
    double systematic_bias_buoy = 0.0;
    double systematic_bias_argo = 0.0;
    double geographic_bias_buoy = 0.0;
    double geographic_bias_argo = 0.0;
    double non_systematic_dispersion = 0.0;
    double stability_trend = 0.0;
    double sensitivity = 0.0;
    double record_years = 0.0;
    bool gds2_format = false;
};

struct CdafRequirements {
    double bias = 0.0;
    double geographic = 0.0;
    double dispersion = 0.0;
    double stability = 0.0;
    double sensitivity = 0.0;
};

bool cdaf_eligible(bool gds2_format, double record_years, double minimum_years);

bool cdaf_accept(const CdafMetrics& metrics, const CdafRequirements& requirements);

double geographic_bias_variation(const std::vector<double>& satellite,
                                 const std::vector<double>& reference,
                                 const std::vector<double>& latitudes,
                                 const std::vector<double>& longitudes, double latitude_bin,
                                 double longitude_bin);

double robust_dispersion(const std::vector<double>& differences,
                         const std::vector<double>& geographic_bias);

double stability_trend(const std::vector<double>& differences);

// ---- CEOS-ARD STAC validation -------------------------------------------

bool ceos_ard_valid(const std::string& type, const std::string& specification,
                    const std::string& version);

}  // namespace climate

#endif
