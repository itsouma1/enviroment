#include "maturity.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double median(std::vector<double> values) {
    require(!values.empty(), "median: empty");
    const std::size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<long>(middle), values.end());
    if (values.size() % 2 == 1) return values[middle];
    const double upper = values[middle];
    std::nth_element(values.begin(), values.begin() + static_cast<long>(middle - 1), values.end());
    return 0.5 * (values[middle - 1] + upper);
}

bool valid_version(const std::string& version) {
    int dots = 0;
    bool previous_digit = false;
    for (char character : version) {
        if (character == '.') {
            if (!previous_digit || ++dots > 2) return false;
            previous_digit = false;
        } else if (character >= '0' && character <= '9') {
            previous_digit = true;
        } else {
            return false;
        }
    }
    return previous_digit;
}

}  // namespace

bool maturity_accept(const MaturityScore& score, int scientific_required,
                     int preservation_required, int societal_required) {
    return score.scientific >= scientific_required &&
           score.preservation >= preservation_required &&
           score.societal >= societal_required;
}

int maturity_imbalance(const MaturityScore& score) {
    const int maximum = std::max({score.scientific, score.preservation, score.societal});
    const int minimum = std::min({score.scientific, score.preservation, score.societal});
    return maximum - minimum;
}

bool cdaf_eligible(bool gds2_format, double record_years, double minimum_years) {
    return gds2_format && record_years >= minimum_years;
}

bool cdaf_accept(const CdafMetrics& metrics, const CdafRequirements& requirements) {
    return cdaf_eligible(metrics.gds2_format, metrics.record_years, 10.0) &&
           std::fabs(metrics.systematic_bias_buoy) <= requirements.bias &&
           std::fabs(metrics.systematic_bias_argo) <= requirements.bias &&
           metrics.geographic_bias_buoy <= requirements.geographic &&
           metrics.geographic_bias_argo <= requirements.geographic &&
           metrics.non_systematic_dispersion <= requirements.dispersion &&
           std::fabs(metrics.stability_trend) <= requirements.stability &&
           metrics.sensitivity <= requirements.sensitivity;
}

double geographic_bias_variation(const std::vector<double>& satellite,
                                 const std::vector<double>& reference,
                                 const std::vector<double>& latitudes,
                                 const std::vector<double>& longitudes, double latitude_bin,
                                 double longitude_bin) {
    require(satellite.size() == reference.size() && satellite.size() == latitudes.size() &&
                satellite.size() == longitudes.size() && !satellite.empty() &&
                latitude_bin > 0.0 && longitude_bin > 0.0,
            "geographic_bias_variation: inputs");
    std::map<std::pair<long long, long long>, std::vector<double>> bins;
    for (std::size_t i = 0; i < satellite.size(); ++i) {
        const long long lat_index =
            static_cast<long long>(std::floor(latitudes[i] / latitude_bin));
        const long long lon_index =
            static_cast<long long>(std::floor(longitudes[i] / longitude_bin));
        bins[{lat_index, lon_index}].push_back(satellite[i] - reference[i]);
    }
    std::vector<double> medians;
    for (auto& entry : bins) medians.push_back(median(entry.second));
    require(medians.size() >= 2, "geographic_bias_variation: too few bins");
    double mean = 0.0;
    for (double value : medians) mean += value;
    mean /= static_cast<double>(medians.size());
    double variance = 0.0;
    for (double value : medians) variance += (value - mean) * (value - mean);
    return std::sqrt(variance / static_cast<double>(medians.size() - 1));
}

double robust_dispersion(const std::vector<double>& differences,
                         const std::vector<double>& geographic_bias) {
    require(differences.size() == geographic_bias.size() && !differences.empty(),
            "robust_dispersion: size");
    std::vector<double> adjusted(differences.size(), 0.0);
    for (std::size_t i = 0; i < differences.size(); ++i)
        adjusted[i] = differences[i] - geographic_bias[i];
    std::vector<double> absolute(adjusted.size(), 0.0);
    for (std::size_t i = 0; i < adjusted.size(); ++i)
        absolute[i] = std::fabs(adjusted[i] - median(adjusted));
    return 1.4826 * median(absolute);
}

double stability_trend(const std::vector<double>& differences) {
    require(differences.size() >= 2, "stability_trend: size");
    const double n = static_cast<double>(differences.size());
    const double mean_t = 0.5 * (n - 1.0);
    double mean_y = 0.0;
    for (double value : differences) mean_y += value;
    mean_y /= n;
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t t = 0; t < differences.size(); ++t) {
        const double dt = static_cast<double>(t) - mean_t;
        numerator += dt * (differences[t] - mean_y);
        denominator += dt * dt;
    }
    return numerator / denominator;
}

bool ceos_ard_valid(const std::string& type, const std::string& specification,
                    const std::string& version) {
    if (!valid_version(version)) return false;
    if (type == "optical") {
        return specification == "SR" || specification == "ST" ||
               specification == "AR" || specification == "NLSR";
    }
    if (type == "radar") {
        return specification == "NRB" || specification == "POL" ||
               specification == "ORB" || specification == "GSLC";
    }
    return false;
}

}  // namespace climate
