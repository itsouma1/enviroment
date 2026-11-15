#include "validation.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double standardised_difference_variance(const std::vector<double>& observed,
                                        const std::vector<double>& reference,
                                        double observation_uncertainty,
                                        double reference_uncertainty,
                                        double representation_error) {
    require(!observed.empty() && observed.size() == reference.size(),
            "standardised_difference_variance: size mismatch");
    double observed_mean = 0.0, reference_mean = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        observed_mean += observed[i];
        reference_mean += reference[i];
    }
    observed_mean /= static_cast<double>(observed.size());
    reference_mean /= static_cast<double>(reference.size());
    const double variance = observation_uncertainty * observation_uncertainty +
                            reference_uncertainty * reference_uncertainty +
                            representation_error * representation_error;
    require(variance > 0.0, "standardised_difference_variance: zero variance");
    double total = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        const double difference = (observed[i] - observed_mean) - (reference[i] - reference_mean);
        total += difference * difference / variance;
    }
    return total / static_cast<double>(observed.size());
}

double quadruple_collocation_covariance(double covariance, double scale_i, double scale_j,
                                        double truth_variance) {
    require(std::fabs(scale_i * scale_j) > 0.0,
            "quadruple_collocation_covariance: zero scale");
    return covariance / (scale_i * scale_j) - truth_variance;
}

double linear_ice_fraction(double ice_length, double total_length) {
    require(total_length > 0.0, "linear_ice_fraction: zero total length");
    require(ice_length >= 0.0 && ice_length <= total_length, "linear_ice_fraction: length");
    return 100.0 * ice_length / total_length;
}

double sea_ice_quadratic_bias(double a0, double a1, double a2, double sea_ice_concentration) {
    return a0 + a1 * sea_ice_concentration + a2 * sea_ice_concentration * sea_ice_concentration;
}

double sea_ice_bias_derivative(double a1, double a2, double sea_ice_concentration) {
    return a1 + 2.0 * a2 * sea_ice_concentration;
}

double inverse_variance_merge(const std::vector<double>& values,
                              const std::vector<double>& variances) {
    require(!values.empty() && values.size() == variances.size(),
            "inverse_variance_merge: size mismatch");
    double weighted = 0.0, weights = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        require(variances[i] > 0.0, "inverse_variance_merge: non-positive variance");
        const double weight = 1.0 / variances[i];
        weighted += weight * values[i];
        weights += weight;
    }
    return weighted / weights;
}

double precipitation_bias_sum(double hit_bias, double missed_bias, double false_bias) {
    return hit_bias + missed_bias + false_bias;
}

Matrix observation_error_covariance(const Matrix& retrieval, const Matrix& representation) {
    require(retrieval.size() == representation.size(), "observation_error_covariance: size");
    Matrix result = retrieval;
    for (std::size_t i = 0; i < retrieval.size(); ++i) {
        require(retrieval[i].size() == representation[i].size(),
                "observation_error_covariance: ragged");
        for (std::size_t j = 0; j < retrieval[i].size(); ++j)
            result[i][j] += representation[i][j];
    }
    return result;
}

double flux_error_decomposition(double observation_error, double representation_error,
                                double transport_error) {
    return std::sqrt(observation_error * observation_error +
                     representation_error * representation_error +
                     transport_error * transport_error);
}

double ocean_colour_variance(const std::vector<double>& jacobian, const Matrix& covariance) {
    const std::size_t n = jacobian.size();
    require(covariance.size() == n, "ocean_colour_variance: size mismatch");
    double variance = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        require(covariance[i].size() == n, "ocean_colour_variance: ragged");
        for (std::size_t j = 0; j < n; ++j)
            variance += jacobian[i] * covariance[i][j] * jacobian[j];
    }
    require(variance >= 0.0, "ocean_colour_variance: negative variance");
    return variance;
}

double measurement_uncertainty(double data_uncertainty, double model_uncertainty) {
    return std::sqrt(data_uncertainty * data_uncertainty +
                     model_uncertainty * model_uncertainty);
}

double partially_linear_effect(const std::vector<double>& residual_outcome,
                               const std::vector<double>& residual_treatment) {
    require(!residual_outcome.empty() && residual_outcome.size() == residual_treatment.size(),
            "partially_linear_effect: size mismatch");
    double covariance = 0.0, variance = 0.0;
    for (std::size_t i = 0; i < residual_outcome.size(); ++i) {
        covariance += residual_treatment[i] * residual_outcome[i];
        variance += residual_treatment[i] * residual_treatment[i];
    }
    require(variance > 0.0, "partially_linear_effect: zero treatment variance");
    return covariance / variance;
}

}  // namespace climate
