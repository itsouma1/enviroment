#include "uncertainty.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double combined_standard_uncertainty(const std::vector<double>& sensitivities,
                                     const std::vector<double>& standard_uncertainties) {
    require(sensitivities.size() == standard_uncertainties.size(),
            "combined_standard_uncertainty: size mismatch");
    double variance = 0.0;
    for (std::size_t i = 0; i < sensitivities.size(); ++i) {
        const double contribution = sensitivities[i] * standard_uncertainties[i];
        variance += contribution * contribution;
    }
    return std::sqrt(variance);
}

double correlated_standard_uncertainty(const std::vector<double>& sensitivities,
                                       const std::vector<double>& standard_uncertainties,
                                       const Matrix& correlation) {
    const std::size_t n = sensitivities.size();
    require(standard_uncertainties.size() == n && correlation.size() == n,
            "correlated_standard_uncertainty: size mismatch");
    double variance = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        require(correlation[i].size() == n, "correlated_standard_uncertainty: ragged");
        for (std::size_t j = 0; j < n; ++j)
            variance += sensitivities[i] * sensitivities[j] * correlation[i][j] *
                        standard_uncertainties[i] * standard_uncertainties[j];
    }
    return std::sqrt(std::max(variance, 0.0));
}

double correlated_mean_variance(const std::vector<double>& standard_uncertainties,
                                const Matrix& correlation) {
    const std::size_t n = standard_uncertainties.size();
    require(n > 0 && correlation.size() == n, "correlated_mean_variance: size mismatch");
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        require(correlation[i].size() == n, "correlated_mean_variance: ragged");
        for (std::size_t j = 0; j < n; ++j)
            total += correlation[i][j] * standard_uncertainties[i] * standard_uncertainties[j];
    }
    return total / (static_cast<double>(n) * static_cast<double>(n));
}

double two_component_error(double retrieval_uncertainty, double sampling_uncertainty) {
    return std::sqrt(retrieval_uncertainty * retrieval_uncertainty +
                     sampling_uncertainty * sampling_uncertainty);
}

double reduced_chi_square(const std::vector<double>& observed,
                          const std::vector<double>& reference,
                          const std::vector<double>& uncertainty) {
    require(!observed.empty() && observed.size() == reference.size() &&
                observed.size() == uncertainty.size(),
            "reduced_chi_square: size mismatch");
    require(observed.size() > 1, "reduced_chi_square: need at least two values");
    double sum = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        require(uncertainty[i] > 0.0, "reduced_chi_square: non-positive uncertainty");
        const double standardised = (observed[i] - reference[i]) / uncertainty[i];
        sum += standardised * standardised;
    }
    return sum / static_cast<double>(observed.size() - 1);
}

double root_mean_square_error(const std::vector<double>& observed,
                              const std::vector<double>& reference) {
    require(!observed.empty() && observed.size() == reference.size(),
            "root_mean_square_error: size mismatch");
    double sum = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        const double difference = observed[i] - reference[i];
        sum += difference * difference;
    }
    return std::sqrt(sum / static_cast<double>(observed.size()));
}

double mean_bias(const std::vector<double>& observed, const std::vector<double>& reference) {
    require(!observed.empty() && observed.size() == reference.size(),
            "mean_bias: size mismatch");
    double sum = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) sum += observed[i] - reference[i];
    return sum / static_cast<double>(observed.size());
}

double normalised_mean_bias(const std::vector<double>& observed,
                            const std::vector<double>& reference) {
    require(!observed.empty() && observed.size() == reference.size(),
            "normalised_mean_bias: size mismatch");
    double reference_sum = 0.0;
    for (double value : reference) reference_sum += value;
    require(std::fabs(reference_sum) > 0.0, "normalised_mean_bias: zero reference sum");
    return (mean_bias(observed, reference) * static_cast<double>(observed.size())) /
           reference_sum;
}

}  // namespace climate
