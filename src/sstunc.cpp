#include "sstunc.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "regional.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Matrix cholesky(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    Matrix lower(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k) sum -= lower[i][k] * lower[j][k];
            if (i == j) {
                require(sum >= -1e-12, "cholesky: not positive semidefinite");
                lower[i][j] = std::sqrt(std::max(0.0, sum));
            } else {
                lower[i][j] = sum / lower[j][j];
            }
        }
    return lower;
}

}  // namespace

double uncorrelated_variance(const std::vector<double>& weights,
                             const std::vector<double>& sigmas) {
    require(weights.size() == sigmas.size() && !weights.empty(),
            "uncorrelated_variance: size");
    double variance = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i)
        variance += weights[i] * weights[i] * sigmas[i] * sigmas[i];
    return variance;
}

double annual_average_variance(const std::vector<double>& monthly_sigmas) {
    require(monthly_sigmas.size() == 12, "annual_average_variance: 12 months");
    double variance = 0.0;
    for (double sigma : monthly_sigmas)
        variance += (1.0 / 12.0) * (1.0 / 12.0) * sigma * sigma;
    return variance;
}

double covariance_variance(const std::vector<double>& weights, const Matrix& covariance) {
    require(!weights.empty() && covariance.size() == weights.size() &&
                covariance[0].size() == weights.size(),
            "covariance_variance: size");
    double variance = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i)
        for (std::size_t j = 0; j < weights.size(); ++j)
            variance += weights[i] * weights[j] * covariance[i][j];
    return variance;
}

double hemispheric_difference_variance(const std::vector<double>& north_weights,
                                       const std::vector<double>& south_weights,
                                       const Matrix& covariance) {
    require(north_weights.size() == south_weights.size() &&
                covariance.size() == north_weights.size(),
            "hemispheric_difference_variance: size");
    std::vector<double> difference(north_weights.size(), 0.0);
    for (std::size_t i = 0; i < difference.size(); ++i)
        difference[i] = north_weights[i] - south_weights[i];
    return covariance_variance(difference, covariance);
}

double constant_covariance_variance(double sigma, double correlation,
                                    const std::vector<double>& weights) {
    require(sigma >= 0.0 && correlation >= -1.0 && correlation <= 1.0,
            "constant_covariance_variance: parameters");
    double sum = 0.0, sum_squares = 0.0;
    for (double weight : weights) {
        sum += weight;
        sum_squares += weight * weight;
    }
    return sigma * sigma *
           (correlation * sum * sum + (1.0 - correlation) * sum_squares);
}

double effective_sample_size(const std::vector<double>& weights, double correlation) {
    require(!weights.empty() && correlation >= 0.0 && correlation < 1.0,
            "effective_sample_size: inputs");
    double sum = 0.0, sum_squares = 0.0;
    for (double weight : weights) {
        sum += weight;
        sum_squares += weight * weight;
    }
    const double correlated = correlation * sum * sum +
                              (1.0 - correlation) * sum_squares;
    require(correlated > 0.0, "effective_sample_size: zero variance");
    return static_cast<double>(weights.size()) * sum_squares / correlated;
}

Matrix sample_correlated_errors(const Matrix& covariance, std::size_t samples,
                                unsigned seed) {
    require(!covariance.empty() && covariance.size() == covariance[0].size() &&
                samples > 0,
            "sample_correlated_errors: inputs");
    const Matrix lower = cholesky(covariance);
    const std::size_t n = covariance.size();
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    Matrix result(samples, std::vector<double>(n, 0.0));
    for (std::size_t s = 0; s < samples; ++s) {
        std::vector<double> draw(n, 0.0);
        for (double& value : draw) value = normal(generator);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k <= i; ++k)
                result[s][i] += lower[i][k] * draw[k];
    }
    return result;
}

double ensemble_standard_deviation(const std::vector<double>& values) {
    require(values.size() >= 2, "ensemble_standard_deviation: size");
    double mean = 0.0;
    for (double value : values) mean += value;
    mean /= static_cast<double>(values.size());
    double variance = 0.0;
    for (double value : values) variance += (value - mean) * (value - mean);
    return std::sqrt(variance / static_cast<double>(values.size() - 1));
}

double combine_uncertainty(double uncorrelated, double simply_correlated, double ensemble) {
    require(uncorrelated >= 0.0 && simply_correlated >= 0.0 && ensemble >= 0.0,
            "combine_uncertainty: components");
    return std::sqrt(uncorrelated * uncorrelated + simply_correlated * simply_correlated +
                     ensemble * ensemble);
}

double uncorrelated_weight_sensitivity(double weight, double sigma) {
    return 2.0 * weight * sigma * sigma;
}

double covariance_entry_sensitivity(double first_weight, double second_weight,
                                    bool diagonal) {
    return diagonal ? first_weight * first_weight : first_weight * second_weight;
}

double correlation_sensitivity(double sigma, const std::vector<double>& weights) {
    double sum_squares = 0.0, sum = 0.0;
    for (double weight : weights) {
        sum_squares += weight * weight;
        sum += weight;
    }
    return sigma * sigma * (sum * sum - sum_squares);
}

double total_component_sensitivity(double component, double total) {
    require(total > 0.0, "total_component_sensitivity: total");
    return component / total;
}

double covariance_condition_number(const Matrix& covariance) {
    require(!covariance.empty() && covariance.size() == covariance[0].size(),
            "covariance_condition_number: square");
    const EigenSystem eigen = symmetric_eigen(covariance);
    const double largest = *std::max_element(eigen.values.begin(), eigen.values.end());
    const double smallest = *std::min_element(eigen.values.begin(), eigen.values.end());
    require(smallest > 0.0, "covariance_condition_number: singular");
    return largest / smallest;
}

}  // namespace climate
