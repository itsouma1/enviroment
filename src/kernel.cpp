#include "kernel.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

std::vector<double> retrieved_trend_profile(const Matrix& kernel,
                                            const std::vector<double>& true_trend) {
    require(!kernel.empty() && kernel.size() == true_trend.size(),
            "retrieved_trend_profile: size");
    std::vector<double> retrieved(kernel.size(), 0.0);
    for (std::size_t i = 0; i < kernel.size(); ++i) {
        require(kernel[i].size() == true_trend.size(),
                "retrieved_trend_profile: ragged kernel");
        for (std::size_t j = 0; j < true_trend.size(); ++j)
            retrieved[i] += kernel[i][j] * true_trend[j];
    }
    return retrieved;
}

Matrix trend_covariance(const Matrix& kernel, const Matrix& trend_covariance,
                        const Matrix& error_covariance) {
    const std::size_t n = kernel.size();
    require(n > 0 && trend_covariance.size() == n && error_covariance.size() == n,
            "trend_covariance: size");
    Matrix product(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t k = 0; k < n; ++k)
                for (std::size_t l = 0; l < n; ++l)
                    product[i][j] += kernel[i][k] * trend_covariance[k][l] *
                                     kernel[j][l];
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) product[i][j] += error_covariance[i][j];
    return product;
}

double toe_relative_error(double toe_with_kernel, double toe_without_kernel,
                          double natural_variability) {
    require(std::fabs(natural_variability) > 0.0, "toe_relative_error: variability");
    return (toe_with_kernel - toe_without_kernel) / natural_variability;
}

double kernel_sensitivity(const Matrix& kernel, std::size_t retrieved_level,
                          std::size_t true_level) {
    require(retrieved_level < kernel.size() && true_level < kernel[retrieved_level].size(),
            "kernel_sensitivity: index");
    return kernel[retrieved_level][true_level];
}

std::vector<double> apriori_trend_sensitivity(const Matrix& kernel,
                                              const std::vector<double>& trend_apriori_derivative) {
    const std::size_t n = kernel.size();
    require(n > 0 && trend_apriori_derivative.size() == n,
            "apriori_trend_sensitivity: size");
    std::vector<double> result(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const double identity = i == j ? 1.0 : 0.0;
            result[i] += (identity - kernel[i][j]) * trend_apriori_derivative[j];
        }
    return result;
}

}  // namespace climate
