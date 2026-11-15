#ifndef CLIMATE_KERNEL_HPP
#define CLIMATE_KERNEL_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

std::vector<double> retrieved_trend_profile(const Matrix& kernel,
                                            const std::vector<double>& true_trend);

Matrix trend_covariance(const Matrix& kernel, const Matrix& trend_covariance,
                        const Matrix& error_covariance);

double toe_relative_error(double toe_with_kernel, double toe_without_kernel,
                          double natural_variability);

double kernel_sensitivity(const Matrix& kernel, std::size_t retrieved_level,
                          std::size_t true_level);

std::vector<double> apriori_trend_sensitivity(const Matrix& kernel,
                                              const std::vector<double>& trend_apriori_derivative);

}  // namespace climate

#endif
