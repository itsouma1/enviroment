#ifndef CLIMATE_MULTITREND_HPP
#define CLIMATE_MULTITREND_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

// ---- Mann-Kendall building blocks ---------------------------------------

double mann_kendall_score(const std::vector<double>& values);

double mann_kendall_variance(std::size_t n);

double mann_kendall_rank(const std::vector<double>& values, std::size_t index);

double mann_kendall_concordance(const std::vector<double>& first,
                                const std::vector<double>& second);

// Cross-covariance Gamma_xy of the two MK scores (Libiseller-Grimvall).
double mann_kendall_cross_covariance(const std::vector<double>& first,
                                     const std::vector<double>& second);

double multivariate_mann_kendall(const std::vector<double>& scores,
                                 const Matrix& covariance);

// S^T Gamma^-1 S, asymptotically chi-square with d degrees of freedom.
double multivariate_mann_kendall_quadratic(const std::vector<double>& scores,
                                           const Matrix& covariance);

// Partial MK: trend in `target` after conditioning on `covariate`.
double partial_mann_kendall(const std::vector<double>& target,
                            const std::vector<double>& covariate);

// Seasonal block scores: each row is one season, columns are years.
std::vector<double> seasonal_kendall_scores(const Matrix& seasonal_data);

double seasonal_kendall_variance(const Matrix& seasonal_data);

// Hirsch-Slack cross-covariance from two seasonal score matrices.
double seasonal_kendall_cross_covariance(const Matrix& first, const Matrix& second,
                                         bool missing_values);

// ---- Moving-window dependence trends ------------------------------------

std::vector<double> moving_window_kendall(const std::vector<double>& first,
                                          const std::vector<double>& second,
                                          std::size_t window);

double dependence_trend(const std::vector<double>& first, const std::vector<double>& second,
                        std::size_t window);

double overall_multivariate_trend(const std::vector<double>& first,
                                  const std::vector<double>& second, std::size_t window);

double dependence_trend_p_value(const std::vector<double>& first,
                                const std::vector<double>& second, std::size_t window,
                                std::size_t resamples, unsigned seed);

}  // namespace climate

#endif
