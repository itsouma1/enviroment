#ifndef CLIMATE_DETECT_HPP
#define CLIMATE_DETECT_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// ---------------------------------------------------------------------------
// Extreme value theory
// ---------------------------------------------------------------------------
// GEV CDF with location mu, scale sigma and shape xi.
double gev_cdf(double value, double location, double scale, double shape);

// Return level for return period 1/p.
double gev_return_level(double probability, double location, double scale, double shape);

// GEV log-likelihood for one observation.
double gev_log_likelihood(double value, double location, double scale, double shape);

// Generalized Pareto CDF for exceedance y above the threshold.
double gp_cdf(double exceedance, double scale, double shape);

// ---------------------------------------------------------------------------
// Early warning (Ornstein-Uhlenbeck)
// ---------------------------------------------------------------------------
double ou_stationary_variance(double noise_amplitude, double damping_rate);

double ou_lag1_autocorrelation(double damping_rate, double time_step);

double ar1_estimate(const std::vector<double>& values);

// Rising variance and AR1 combined into one operational score (sum of the two
// standardised trends).
double rate_mismatch_criterion(double variance_slope, double autocorrelation_slope,
                               double variance_scale, double autocorrelation_scale);

// ---------------------------------------------------------------------------
// Regime shifts
// ---------------------------------------------------------------------------
// RSI = (1/(m sigma)) sum_{i=0}^{m} (x_{t+i} - mean_new).
double regime_shift_index(const std::vector<double>& values, double new_mean,
                          double regime_sigma);

// Maximum of the normalised cumulative sum.
double cusum_statistic(const std::vector<double>& values);

// ---------------------------------------------------------------------------
// Optimal fingerprinting
// ---------------------------------------------------------------------------
// beta = (X^T C^-1 X)^-1 X^T C^-1 y.
std::vector<double> optimal_fingerprint(const std::vector<std::vector<double>>& responses,
                                        const std::vector<double>& observation,
                                        const std::vector<std::vector<double>>& noise_covariance);

// ---------------------------------------------------------------------------
// Random matrix theory
// ---------------------------------------------------------------------------
double marchenko_pastur_bounds_lower(double variance, double aspect_ratio);
double marchenko_pastur_bounds_upper(double variance, double aspect_ratio);
double marchenko_pastur_density(double eigenvalue, double variance, double aspect_ratio);

// Number of eigenvalues above the Marchenko-Pastur upper bound.
int significant_eigenvalue_count(const std::vector<double>& eigenvalues, double variance,
                                 double aspect_ratio);

// C = sum_gamma lambda_gamma u_gamma u_gamma^T; eigenvectors[gamma] is u_gamma.
Matrix rmt_reconstruction(const std::vector<double>& eigenvalues,
                          const std::vector<std::vector<double>>& eigenvectors);

// Total least squares: smallest right singular vector of [X | y].
std::vector<double> total_least_squares(const Matrix& predictors,
                                        const std::vector<double>& response);

}  // namespace climate

#endif
