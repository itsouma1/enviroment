#ifndef CLIMATE_CHANGE_HPP
#define CLIMATE_CHANGE_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// ---------------------------------------------------------------------------
// Piecewise linear changepoints
// ---------------------------------------------------------------------------
double piecewise_linear_value(double intercept, double rate, double rate_change, double time,
                              double change_point);

double piecewise_rate_before(double rate);

double piecewise_rate_after(double rate, double rate_change);

// ---------------------------------------------------------------------------
// Non-stationary extremes
// ---------------------------------------------------------------------------
// G^-1(1-p) = mu + sigma/xi * {[-ln(1-p)]^-xi - 1}.
double effective_return_level(double probability, double location, double scale, double shape);

// dT/dt = Delta tau / (1 - F)^2 * dF/dt.
double return_period_derivative(double observation_interval, double cdf,
                                double cdf_derivative);

// ---------------------------------------------------------------------------
// AMOC stability and heat transport
// ---------------------------------------------------------------------------
double ar2_step(double first_lag, double second_lag, double phi1, double phi2, double shock);

double ar2_stability(double phi1, double phi2);

double oht_amplitude_term(double velocity_prime, double temperature_mean);

double oht_circulation_term(double velocity_mean, double temperature_prime);

double oht_covariance_term(double velocity_prime, double temperature_prime);

// ---------------------------------------------------------------------------
// Information geometry
// ---------------------------------------------------------------------------
double fisher_information_gaussian_mean(double variance);

double information_speed(double fisher_information, double parameter_rate);

// ---------------------------------------------------------------------------
// Sparse change detection
// ---------------------------------------------------------------------------
double clipped_standardised_value(double value, double median, double median_absolute_deviation,
                                  double threshold);

double local_contrast(double left_mean, double right_mean, double left_variance,
                      double right_variance, double regularisation, double half_width);

double max_cusum_statistic(const std::vector<double>& contrasts);

double higher_criticism(const std::vector<double>& p_values);

// ---------------------------------------------------------------------------
// Distributional and covariance change
// ---------------------------------------------------------------------------
double gaussian_covariance_log_likelihood_ratio(int first_count, int second_count,
                                               const Matrix& pooled_covariance,
                                               const Matrix& first_covariance,
                                               const Matrix& second_covariance);

double optimal_transport_cost_1d(const std::vector<double>& from, const std::vector<double>& to);

// ---------------------------------------------------------------------------
// Functional slope
// ---------------------------------------------------------------------------
std::vector<double> functional_slope(const Matrix& curves, const std::vector<double>& times);

double functional_slope_norm(const std::vector<double>& slope, double spacing);

}  // namespace climate

#endif
