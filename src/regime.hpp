#ifndef CLIMATE_REGIME_HPP
#define CLIMATE_REGIME_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// ---------------------------------------------------------------------------
// Continuous-time processes
// ---------------------------------------------------------------------------
// Symmetric alpha-stable increment given one uniform and one exponential draw.
double symmetric_stable_sample(double alpha, double uniform, double exponential);

// ---------------------------------------------------------------------------
// Climate feedback manifold
// ---------------------------------------------------------------------------
double climate_feedback_manifold(double lambda, double kappa, double temperature);

double temperature_tendency(double forcing, double lambda, double kappa, double temperature);

// ---------------------------------------------------------------------------
// Smooth transition regression
// ---------------------------------------------------------------------------
double logistic_transition(double covariate, double gamma, double threshold);

// Exponential (LSTR2) transition, bounded in [0, 1).
double exponential_transition(double covariate, double gamma, double threshold);

// y = base + regime1 G1 + regime2 G2.
double multiple_lstr(double base, double regime1, double regime2, double transition1,
                     double transition2);

// ---------------------------------------------------------------------------
// Stochastic volatility
// ---------------------------------------------------------------------------
double stochastic_volatility_observation(double log_volatility, double shock);

double stochastic_volatility_step(double log_volatility, double long_run_mean, double persistence,
                                  double volatility_of_volatility, double shock);

// Discretised-variance estimates: K = -ln(phi)/dt, gamma0 = theta/(1 - phi).
double stochastic_volatility_mean_reversion(double autoregressive_coefficient, double time_step);

double stochastic_volatility_intercept(double intercept, double autoregressive_coefficient);

// ---------------------------------------------------------------------------
// GARCH family
// ---------------------------------------------------------------------------
double gjr_garch_variance(double constant, double shock_coefficient,
                          double asymmetry_coefficient, double persistence,
                          double previous_variance, double previous_shock);

double beta_weight(int lag, int maximum_lag, double omega1, double omega2);

double long_run_log_variance(double intercept, double slope,
                             const std::vector<double>& weights,
                             const std::vector<double>& covariate_lags);

double garch_midas_variance(double previous_component, double previous_return,
                            double previous_long_run, double alpha, double beta);

// ---------------------------------------------------------------------------
// Quantile regimes and seasonal trends
// ---------------------------------------------------------------------------
double asymmetric_laplace_density(double value, double location, double scale,
                                  double quantile);

double average_seasonal_trends(const std::vector<double>& seasonal_trends);

double seasonal_trend_standard_error(const Matrix& covariance);

// ---------------------------------------------------------------------------
// BEAST-style decomposition
// ---------------------------------------------------------------------------
double beast_seasonal(double time, double period, const std::vector<double>& cosines,
                      const std::vector<double>& sines);

double beast_steps(double time, const std::vector<double>& change_points,
                   const std::vector<double>& magnitudes);

}  // namespace climate

#endif
