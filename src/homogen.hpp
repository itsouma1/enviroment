#ifndef CLIMATE_HOMOGEN_HPP
#define CLIMATE_HOMOGEN_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// ---------------------------------------------------------------------------
// Benchmark validation
// ---------------------------------------------------------------------------
double absolute_mean_bias(const std::vector<double>& homogenised,
                          const std::vector<double>& truth);

double factor_of_exceedance(const std::vector<double>& homogenised,
                            const std::vector<double>& corrupted);

double relative_rmse(const std::vector<double>& homogenised,
                     const std::vector<double>& raw, const std::vector<double>& truth);

// Bivariate ACMANT objective with weighting c0^2.
double acmant_objective(const std::vector<double>& first, const std::vector<double>& second,
                        const std::vector<int>& break_indices, double weight);

// Ensemble mean and variance of homogenisation members.
void homogenisation_ensemble(const Matrix& members, std::vector<double>& mean,
                             std::vector<double>& variance);

// ---------------------------------------------------------------------------
// BYS-CUSUM
// ---------------------------------------------------------------------------
double residual_cusum(const std::vector<double>& residuals, std::size_t candidate);

double student_t_log_likelihood(const std::vector<double>& residuals, double degrees_of_freedom,
                                double scale);

// ---------------------------------------------------------------------------
// Heavy-tailed and regime models
// ---------------------------------------------------------------------------
double singh_maddala_density(double value, double a, double b, double q);

double sparse_jump_objective(const Matrix& observations, const std::vector<int>& states,
                             const Matrix& state_parameters, double jump_penalty);

// ---------------------------------------------------------------------------
// Selection and spatial diagnostics
// ---------------------------------------------------------------------------
double cautious_correction(double reported, double estimated_bias);

double surface_urban_heat_island(double baseline, double nonurban_minimum);

double excursion_probability(double trend, double domain_trend, double standard_error);

double variational_lower_bound(double expected_log_likelihood, double kl_divergence);

}  // namespace climate

#endif
