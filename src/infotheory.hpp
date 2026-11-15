#ifndef CLIMATE_INFOTHEORY_HPP
#define CLIMATE_INFOTHEORY_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

double digamma(double value);

// Histogram mutual information with equal-width bins.
double mutual_information_plugin(const std::vector<double>& x, const std::vector<double>& y,
                                 std::size_t bins);

double miller_madow_correction(std::size_t x_bins, std::size_t y_bins, std::size_t samples);

// Conditional MI by histogram stratification on the conditioning columns.
double partial_mutual_information_plugin(const std::vector<double>& x,
                                         const std::vector<double>& y,
                                         const Matrix& conditioning, std::size_t bins);

// Kraskov-Stogbauer-Grassberger k-NN estimator (max norm).
double ksg_mutual_information(const std::vector<double>& x, const std::vector<double>& y,
                              std::size_t neighbours);

// Chi-square calibrated PMI threshold for discrete counts.
double pmi_threshold(double alpha, std::size_t conditioning_states, std::size_t x_bins,
                     std::size_t y_bins, std::size_t samples);

double chi_square_quantile(double probability, double degrees_of_freedom);

double optimal_bin_width(std::size_t samples, std::size_t dimensions, double c1, double c2);

double mutual_information_influence(double log_density_ratio, double mutual_information);

double mutual_information_influence_variance(const std::vector<double>& log_density_ratios);

double gaussian_partial_mutual_information(double delta_r_squared);

// Redundancy = sum of marginal MI minus joint MI.
double information_redundancy(const std::vector<double>& marginal_information,
                              double joint_information);

// Forward selection with re-testing of excluded variables at every step.
std::vector<std::size_t> pmi_forward_selection(const std::vector<double>& response,
                                               const Matrix& predictors, double threshold,
                                               std::size_t bins);

}  // namespace climate

#endif
