#ifndef CLIMATE_UNCERTAINTY_HPP
#define CLIMATE_UNCERTAINTY_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// Independent propagation: u_c = sqrt(sum (c_i u_i)^2).
double combined_standard_uncertainty(const std::vector<double>& sensitivities,
                                     const std::vector<double>& standard_uncertainties);

// Correlated propagation: u_c = sqrt(sum_ij c_i c_j rho_ij u_i u_j).
double correlated_standard_uncertainty(const std::vector<double>& sensitivities,
                                       const std::vector<double>& standard_uncertainties,
                                       const Matrix& correlation);

// Variance of a spatial mean with correlated per-point errors.
double correlated_mean_variance(const std::vector<double>& standard_uncertainties,
                                const Matrix& correlation);

// Retrieval plus sampling: u_total = sqrt(u_retrieval^2 + u_sampling^2).
double two_component_error(double retrieval_uncertainty, double sampling_uncertainty);

// chi-square / degrees of freedom.
double reduced_chi_square(const std::vector<double>& observed,
                          const std::vector<double>& reference,
                          const std::vector<double>& uncertainty);

double root_mean_square_error(const std::vector<double>& observed,
                              const std::vector<double>& reference);

double mean_bias(const std::vector<double>& observed, const std::vector<double>& reference);

double normalised_mean_bias(const std::vector<double>& observed,
                            const std::vector<double>& reference);

}  // namespace climate

#endif
