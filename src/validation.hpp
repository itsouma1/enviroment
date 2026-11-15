#ifndef CLIMATE_VALIDATION_HPP
#define CLIMATE_VALIDATION_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// Mean squared standardised difference against an independent reference.
double standardised_difference_variance(const std::vector<double>& observed,
                                        const std::vector<double>& reference,
                                        double observation_uncertainty,
                                        double reference_uncertainty,
                                        double representation_error);

// Quadruple-collocation closure: e_ij = C_ij / (a_i a_j) - T.
double quadruple_collocation_covariance(double covariance, double scale_i, double scale_j,
                                        double truth_variance);

double linear_ice_fraction(double ice_length, double total_length);

double sea_ice_quadratic_bias(double a0, double a1, double a2, double sea_ice_concentration);

double sea_ice_bias_derivative(double a1, double a2, double sea_ice_concentration);

// Inverse-variance merge.
double inverse_variance_merge(const std::vector<double>& values,
                              const std::vector<double>& variances);

double precipitation_bias_sum(double hit_bias, double missed_bias, double false_bias);

// R = R_retrieval + R_representation.
Matrix observation_error_covariance(const Matrix& retrieval, const Matrix& representation);

// RMSE^2 = ERR_OBS^2 + ERR_REP^2 + ERR_MIP^2.
double flux_error_decomposition(double observation_error, double representation_error,
                                double transport_error);

// u^2(y) = J^T C J.
double ocean_colour_variance(const std::vector<double>& jacobian, const Matrix& covariance);

double measurement_uncertainty(double data_uncertainty, double model_uncertainty);

// Partially linear CATE residual estimator.
double partially_linear_effect(const std::vector<double>& residual_outcome,
                               const std::vector<double>& residual_treatment);

}  // namespace climate

#endif
