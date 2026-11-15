#ifndef CLIMATE_ATTRIB2_HPP
#define CLIMATE_ATTRIB2_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

// ---- Optimal fingerprinting ---------------------------------------------

struct FingerprintFit {
    std::vector<double> scaling;
    Matrix covariance;
    std::vector<double> standard_error;
};

FingerprintFit gls_fingerprint(const std::vector<double>& observations,
                               const Matrix& signals, const Matrix& noise_precision);

double detection_z(double scaling, double standard_error);

bool attribution_consistent(double scaling, double standard_error, double alpha);

// Returns the uniquely attributed index, or -1 when the criterion fails.
int attribution_decision(const std::vector<double>& confidence_levels, double critical);

double fingerprint_signal_to_noise(const std::vector<double>& fingerprint,
                                   const Matrix& noise_precision);

double joint_signal_to_noise(const std::vector<double>& first,
                             const std::vector<double>& second, const Matrix& first_precision,
                             const Matrix& second_precision, const Matrix& cross_covariance);

double joint_signal_to_noise_schur(const std::vector<double>& first,
                                   const std::vector<double>& second,
                                   const Matrix& first_covariance,
                                   const Matrix& second_covariance,
                                   const Matrix& cross_covariance);

// ---- Errors-in-variables estimating equations ---------------------------

struct EeFit {
    std::vector<double> scaling;
    Matrix precision;  // A_T
    double condition_number = 0.0;
};

EeFit estimating_equation_fit(const std::vector<double>& observations,
                              const Matrix& signals, const Matrix& noise_precision,
                              double scale_parameter, double spatial_dimension,
                              const std::vector<double>& model_runs);

double ee_bias_correction(double scale_parameter, double spatial_dimension,
                          const std::vector<double>& model_runs, std::size_t coefficient);

double ee_pseudo_residual_scale(const std::vector<double>& scaling,
                                const std::vector<double>& model_runs,
                                double scale_parameter);

double confidence_interval_width(double standard_error, double alpha);

// ---- Box-Cox nonlinear co-trending --------------------------------------

double box_cox(double value, double lambda);

double box_cox_lambda_derivative(double value, double lambda);

double box_cox_marginal_effect(double beta, double value, double lambda);

double box_cox_marginal_effect_derivative(double beta, double value, double lambda);

// ---- Bierens nonparametric co-trending ----------------------------------

struct BierensMatrices {
    Matrix partial_sums;
    Matrix moving_average;
};

BierensMatrices bierens_matrices(const Matrix& series, double alpha);

std::vector<double> bierens_eigenvalues(const Matrix& series, double alpha);

double bierens_statistic(const Matrix& series, double alpha, std::size_t r);

double bierens_eigenvalue_sensitivity(const std::vector<double>& right_vector,
                                      const std::vector<double>& left_vector,
                                      const Matrix& perturbed_first,
                                      const Matrix& perturbed_second, double eigenvalue,
                                      const Matrix& second_matrix);

}  // namespace climate

#endif
