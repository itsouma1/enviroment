#pragma once

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

using Weights = std::vector<std::vector<double>>;

// ---------------------------------------------------------------------------
// Local spatial statistics
// ---------------------------------------------------------------------------
// Getis-Ord Gi*:
//   Gi* = (sum_j w_ij x_j - xbar sum_j w_ij) /
//         (S sqrt[(n sum_j w_ij^2 - (sum_j w_ij)^2) / (n - 1)])
// with S^2 = sum_j x_j^2 / n - xbar^2. The weights normally include the
// location itself (i = j).
double getis_ord_gi(const std::vector<double>& values, const Weights& weights, std::size_t index);

// Local Moran's I:
//   I_i = (x_i - xbar) / m2 * sum_j w_ij (x_j - xbar),  m2 = sum_j (x_j - xbar)^2 / n
double local_moran_i(const std::vector<double>& values, const Weights& weights,
                     std::size_t index);

// Geo-Detector q-statistic:
//   q = 1 - sum_h n_h sigma_h^2 / (N sigma^2)
// Strata labels are arbitrary integers.
double geodetector_q(const std::vector<double>& values, const std::vector<int>& strata);

// LISA regime from the sign of the local statistic and its neighbourhood.
enum class LisaClass { HighHigh, LowLow, HighLow, LowHigh };
LisaClass lisa_classification(const std::vector<double>& values, const Weights& weights,
                              std::size_t index);

// Gaussian kernel weight for a distance and bandwidth.
double gaussian_kernel(double distance, double bandwidth);

// Geographically weighted regression coefficients at one location:
// beta = (X^T W X)^-1 X^T W y.
std::vector<double> gwr_coefficients(const Matrix& predictors, const std::vector<double>& response,
                                     const std::vector<double>& distances, double bandwidth);

// Matern 3/2 covariance.
double matern32_covariance(double distance, double variance, double length_scale);

// ---------------------------------------------------------------------------
// Spatially varying changepoints
// ---------------------------------------------------------------------------
// P(lower < tau <= upper) for a probit threshold model with latent mean mu.
double probit_category_probability(double latent_mean, double lower_threshold,
                                   double upper_threshold, double sigma);

// Hotelling T^2 for a functional scan window.
double hotelling_functional_scan(int inside_count, int outside_count,
                                 const std::vector<double>& inside_mean,
                                 const std::vector<double>& outside_mean,
                                 const Matrix& pooled_covariance);

// CRISP CUSUM contrast over a candidate region (published form).
double crisp_cusum_statistic(double inside_sum, double outside_sum, std::size_t inside_count,
                             std::size_t outside_count);

std::size_t crisp_argmax(const std::vector<double>& statistics);

// CRISP objective: within-region and outside sums of squares plus a penalty.
double crisp_criterion(const std::vector<double>& values, const std::vector<int>& region_mask,
                       double penalty, double lambda);

// Log-Gaussian Cox process log-intensity log lambda = x^T beta + omega.
double log_gaussian_cox_intensity(double intercept, const std::vector<double>& covariates,
                                  const std::vector<double>& coefficients,
                                  double spatial_effect);

// DIC = p_D + Dbar with p_D = Dbar - D(theta_bar).
double deviance_information_criterion(double mean_deviance, double deviance_at_mean);

// ---------------------------------------------------------------------------
// Field significance
// ---------------------------------------------------------------------------
// Walker threshold: p = 1 - (1 - alpha_global)^(1/K).
double walker_threshold(double alpha_global, int hypothesis_count);

// Benjamini-Hochberg adjusted p-values (q-values); monotone from the largest p.
std::vector<double> benjamini_hochberg(const std::vector<double>& p_values);

// Indices of hypotheses rejected by the BH procedure at the given FDR.
std::vector<std::size_t> benjamini_hochberg_rejections(const std::vector<double>& p_values,
                                                       double false_discovery_rate);

// ---------------------------------------------------------------------------
// Kulldorff spatial scan statistic (Poisson model)
// ---------------------------------------------------------------------------
// lambda_Z = (n_Z / mu_Z)^n_Z * ((N - n_Z) / (N - mu_Z))^(N - n_Z), returned
// as a log-likelihood ratio.
double poisson_scan_log_likelihood(int inside_count, double inside_expected, int total_count,
                                   double total_expected);

}  // namespace climate
