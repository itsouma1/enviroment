#pragma once

#include <cstddef>
#include <vector>

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
