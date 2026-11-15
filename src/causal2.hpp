#ifndef CLIMATE_CAUSAL2_HPP
#define CLIMATE_CAUSAL2_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

struct CausalLink {
    std::size_t driver = 0;
    std::size_t target = 0;
    std::size_t lag = 1;
    double partial_correlation = 0.0;
};

double partial_correlation(const std::vector<double>& x, const std::vector<double>& y,
                           const Matrix& conditioning);

double fisher_z_p_value(double correlation, std::size_t samples);

// Lagged PC skeleton with partial-correlation CI tests.
std::vector<CausalLink> pcmci_lagged_links(const Matrix& series, std::size_t max_lag,
                                           double alpha);

double mci_partial_correlation(double coupling, double driver_sd, double target_sd);

// Network-level MK sum for a T x N field.
double network_mann_kendall_sum(const Matrix& field);

double rcamk_p_value(const Matrix& field, std::size_t resamples, double ar1,
                     unsigned seed);

double gaussian_vine_full_loglik(const Matrix& correlation, std::size_t samples);

double gaussian_vine_truncated_loglik(const Matrix& correlation, std::size_t samples);

double higher_tree_diagnostic(const Matrix& correlation, std::size_t samples);

}  // namespace climate

#endif
