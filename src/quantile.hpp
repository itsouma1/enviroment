#ifndef CLIMATE_QUANTILE_HPP
#define CLIMATE_QUANTILE_HPP

#include <cstddef>
#include <random>
#include <vector>

namespace climate {

double check_loss(double residual, double probability);

double sample_quantile(std::vector<double> values, double probability);

struct QuantileFit {
    double slope = 0.0;
    double intercept = 0.0;
};

// Iteratively reweighted least squares for the check-loss regression.
QuantileFit quantile_trend(const std::vector<double>& times, const std::vector<double>& values,
                           double probability);

double ald_density(double residual, double probability, double scale);

double ald_log_likelihood(const std::vector<double>& residuals, double probability,
                          double scale);

double ald_subgradient(double residual, double probability);

// Potts conditional probabilities for one site given neighbour counts.
std::vector<double> potts_conditional_probabilities(const std::vector<double>& neighbour_counts,
                                                    const std::vector<double>& alpha,
                                                    const std::vector<double>& beta);

double potts_log_odds_neighbour_derivative(double beta);

// One sequential Gibbs sweep over a rows x cols label grid.
void potts_gibbs_sweep(std::vector<int>& labels, std::size_t rows, std::size_t cols,
                       const std::vector<double>& alpha, const std::vector<double>& beta,
                       std::mt19937& generator);

// One Swendsen-Wang global update with edge activation 1 - exp(-beta[label]).
void potts_swendsen_wang_step(std::vector<int>& labels, std::size_t rows, std::size_t cols,
                              const std::vector<double>& alpha, const std::vector<double>& beta,
                              std::mt19937& generator);

}  // namespace climate

#endif
