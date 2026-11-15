#include "quantile.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double weighted_slope(const std::vector<double>& t, const std::vector<double>& y,
                      const std::vector<double>& w) {
    double sum_w = 0.0, mean_t = 0.0, mean_y = 0.0;
    for (std::size_t i = 0; i < t.size(); ++i) {
        sum_w += w[i];
        mean_t += w[i] * t[i];
        mean_y += w[i] * y[i];
    }
    mean_t /= sum_w;
    mean_y /= sum_w;
    double covariance = 0.0, variance = 0.0;
    for (std::size_t i = 0; i < t.size(); ++i) {
        covariance += w[i] * (t[i] - mean_t) * (y[i] - mean_y);
        variance += w[i] * (t[i] - mean_t) * (t[i] - mean_t);
    }
    require(variance > 0.0, "quantile_trend: zero time variance");
    return covariance / variance;
}

std::vector<double> neighbour_counts(const std::vector<int>& labels, std::size_t rows,
                                     std::size_t cols, std::size_t index, int clusters) {
    std::vector<double> counts(static_cast<std::size_t>(clusters), 0.0);
    const std::size_t row = index / cols, column = index % cols;
    const int dr[4] = {-1, 1, 0, 0};
    const int dc[4] = {0, 0, -1, 1};
    for (int k = 0; k < 4; ++k) {
        const int r = static_cast<int>(row) + dr[k];
        const int c = static_cast<int>(column) + dc[k];
        if (r < 0 || c < 0 || r >= static_cast<int>(rows) || c >= static_cast<int>(cols))
            continue;
        const int label = labels[static_cast<std::size_t>(r) * cols +
                                 static_cast<std::size_t>(c)];
        if (label >= 0 && label < clusters) counts[static_cast<std::size_t>(label)] += 1.0;
    }
    return counts;
}

std::size_t find_root(std::vector<std::size_t>& parent, std::size_t index) {
    while (parent[index] != index) {
        parent[index] = parent[parent[index]];
        index = parent[index];
    }
    return index;
}

void unite(std::vector<std::size_t>& parent, std::size_t a, std::size_t b) {
    a = find_root(parent, a);
    b = find_root(parent, b);
    if (a != b) parent[b] = a;
}

}  // namespace

double check_loss(double residual, double probability) {
    require(probability > 0.0 && probability < 1.0, "check_loss: probability");
    return residual * (probability - (residual < 0.0 ? 1.0 : 0.0));
}

double sample_quantile(std::vector<double> values, double probability) {
    require(!values.empty() && probability > 0.0 && probability < 1.0,
            "sample_quantile: inputs");
    std::sort(values.begin(), values.end());
    const double position = probability * static_cast<double>(values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = std::min(values.size() - 1, lower + 1);
    const double weight = position - static_cast<double>(lower);
    return values[lower] * (1.0 - weight) + values[upper] * weight;
}

QuantileFit quantile_trend(const std::vector<double>& times, const std::vector<double>& values,
                           double probability) {
    require(times.size() == values.size() && times.size() > 1,
            "quantile_trend: size mismatch");
    require(probability > 0.0 && probability < 1.0, "quantile_trend: probability");
    const std::size_t n = times.size();
    std::vector<double> weights(n, 1.0);
    double slope = 0.0, intercept = 0.0;
    for (int iteration = 0; iteration < 100; ++iteration) {
        slope = weighted_slope(times, values, weights);
        std::vector<double> residuals(n);
        for (std::size_t i = 0; i < n; ++i) residuals[i] = values[i] - slope * times[i];
        intercept = sample_quantile(residuals, probability);
        double change = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double residual = values[i] - intercept - slope * times[i];
            const double weight = residual >= 0.0 ? probability : (1.0 - probability);
            const double magnitude = std::max(std::fabs(residual), 1e-8);
            const double updated = weight / magnitude;
            change = std::max(change, std::fabs(updated - weights[i]));
            weights[i] = updated;
        }
        if (change < 1e-10) break;
    }
    QuantileFit fit;
    fit.slope = slope;
    fit.intercept = intercept;
    return fit;
}

double ald_density(double residual, double probability, double scale) {
    require(scale > 0.0 && probability > 0.0 && probability < 1.0, "ald_density: inputs");
    return probability * (1.0 - probability) / scale *
           std::exp(-check_loss(residual / scale, probability));
}

double ald_log_likelihood(const std::vector<double>& residuals, double probability,
                          double scale) {
    require(scale > 0.0, "ald_log_likelihood: scale");
    double total = 0.0;
    for (double residual : residuals)
        total += std::log(ald_density(residual, probability, scale));
    return total;
}

double ald_subgradient(double residual, double probability) {
    require(probability > 0.0 && probability < 1.0, "ald_subgradient: probability");
    if (residual < 0.0) return probability - 1.0;
    if (residual > 0.0) return probability;
    return 0.0;
}

std::vector<double> potts_conditional_probabilities(const std::vector<double>& neighbour_counts,
                                                    const std::vector<double>& alpha,
                                                    const std::vector<double>& beta) {
    require(!alpha.empty() && alpha.size() == beta.size() &&
                neighbour_counts.size() == alpha.size(),
            "potts_conditional_probabilities: size");
    std::vector<double> logits(alpha.size());
    double maximum = -1e300;
    for (std::size_t k = 0; k < alpha.size(); ++k) {
        logits[k] = alpha[k] + beta[k] * neighbour_counts[k];
        maximum = std::max(maximum, logits[k]);
    }
    double total = 0.0;
    for (double& value : logits) {
        value = std::exp(value - maximum);
        total += value;
    }
    for (double& value : logits) value /= total;
    return logits;
}

double potts_log_odds_neighbour_derivative(double beta) { return beta; }

void potts_gibbs_sweep(std::vector<int>& labels, std::size_t rows, std::size_t cols,
                       const std::vector<double>& alpha, const std::vector<double>& beta,
                       std::mt19937& generator) {
    require(labels.size() == rows * cols && alpha.size() == beta.size(),
            "potts_gibbs_sweep: size");
    const int clusters = static_cast<int>(alpha.size());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const std::vector<double> counts =
            neighbour_counts(labels, rows, cols, index, clusters);
        const std::vector<double> probabilities =
            potts_conditional_probabilities(counts, alpha, beta);
        const double draw = uniform(generator);
        double cumulative = 0.0;
        int chosen = clusters - 1;
        for (int k = 0; k < clusters; ++k) {
            cumulative += probabilities[static_cast<std::size_t>(k)];
            if (draw <= cumulative) {
                chosen = k;
                break;
            }
        }
        labels[index] = chosen;
    }
}

void potts_swendsen_wang_step(std::vector<int>& labels, std::size_t rows, std::size_t cols,
                              const std::vector<double>& alpha, const std::vector<double>& beta,
                              std::mt19937& generator) {
    require(labels.size() == rows * cols && alpha.size() == beta.size(),
            "potts_swendsen_wang_step: size");
    const std::size_t n = labels.size();
    const int clusters = static_cast<int>(alpha.size());
    std::vector<std::size_t> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (std::size_t row = 0; row < rows; ++row)
        for (std::size_t column = 0; column < cols; ++column) {
            const std::size_t index = row * cols + column;
            if (column + 1 < cols) {
                const std::size_t right = index + 1;
                if (labels[index] == labels[right] && labels[index] >= 0 &&
                    labels[index] < clusters) {
                    const double probability =
                        1.0 - std::exp(-beta[static_cast<std::size_t>(labels[index])]);
                    if (uniform(generator) < probability) unite(parent, index, right);
                }
            }
            if (row + 1 < rows) {
                const std::size_t down = index + cols;
                if (labels[index] == labels[down] && labels[index] >= 0 &&
                    labels[index] < clusters) {
                    const double probability =
                        1.0 - std::exp(-beta[static_cast<std::size_t>(labels[index])]);
                    if (uniform(generator) < probability) unite(parent, index, down);
                }
            }
        }
    std::vector<int> component_label(n, -1);
    std::vector<double> weights(static_cast<std::size_t>(clusters));
    double normalisation = 0.0;
    for (int k = 0; k < clusters; ++k) {
        weights[static_cast<std::size_t>(k)] = std::exp(alpha[static_cast<std::size_t>(k)]);
        normalisation += weights[static_cast<std::size_t>(k)];
    }
    for (std::size_t index = 0; index < n; ++index) {
        const std::size_t root = find_root(parent, index);
        if (component_label[root] >= 0) {
            labels[index] = component_label[root];
            continue;
        }
        const double draw = uniform(generator) * normalisation;
        double cumulative = 0.0;
        int chosen = clusters - 1;
        for (int k = 0; k < clusters; ++k) {
            cumulative += weights[static_cast<std::size_t>(k)];
            if (draw <= cumulative) {
                chosen = k;
                break;
            }
        }
        component_label[root] = chosen;
        labels[index] = chosen;
    }
}

double quantile_influence(double time, double mean_time, double time_sum_squares,
                          double probability, bool above_quantile, double conditional_density) {
    require(time_sum_squares > 0.0 && conditional_density > 0.0,
            "quantile_influence: inputs");
    const double indicator = above_quantile ? 0.0 : 1.0;
    return (time - mean_time) / time_sum_squares *
           (probability - indicator) / conditional_density;
}

double fused_trend(const std::vector<double>& weights, const std::vector<double>& trends) {
    require(!weights.empty() && weights.size() == trends.size(),
            "fused_trend: size");
    double weighted = 0.0, total = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) {
        weighted += weights[i] * trends[i];
        total += weights[i];
    }
    require(total > 0.0, "fused_trend: zero weight");
    return weighted / total;
}

double fused_weight_sensitivity(const std::vector<double>& trends, double fused,
                                std::size_t source) {
    require(source < trends.size(), "fused_weight_sensitivity: source");
    return trends[source] - fused;
}

double vertical_trend_difference(double free_troposphere, double boundary_layer) {
    return free_troposphere - boundary_layer;
}

double trend_signal_to_noise(double trend, double standard_error) {
    require(standard_error > 0.0, "trend_signal_to_noise: standard error");
    return trend / standard_error;
}

double potts_energy(const std::vector<int>& labels, std::size_t rows, std::size_t cols,
                    const std::vector<double>& alpha, const std::vector<double>& beta) {
    require(labels.size() == rows * cols && alpha.size() == beta.size(),
            "potts_energy: size");
    const std::size_t clusters = alpha.size();
    std::vector<double> counts(clusters, 0.0);
    for (int label : labels) {
        require(label >= 0 && static_cast<std::size_t>(label) < clusters,
                "potts_energy: label");
        counts[static_cast<std::size_t>(label)] += 1.0;
    }
    double energy = 0.0;
    for (std::size_t k = 0; k < clusters; ++k) energy += alpha[k] * counts[k];
    for (std::size_t row = 0; row < rows; ++row)
        for (std::size_t column = 0; column < cols; ++column) {
            const std::size_t index = row * cols + column;
            const int label = labels[index];
            if (column + 1 < cols && labels[index + 1] == label)
                energy += beta[static_cast<std::size_t>(label)];
            if (row + 1 < rows && labels[index + cols] == label)
                energy += beta[static_cast<std::size_t>(label)];
        }
    return energy;
}

double parallel_tempering_swap_log_ratio(double energy_cold, double energy_hot,
                                         double cold_temperature, double hot_temperature) {
    require(cold_temperature > 0.0 && hot_temperature > 0.0,
            "parallel_tempering_swap_log_ratio: temperatures");
    return (1.0 / cold_temperature - 1.0 / hot_temperature) *
           (energy_hot - energy_cold);
}

bool potts_exchange_update(
    std::vector<int>& labels, std::size_t rows, std::size_t cols,
    std::vector<double>& alpha, std::vector<double>& beta,
    const std::function<double(const std::vector<int>&)>& log_likelihood,
    double proposal_sd, std::mt19937& generator) {
    require(labels.size() == rows * cols && alpha.size() == beta.size() && proposal_sd >= 0.0,
            "potts_exchange_update: inputs");
    require(alpha.size() == 1,
            "potts_exchange_update: single cluster exchange is degenerate");
    std::normal_distribution<double> normal(0.0, proposal_sd);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::vector<double> proposed_alpha = alpha, proposed_beta = beta;
    proposed_alpha[0] += normal(generator);
    proposed_beta[0] = std::max(0.0, proposed_beta[0] + normal(generator));
    std::vector<int> auxiliary = labels;
    for (int sweep = 0; sweep < 10; ++sweep)
        potts_swendsen_wang_step(auxiliary, rows, cols, proposed_alpha, proposed_beta,
                                 generator);
    const double log_ratio =
        potts_energy(auxiliary, rows, cols, proposed_alpha, proposed_beta) -
        potts_energy(labels, rows, cols, alpha, beta) +
        log_likelihood(auxiliary) - log_likelihood(labels);
    if (uniform(generator) < std::exp(std::min(0.0, log_ratio))) {
        labels = auxiliary;
        alpha = proposed_alpha;
        beta = proposed_beta;
        return true;
    }
    return false;
}

}  // namespace climate
