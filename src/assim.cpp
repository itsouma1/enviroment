#include "assim.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

using Matrix = std::vector<std::vector<double>>;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Matrix identity(std::size_t n) {
    Matrix result(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) result[i][i] = 1.0;
    return result;
}

Matrix transpose(const Matrix& matrix) {
    if (matrix.empty()) return {};
    Matrix result(matrix[0].size(), std::vector<double>(matrix.size(), 0.0));
    for (std::size_t i = 0; i < matrix.size(); ++i)
        for (std::size_t j = 0; j < matrix[i].size(); ++j) result[j][i] = matrix[i][j];
    return result;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
    require(!a.empty() && !b.empty(), "multiply: empty matrix");
    require(a[0].size() == b.size(), "multiply: shape mismatch");
    Matrix result(a.size(), std::vector<double>(b[0].size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t k = 0; k < b.size(); ++k)
            for (std::size_t j = 0; j < b[0].size(); ++j)
                result[i][j] += a[i][k] * b[k][j];
    return result;
}

std::vector<double> multiply(const Matrix& a, const std::vector<double>& x) {
    require(a.size() == x.size(), "multiply: vector shape mismatch");
    std::vector<double> result(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < x.size(); ++j) result[i] += a[i][j] * x[j];
    return result;
}

Matrix cholesky_lower(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    Matrix lower(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        if (matrix[i].size() != n) throw std::runtime_error("assim: matrix not square");
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k) sum -= lower[i][k] * lower[j][k];
            if (i == j) {
                require(sum > 0.0, "assim: matrix not positive definite");
                lower[i][j] = std::sqrt(sum);
            } else {
                lower[i][j] = sum / lower[j][j];
            }
        }
    }
    return lower;
}

std::vector<double> forward_substitute(const Matrix& lower, const std::vector<double>& rhs) {
    const std::size_t n = lower.size();
    std::vector<double> solution(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = rhs[i];
        for (std::size_t k = 0; k < i; ++k) sum -= lower[i][k] * solution[k];
        solution[i] = sum / lower[i][i];
    }
    return solution;
}


Matrix inverse(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    Matrix result(n, std::vector<double>(n, 0.0));
    // Gauss-Jordan with partial pivoting.
    Matrix a = matrix;
    Matrix b = identity(n);
    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < n; ++row)
            if (std::fabs(a[row][column]) > std::fabs(a[pivot][column])) pivot = row;
        require(std::fabs(a[pivot][column]) > 0.0, "assim: singular matrix");
        std::swap(a[column], a[pivot]);
        std::swap(b[column], b[pivot]);
        const double scale = a[column][column];
        for (std::size_t j = 0; j < n; ++j) {
            a[column][j] /= scale;
            b[column][j] /= scale;
        }
        for (std::size_t row = 0; row < n; ++row) {
            if (row == column) continue;
            const double factor = a[row][column];
            for (std::size_t j = 0; j < n; ++j) {
                a[row][j] -= factor * a[column][j];
                b[row][j] -= factor * b[column][j];
            }
        }
    }
    return b;
}

Matrix covariance_of(const std::vector<std::vector<double>>& members,
                     const std::vector<double>& mean) {
    const std::size_t m = members.size();
    const std::size_t n = mean.size();
    require(m >= 2, "assim: at least two members required");
    Matrix covariance(n, std::vector<double>(n, 0.0));
    for (const auto& member : members) {
        require(member.size() == n, "assim: member size mismatch");
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                covariance[i][j] += (member[i] - mean[i]) * (member[j] - mean[j]);
    }
    for (auto& row : covariance)
        for (double& value : row) value /= static_cast<double>(m - 1);
    return covariance;
}

}  // namespace

std::vector<double> hybrid_analysis(const std::vector<double>& enkf_analysis,
                                    const std::vector<double>& enviar_analysis, double gamma) {
    require(enkf_analysis.size() == enviar_analysis.size(), "hybrid_analysis: size mismatch");
    if (gamma < 0.0 || gamma > 1.0)
        throw std::runtime_error("hybrid_analysis: gamma outside [0, 1]");
    std::vector<double> analysis(enkf_analysis.size(), 0.0);
    for (std::size_t i = 0; i < analysis.size(); ++i)
        analysis[i] = gamma * enkf_analysis[i] + (1.0 - gamma) * enviar_analysis[i];
    return analysis;
}

std::vector<double> hybrid_analysis_sensitivity(const std::vector<double>& enkf_analysis,
                                                const std::vector<double>& enviar_analysis) {
    require(enkf_analysis.size() == enviar_analysis.size(),
            "hybrid_analysis_sensitivity: size mismatch");
    std::vector<double> sensitivity(enkf_analysis.size(), 0.0);
    for (std::size_t i = 0; i < sensitivity.size(); ++i)
        sensitivity[i] = enkf_analysis[i] - enviar_analysis[i];
    return sensitivity;
}

std::vector<std::vector<double>> hybrid_gain(const std::vector<std::vector<double>>& enkf_gain,
                                             const std::vector<std::vector<double>>& enviar_gain,
                                             double gamma) {
    require(enkf_gain.size() == enviar_gain.size(), "hybrid_gain: row mismatch");
    std::vector<std::vector<double>> gain(enkf_gain.size());
    for (std::size_t i = 0; i < enkf_gain.size(); ++i) {
        require(enkf_gain[i].size() == enviar_gain[i].size(), "hybrid_gain: column mismatch");
        gain[i].resize(enkf_gain[i].size());
        for (std::size_t j = 0; j < gain[i].size(); ++j)
            gain[i][j] = gamma * enkf_gain[i][j] + (1.0 - gamma) * enviar_gain[i][j];
    }
    return gain;
}

EnsembleAnalysis ensemble_adjustment_filter(
    const std::vector<std::vector<double>>& members, const std::vector<double>& observation,
    const std::vector<std::vector<double>>& observation_operator,
    const std::vector<std::vector<double>>& observation_covariance) {
    require(!members.empty(), "eakf: empty ensemble");
    const std::size_t m = observation.size();
    const std::size_t n = members[0].size();
    require(observation_operator.size() == m, "eakf: H row mismatch");
    for (const auto& row : observation_operator)
        require(row.size() == n, "eakf: H column mismatch");

    std::vector<double> mean(n, 0.0);
    for (const auto& member : members) {
        require(member.size() == n, "eakf: member size mismatch");
        for (std::size_t i = 0; i < n; ++i) mean[i] += member[i];
    }
    for (double& value : mean) value /= static_cast<double>(members.size());
    const Matrix prior_covariance = covariance_of(members, mean);

    const Matrix h = observation_operator;
    const Matrix ht = transpose(h);
    const Matrix r_inverse = inverse(observation_covariance);
    const Matrix prior_inverse = inverse(prior_covariance);
    const Matrix ht_r_inv = multiply(ht, r_inverse);
    const Matrix analysis_precision = [&] {
        Matrix result = multiply(ht_r_inv, h);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) result[i][j] += prior_inverse[i][j];
        return result;
    }();
    const Matrix posterior_covariance = inverse(analysis_precision);

    std::vector<double> rhs = multiply(prior_inverse, mean);
    const std::vector<double> ht_r_inv_y = multiply(ht_r_inv, observation);
    for (std::size_t i = 0; i < n; ++i) rhs[i] += ht_r_inv_y[i];
    const std::vector<double> posterior_mean = multiply(posterior_covariance, rhs);

    const Matrix lower_prior = cholesky_lower(prior_covariance);
    const Matrix lower_posterior = cholesky_lower(posterior_covariance);
    const Matrix lower_posterior_t = transpose(lower_posterior);

    EnsembleAnalysis result;
    result.mean = posterior_mean;
    result.covariance = posterior_covariance;
    result.members.reserve(members.size());
    for (const auto& member : members) {
        std::vector<double> deviation(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) deviation[i] = member[i] - mean[i];
        const std::vector<double> whitened = forward_substitute(lower_prior, deviation);
        std::vector<double> transformed = multiply(lower_posterior_t, whitened);
        for (std::size_t i = 0; i < n; ++i) transformed[i] += posterior_mean[i];
        result.members.push_back(std::move(transformed));
    }
    return result;
}

std::vector<double> network_forward(const FeedForwardNetwork& network,
                                    const std::vector<double>& input) {
    require(!network.weights.empty(), "network_forward: empty network");
    require(network.weights.size() == network.biases.size(), "network_forward: layer mismatch");
    std::vector<double> activations = input;
    for (std::size_t layer = 0; layer < network.weights.size(); ++layer) {
        const auto& weights = network.weights[layer];
        const auto& biases = network.biases[layer];
        require(weights.size() == biases.size(), "network_forward: bias mismatch");
        std::vector<double> next(weights.size(), 0.0);
        for (std::size_t out = 0; out < weights.size(); ++out) {
            require(weights[out].size() == activations.size(),
                    "network_forward: weight shape mismatch");
            double value = biases[out];
            for (std::size_t in = 0; in < activations.size(); ++in)
                value += weights[out][in] * activations[in];
            if (layer < network.activation.size() && network.activation[layer] == 1)
                value = std::max(value, 0.0);
            next[out] = value;
        }
        activations.swap(next);
    }
    return activations;
}

std::vector<std::vector<double>> network_jacobian(const FeedForwardNetwork& network,
                                                  const std::vector<double>& input) {
    require(!network.weights.empty(), "network_jacobian: empty network");
    // Forward-mode propagation of d(activation)/d(input).
    std::vector<std::vector<double>> jacobian(input.size(),
                                              std::vector<double>(input.size(), 0.0));
    for (std::size_t i = 0; i < input.size(); ++i) jacobian[i][i] = 1.0;
    std::vector<double> activations = input;
    for (std::size_t layer = 0; layer < network.weights.size(); ++layer) {
        const auto& weights = network.weights[layer];
        const auto& biases = network.biases[layer];
        require(weights.size() == biases.size(), "network_jacobian: bias mismatch");
        std::vector<std::vector<double>> next_jacobian(
            weights.size(), std::vector<double>(input.size(), 0.0));
        std::vector<double> next(weights.size(), 0.0);
        for (std::size_t out = 0; out < weights.size(); ++out) {
            require(weights[out].size() == activations.size(),
                    "network_jacobian: weight shape mismatch");
            double value = biases[out];
            for (std::size_t in = 0; in < activations.size(); ++in)
                value += weights[out][in] * activations[in];
            double slope = 1.0;
            if (layer < network.activation.size() && network.activation[layer] == 1) {
                slope = value > 0.0 ? 1.0 : 0.0;
                value = std::max(value, 0.0);
            }
            for (std::size_t in = 0; in < activations.size(); ++in)
                for (std::size_t j = 0; j < input.size(); ++j)
                    next_jacobian[out][j] += slope * weights[out][in] * jacobian[in][j];
            next[out] = value;
        }
        activations.swap(next);
        jacobian.swap(next_jacobian);
    }
    return jacobian;
}

}  // namespace climate
