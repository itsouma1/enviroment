#include "relevance.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "linalg.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double r_squared(const std::vector<double>& response, const Matrix& predictors) {
    if (predictors.empty() || predictors[0].empty()) return 0.0;
    require(predictors.size() == response.size(),
            "r_squared: response and predictor rows must match");
    const std::size_t n = response.size(), columns = predictors[0].size() + 1;
    Matrix design(n, std::vector<double>(columns, 1.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < predictors[i].size(); ++j)
            design[i][j + 1] = predictors[i][j];
    std::vector<double> coefficients;
    if (!least_squares(design, response, coefficients)) return 0.0;
    double mean = 0.0;
    for (double value : response) mean += value;
    mean /= static_cast<double>(n);
    double total = 0.0, residual = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double fitted = 0.0;
        for (std::size_t j = 0; j < coefficients.size(); ++j)
            fitted += coefficients[j] * design[i][j];
        total += (response[i] - mean) * (response[i] - mean);
        residual += (response[i] - fitted) * (response[i] - fitted);
    }
    if (total <= 0.0) return 0.0;
    return std::max(0.0, 1.0 - residual / total);
}

Matrix columns(const Matrix& matrix, const std::vector<std::size_t>& selected) {
    Matrix result(matrix.size(), std::vector<double>(selected.size(), 0.0));
    for (std::size_t i = 0; i < matrix.size(); ++i)
        for (std::size_t j = 0; j < selected.size(); ++j)
            result[i][j] = matrix[i][selected[j]];
    return result;
}

}  // namespace

double variance_explained(const std::vector<double>& response, const Matrix& predictors) {
    return r_squared(response, predictors);
}

double gaussian_mutual_information(double explained) {
    require(explained >= 0.0 && explained < 1.0, "gaussian_mutual_information: range");
    return -0.5 * std::log(1.0 - explained);
}

double relevance_score(double variance_explained, double mutual_information,
                       double penalty, double lambda_g, double lambda_i,
                       double lambda_penalty) {
    return lambda_g * variance_explained * lambda_i * mutual_information -
           lambda_penalty * penalty;
}

double effective_relevance_ratio(double lambda_g, double lambda_i, double lambda_penalty) {
    require(lambda_penalty > 0.0, "effective_relevance_ratio: penalty");
    return lambda_g * lambda_i / lambda_penalty;
}

double marginal_relevance_gain(double variance_explained, double mutual_information,
                              double delta_variance, double delta_information,
                              double cost, double lambda_g, double lambda_i,
                              double lambda_penalty) {
    return lambda_g * lambda_i *
               (variance_explained * delta_information +
                mutual_information * delta_variance +
                delta_variance * delta_information) -
           lambda_penalty * cost;
}

double aic_penalty(std::size_t variables) { return static_cast<double>(variables); }

double bic_penalty(std::size_t variables, std::size_t samples) {
    require(samples > 1, "bic_penalty: samples");
    return 0.5 * static_cast<double>(variables) * std::log(static_cast<double>(samples));
}

std::vector<double> variance_inflation_factors(const Matrix& predictors) {
    require(!predictors.empty() && predictors[0].size() > 0,
            "variance_inflation_factors: empty");
    const std::size_t variables = predictors[0].size();
    std::vector<double> factors(variables, 0.0);
    for (std::size_t j = 0; j < variables; ++j) {
        std::vector<double> response(predictors.size(), 0.0);
        std::vector<std::size_t> others;
        for (std::size_t k = 0; k < variables; ++k)
            if (k != j) others.push_back(k);
        for (std::size_t i = 0; i < predictors.size(); ++i) response[i] = predictors[i][j];
        const double explained = r_squared(response, columns(predictors, others));
        factors[j] = 1.0 / std::max(1e-12, 1.0 - explained);
    }
    return factors;
}

std::vector<std::size_t> greedy_relevance_selection(const std::vector<double>& response,
                                                    const Matrix& predictors, double kappa,
                                                    std::size_t penalty_form) {
    require(!predictors.empty() && predictors.size() == response.size() && kappa >= 0.0,
            "greedy_relevance_selection: inputs");
    std::vector<std::size_t> selected;
    std::vector<std::size_t> remaining(predictors[0].size());
    for (std::size_t i = 0; i < remaining.size(); ++i) remaining[i] = i;
    double current_g = 0.0, current_i = 0.0;
    while (!remaining.empty()) {
        double best_gain = -1e300;
        std::size_t best = 0, best_index = 0;
        for (std::size_t k = 0; k < remaining.size(); ++k) {
            std::vector<std::size_t> candidate = selected;
            candidate.push_back(remaining[k]);
            const double g = r_squared(response, columns(predictors, candidate));
            const double information = g <= 0.0 ? 0.0 : gaussian_mutual_information(std::min(g, 0.999999));
            const double dg = g - current_g, di = information - current_i;
            const double penalty = penalty_form == 1
                                       ? bic_penalty(candidate.size(), response.size())
                                       : aic_penalty(candidate.size());
            const double gain = kappa * (current_g * di + current_i * dg + dg * di) - penalty;
            if (gain > best_gain) {
                best_gain = gain;
                best = remaining[k];
                best_index = k;
            }
        }
        if (best_gain <= 0.0) break;
        selected.push_back(best);
        remaining.erase(remaining.begin() + static_cast<long>(best_index));
        current_g = r_squared(response, columns(predictors, selected));
        current_i = current_g <= 0.0 ? 0.0 : gaussian_mutual_information(std::min(current_g, 0.999999));
    }
    return selected;
}

}  // namespace climate
