#include "multitrend.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

#include "linalg.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double sign(double value) { return value > 0.0 ? 1.0 : (value < 0.0 ? -1.0 : 0.0); }

std::vector<double> time_order(std::size_t n) {
    std::vector<double> times(n);
    std::iota(times.begin(), times.end(), 1.0);
    return times;
}

Matrix invert_matrix(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    Matrix inverse_matrix(n, std::vector<double>(n, 0.0));
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<double> unit(n, 0.0), solution;
        unit[column] = 1.0;
        if (!solve_linear_system(matrix, unit, solution))
            throw std::runtime_error("invert_matrix: singular matrix");
        for (std::size_t row = 0; row < n; ++row) inverse_matrix[row][column] = solution[row];
    }
    return inverse_matrix;
}

double kendall_tau(const std::vector<double>& first, const std::vector<double>& second) {
    require(first.size() == second.size() && first.size() >= 2, "kendall_tau: size");
    double sum = 0.0;
    for (std::size_t i = 0; i < first.size(); ++i)
        for (std::size_t j = i + 1; j < first.size(); ++j)
            sum += sign(first[j] - first[i]) * sign(second[j] - second[i]);
    const double pairs = 0.5 * static_cast<double>(first.size()) *
                         static_cast<double>(first.size() - 1);
    return sum / pairs;
}

}  // namespace

double mann_kendall_score(const std::vector<double>& values) {
    double score = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i)
        for (std::size_t j = i + 1; j < values.size(); ++j)
            score += sign(values[j] - values[i]);
    return score;
}

double mann_kendall_variance(std::size_t n) {
    require(n >= 2, "mann_kendall_variance: n");
    return static_cast<double>(n) * static_cast<double>(n - 1) *
           static_cast<double>(2 * n + 5) / 18.0;
}

double mann_kendall_rank(const std::vector<double>& values, std::size_t index) {
    require(index < values.size(), "mann_kendall_rank: index");
    double sum = 0.0;
    for (double value : values) sum += sign(values[index] - value);
    return 0.5 * (static_cast<double>(values.size()) + 1.0 + sum);
}

double mann_kendall_concordance(const std::vector<double>& first,
                                const std::vector<double>& second) {
    require(first.size() == second.size(), "mann_kendall_concordance: size");
    double sum = 0.0;
    for (std::size_t i = 0; i < first.size(); ++i)
        for (std::size_t j = i + 1; j < first.size(); ++j)
            sum += sign((first[j] - first[i]) * (second[j] - second[i]));
    return sum;
}

double mann_kendall_cross_covariance(const std::vector<double>& first,
                                     const std::vector<double>& second) {
    require(first.size() == second.size() && first.size() >= 2,
            "mann_kendall_cross_covariance: size");
    const std::size_t n = first.size();
    const double concordance = mann_kendall_concordance(first, second);
    double rank_product = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        rank_product += mann_kendall_rank(first, i) * mann_kendall_rank(second, i);
    const double total = static_cast<double>(n) * std::pow(static_cast<double>(n) + 1.0, 2.0);
    return (concordance + 4.0 * rank_product - total) / 3.0;
}

double multivariate_mann_kendall(const std::vector<double>& scores,
                                 const Matrix& covariance) {
    require(!scores.empty() && covariance.size() == scores.size(),
            "multivariate_mann_kendall: size");
    double total = 0.0;
    for (double score : scores) total += score;
    double variance = 0.0;
    for (std::size_t i = 0; i < scores.size(); ++i)
        for (std::size_t j = 0; j < scores.size(); ++j) variance += covariance[i][j];
    require(variance > 0.0, "multivariate_mann_kendall: zero variance");
    return total / std::sqrt(variance);
}

double multivariate_mann_kendall_quadratic(const std::vector<double>& scores,
                                           const Matrix& covariance) {
    require(!scores.empty() && covariance.size() == scores.size(),
            "multivariate_mann_kendall_quadratic: size");
    const Matrix inverse_covariance = invert_matrix(covariance);
    double quadratic = 0.0;
    for (std::size_t i = 0; i < scores.size(); ++i)
        for (std::size_t j = 0; j < scores.size(); ++j)
            quadratic += scores[i] * inverse_covariance[i][j] * scores[j];
    return quadratic;
}

double partial_mann_kendall(const std::vector<double>& target,
                            const std::vector<double>& covariate) {
    require(target.size() == covariate.size() && target.size() >= 2,
            "partial_mann_kendall: size");
    const double target_score = mann_kendall_score(target);
    const double covariate_score = mann_kendall_score(covariate);
    const double variance = mann_kendall_variance(target.size());
    const double covariance = mann_kendall_cross_covariance(target, covariate);
    const double correlation = covariance / variance;
    if (std::fabs(correlation) >= 1.0 - 1e-12) return 0.0;
    return (target_score - correlation * covariate_score) /
           std::sqrt((1.0 - correlation * correlation) * variance);
}

std::vector<double> seasonal_kendall_scores(const Matrix& seasonal_data) {
    require(!seasonal_data.empty(), "seasonal_kendall_scores: empty");
    std::vector<double> scores(seasonal_data.size(), 0.0);
    for (std::size_t season = 0; season < seasonal_data.size(); ++season)
        scores[season] = mann_kendall_score(seasonal_data[season]);
    return scores;
}

double seasonal_kendall_variance(const Matrix& seasonal_data) {
    require(!seasonal_data.empty() && !seasonal_data[0].empty(),
            "seasonal_kendall_variance: empty");
    double variance = 0.0;
    for (const auto& season : seasonal_data) variance += mann_kendall_variance(season.size());
    return variance;
}

double seasonal_kendall_cross_covariance(const Matrix& first, const Matrix& second,
                                         bool missing_values) {
    require(first.size() == second.size() && !first.empty(), "seasonal cov: size");
    double covariance = 0.0;
    for (std::size_t season = 0; season < first.size(); ++season) {
        const std::size_t n = first[season].size();
        require(second[season].size() == n, "seasonal cov: season size");
        const double concordance = mann_kendall_concordance(first[season], second[season]);
        double rank_product = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            rank_product += mann_kendall_rank(first[season], i) *
                            mann_kendall_rank(second[season], i);
        const double baseline = missing_values
                                    ? static_cast<double>(n) *
                                          (static_cast<double>(n) + 1.0) *
                                          (static_cast<double>(n) + 1.0)
                                    : static_cast<double>(n) *
                                          std::pow(static_cast<double>(n) + 1.0, 2.0);
        covariance += (concordance + 4.0 * rank_product - baseline) / 3.0;
    }
    return covariance;
}

std::vector<double> moving_window_kendall(const std::vector<double>& first,
                                          const std::vector<double>& second,
                                          std::size_t window) {
    require(first.size() == second.size() && window >= 2 && window <= first.size(),
            "moving_window_kendall: inputs");
    const std::size_t windows = first.size() - window + 1;
    std::vector<double> sequence(windows, 0.0);
    for (std::size_t start = 0; start < windows; ++start) {
        std::vector<double> x(first.begin() + static_cast<long>(start),
                              first.begin() + static_cast<long>(start + window));
        std::vector<double> y(second.begin() + static_cast<long>(start),
                              second.begin() + static_cast<long>(start + window));
        sequence[start] = kendall_tau(x, y);
    }
    return sequence;
}

double dependence_trend(const std::vector<double>& first, const std::vector<double>& second,
                        std::size_t window) {
    const std::vector<double> sequence = moving_window_kendall(first, second, window);
    return kendall_tau(sequence, time_order(sequence.size()));
}

double overall_multivariate_trend(const std::vector<double>& first,
                                  const std::vector<double>& second, std::size_t window) {
    require(first.size() == second.size() && first.size() >= 2,
            "overall_multivariate_trend: size");
    const std::vector<double> times = time_order(first.size());
    const double first_trend = kendall_tau(first, times);
    const double second_trend = kendall_tau(second, times);
    const double dependence = dependence_trend(first, second, window);
    return (first_trend * first_trend + second_trend * second_trend +
            dependence * dependence) / 3.0;
}

double dependence_trend_p_value(const std::vector<double>& first,
                                const std::vector<double>& second, std::size_t window,
                                std::size_t resamples, unsigned seed) {
    require(first.size() == second.size() && resamples > 0,
            "dependence_trend_p_value: inputs");
    const std::vector<double> sequence = moving_window_kendall(first, second, window);
    const std::vector<double> times = time_order(sequence.size());
    const double observed = std::fabs(kendall_tau(sequence, times));
    std::vector<double> shuffled = times;
    std::mt19937 generator(seed);
    std::size_t exceedances = 0;
    for (std::size_t resample = 0; resample < resamples; ++resample) {
        std::shuffle(shuffled.begin(), shuffled.end(), generator);
        if (std::fabs(kendall_tau(sequence, shuffled)) >= observed) ++exceedances;
    }
    return static_cast<double>(exceedances) / static_cast<double>(resamples);
}

}  // namespace climate
