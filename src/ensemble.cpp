#include "ensemble.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

#include "lagrangian.hpp"
#include "regional.hpp"
#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

Matrix spatial_correlation_matrix(std::size_t rows, std::size_t columns,
                                  double neighbour_correlation) {
    require(rows > 0 && columns > 0, "spatial_correlation_matrix: dimensions");
    require(neighbour_correlation > 0.0 && neighbour_correlation < 1.0,
            "spatial_correlation_matrix: correlation");
    const std::size_t n = rows * columns;
    Matrix correlation(n, std::vector<double>(n, 1.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const double dr = static_cast<double>(i / columns) -
                              static_cast<double>(j / columns);
            const double dc = static_cast<double>(i % columns) -
                              static_cast<double>(j % columns);
            correlation[i][j] = std::pow(neighbour_correlation, std::hypot(dr, dc));
        }
    return correlation;
}

bool cholesky_factor(const Matrix& matrix, Matrix& factor) {
    if (matrix.empty() || matrix.size() != matrix[0].size()) return false;
    const std::size_t n = matrix.size();
    factor.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k) sum -= factor[i][k] * factor[j][k];
            if (i == j) {
                if (sum <= 1e-12) return false;
                factor[i][j] = std::sqrt(sum);
            } else {
                factor[i][j] = sum / factor[j][j];
            }
        }
    }
    return true;
}

double empirical_sign_test_p_value(std::size_t positive, std::size_t negative,
                                   std::size_t rows, std::size_t columns,
                                   std::size_t years, double spatial_correlation,
                                   double temporal_correlation, std::size_t realizations,
                                   unsigned seed) {
    require(rows > 0 && columns > 0 && years > 2 && realizations > 0,
            "empirical_sign_test_p_value: inputs");
    require(positive + negative > 0, "empirical_sign_test_p_value: no sites");
    const std::size_t observed = std::max(positive, negative);
    const Matrix correlation = spatial_correlation_matrix(rows, columns, spatial_correlation);
    Matrix factor;
    require(cholesky_factor(correlation, factor),
            "empirical_sign_test_p_value: invalid correlation");
    const std::size_t sites = correlation.size();
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::size_t exceedances = 0;
    Matrix field(years, std::vector<double>(sites, 0.0));
    std::vector<double> innovation(sites, 0.0), persistent(sites, 0.0);
    for (std::size_t realization = 0; realization < realizations; ++realization) {
        for (std::size_t t = 0; t < years; ++t) {
            for (double& value : innovation) value = normal(generator);
            for (std::size_t site = 0; site < sites; ++site) {
                double sum = 0.0;
                for (std::size_t k = 0; k <= site; ++k) sum += factor[site][k] * innovation[k];
                persistent[site] = t == 0 ? sum
                                          : temporal_correlation * persistent[site] +
                                                std::sqrt(std::max(0.0, 1.0 -
                                                    temporal_correlation * temporal_correlation)) * sum;
                field[t][site] = persistent[site];
            }
        }
        detrend_field(field);
        std::size_t sign_positive = 0;
        for (std::size_t site = 0; site < sites; ++site) {
            std::vector<double> series(years, 0.0);
            for (std::size_t t = 0; t < years; ++t) series[t] = field[t][site];
            if (linear_trend(series) > 0.0) ++sign_positive;
        }
        if (std::max(sign_positive, sites - sign_positive) >= observed) ++exceedances;
    }
    return static_cast<double>(exceedances) / static_cast<double>(realizations);
}

double wilks_correlation(double distance, double d1, double d2) {
    require(distance >= 0.0, "wilks_correlation: distance");
    return std::exp(d1 * std::pow(distance, d2));
}

Matrix simulate_occurrence_field(std::size_t rows, std::size_t columns, std::size_t days,
                                 double neighbour_correlation, double occurrence_probability,
                                 unsigned seed) {
    require(occurrence_probability >= 0.0 && occurrence_probability <= 1.0,
            "simulate_occurrence_field: probability");
    const Matrix correlation =
        spatial_correlation_matrix(rows, columns, neighbour_correlation);
    Matrix factor;
    require(cholesky_factor(correlation, factor),
            "simulate_occurrence_field: invalid correlation matrix");
    const std::size_t n = correlation.size();
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    Matrix occurrence(days, std::vector<double>(n, 0.0));
    const double threshold = normal_quantile(occurrence_probability);
    std::vector<double> draw(n, 0.0), latent(n, 0.0);
    for (std::size_t t = 0; t < days; ++t) {
        for (double& value : draw) value = normal(generator);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = 0; k <= i; ++k) sum += factor[i][k] * draw[k];
            latent[i] = sum;
            occurrence[t][i] = latent[i] <= threshold ? 1.0 : 0.0;
        }
    }
    return occurrence;
}

Matrix simulate_gamma_field(std::size_t rows, std::size_t columns, double shape,
                            double rate, unsigned seed) {
    require(shape > 0.0 && rate > 0.0, "simulate_gamma_field: parameters");
    std::mt19937 generator(seed);
    std::gamma_distribution<double> gamma(shape, 1.0 / rate);
    Matrix field(rows, std::vector<double>(columns, 0.0));
    for (auto& row : field)
        for (double& value : row) value = gamma(generator);
    return field;
}

void detrend_field(Matrix& field) {
    require(!field.empty() && !field[0].empty(), "detrend_field: empty field");
    const std::size_t sites = field[0].size();
    for (std::size_t site = 0; site < sites; ++site) {
        std::vector<double> series(field.size());
        for (std::size_t t = 0; t < field.size(); ++t) series[t] = field[t][site];
        const double slope = linear_trend(series);
        for (std::size_t t = 0; t < field.size(); ++t)
            field[t][site] -= slope * static_cast<double>(t);
    }
}

Matrix simulate_zero_trend_field(std::size_t rows, std::size_t columns, std::size_t years,
                                 double neighbour_correlation, double temporal_correlation,
                                 unsigned seed) {
    require(years > 0 && temporal_correlation > -1.0 && temporal_correlation < 1.0,
            "simulate_zero_trend_field: inputs");
    const Matrix correlation =
        spatial_correlation_matrix(rows, columns, neighbour_correlation);
    Matrix factor;
    require(cholesky_factor(correlation, factor),
            "simulate_zero_trend_field: correlation matrix is not positive definite");
    const std::size_t n = correlation.size();
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    Matrix field(years, std::vector<double>(n, 0.0));
    std::vector<double> previous(n, 0.0), innovation(n, 0.0), draw(n, 0.0);
    const double memory = std::sqrt(1.0 - temporal_correlation * temporal_correlation);
    for (std::size_t t = 0; t < years; ++t) {
        for (double& value : draw) value = normal(generator);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = 0; k <= i; ++k) sum += factor[i][k] * draw[k];
            innovation[i] = sum;
        }
        for (std::size_t i = 0; i < n; ++i) {
            previous[i] = temporal_correlation * previous[i] + memory * innovation[i];
            field[t][i] = previous[i];
        }
    }
    return field;
}

std::size_t positive_trend_count(const Matrix& field) {
    require(!field.empty() && !field[0].empty(), "positive_trend_count: empty field");
    const std::size_t sites = field[0].size();
    std::size_t count = 0;
    for (std::size_t s = 0; s < sites; ++s) {
        std::vector<double> series(field.size());
        for (std::size_t t = 0; t < field.size(); ++t) series[t] = field[t][s];
        if (linear_trend(series) > 0.0) ++count;
    }
    return count;
}

std::vector<std::size_t> sign_test_null_ensemble(std::size_t rows, std::size_t columns,
                                                 std::size_t years,
                                                 double neighbour_correlation,
                                                 double temporal_correlation, int replicates,
                                                 unsigned seed) {
    require(replicates > 0, "sign_test_null_ensemble: replicates");
    std::vector<std::size_t> counts(static_cast<std::size_t>(replicates), 0);
    for (int r = 0; r < replicates; ++r) {
        const Matrix field = simulate_zero_trend_field(
            rows, columns, years, neighbour_correlation, temporal_correlation,
            seed + static_cast<unsigned>(r) * 7919u);
        counts[static_cast<std::size_t>(r)] = positive_trend_count(field);
    }
    return counts;
}

void inject_uniform_trend(Matrix& field, double trend) {
    for (std::size_t t = 0; t < field.size(); ++t)
        for (double& value : field[t]) value += trend * static_cast<double>(t);
}

double sign_test_power(std::size_t rows, std::size_t columns, std::size_t years,
                       double neighbour_correlation, double temporal_correlation,
                       double trend, double alpha, int null_replicates,
                       int signal_replicates, unsigned seed) {
    require(null_replicates > 0 && signal_replicates > 0, "sign_test_power: replicates");
    const std::vector<std::size_t> null_counts = sign_test_null_ensemble(
        rows, columns, years, neighbour_correlation, temporal_correlation, null_replicates,
        seed);
    const double critical = sign_test_critical_value(
        static_cast<double>(rows * columns), neighbour_correlation, alpha);
    const double threshold = 0.5 * (1.0 + critical) *
                             static_cast<double>(rows * columns);
    int rejections = 0;
    for (int r = 0; r < signal_replicates; ++r) {
        Matrix field = simulate_zero_trend_field(
            rows, columns, years, neighbour_correlation, temporal_correlation,
            seed + 104729u + static_cast<unsigned>(r) * 7919u);
        inject_uniform_trend(field, trend);
        if (static_cast<double>(positive_trend_count(field)) >= threshold) ++rejections;
    }
    return static_cast<double>(rejections) / static_cast<double>(signal_replicates);
}

double central_difference(const std::function<double(double)>& function, double point,
                          double step) {
    require(step > 0.0, "central_difference: step");
    return (function(point + step) - function(point - step)) / (2.0 * step);
}

}  // namespace climate
