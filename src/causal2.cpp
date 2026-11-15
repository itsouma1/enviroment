#include "causal2.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "linalg.hpp"
#include "multitrend.hpp"
#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<double> residuals(const std::vector<double>& values, const Matrix& predictors) {
    require(predictors.size() == values.size(), "residuals: size");
    const std::size_t n = values.size();
    const std::size_t columns = predictors[0].size() + 1;
    Matrix design(n, std::vector<double>(columns, 1.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < predictors[i].size(); ++j)
            design[i][j + 1] = predictors[i][j];
    std::vector<double> coefficients;
    require(least_squares(design, values, coefficients), "residuals: rank deficient");
    std::vector<double> result(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = values[i];
        for (std::size_t j = 0; j < coefficients.size(); ++j)
            result[i] -= coefficients[j] * design[i][j];
    }
    return result;
}

double correlation(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = x.size();
    double mean_x = 0.0, mean_y = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }
    mean_x /= static_cast<double>(n);
    mean_y /= static_cast<double>(n);
    double covariance = 0.0, variance_x = 0.0, variance_y = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        covariance += (x[i] - mean_x) * (y[i] - mean_y);
        variance_x += (x[i] - mean_x) * (x[i] - mean_x);
        variance_y += (y[i] - mean_y) * (y[i] - mean_y);
    }
    require(variance_x > 0.0 && variance_y > 0.0, "correlation: zero variance");
    return covariance / std::sqrt(variance_x * variance_y);
}

Matrix correlation_matrix(const Matrix& field) {
    const std::size_t sites = field[0].size();
    Matrix result(sites, std::vector<double>(sites, 0.0));
    std::vector<double> mean(sites, 0.0);
    for (const auto& row : field)
        for (std::size_t s = 0; s < sites; ++s) mean[s] += row[s];
    for (double& value : mean) value /= static_cast<double>(field.size());
    for (std::size_t i = 0; i < sites; ++i)
        for (std::size_t j = 0; j < sites; ++j) {
            double covariance = 0.0, variance_i = 0.0, variance_j = 0.0;
            for (std::size_t t = 0; t < field.size(); ++t) {
                covariance += (field[t][i] - mean[i]) * (field[t][j] - mean[j]);
                variance_i += (field[t][i] - mean[i]) * (field[t][i] - mean[i]);
                variance_j += (field[t][j] - mean[j]) * (field[t][j] - mean[j]);
            }
            result[i][j] = covariance / std::sqrt(std::max(1e-300, variance_i * variance_j));
        }
    return result;
}

double determinant(const Matrix& matrix) {
    Matrix a = matrix;
    const std::size_t n = a.size();
    double det = 1.0;
    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < n; ++row)
            if (std::fabs(a[row][column]) > std::fabs(a[pivot][column])) pivot = row;
        if (std::fabs(a[pivot][column]) < 1e-15) return 0.0;
        if (pivot != column) {
            std::swap(a[pivot], a[column]);
            det = -det;
        }
        det *= a[column][column];
        for (std::size_t row = column + 1; row < n; ++row) {
            const double factor = a[row][column] / a[column][column];
            for (std::size_t k = column; k < n; ++k) a[row][k] -= factor * a[column][k];
        }
    }
    return det;
}

}  // namespace

double partial_correlation(const std::vector<double>& x, const std::vector<double>& y,
                           const Matrix& conditioning) {
    require(x.size() == y.size(), "partial_correlation: size");
    if (conditioning.empty() || conditioning[0].empty()) return correlation(x, y);
    require(conditioning.size() == x.size(), "partial_correlation: conditioning size");
    return correlation(residuals(x, conditioning), residuals(y, conditioning));
}

double fisher_z_p_value(double correlation, std::size_t samples) {
    require(samples > 3 && std::fabs(correlation) < 1.0, "fisher_z_p_value: inputs");
    const double statistic =
        std::fabs(0.5 * std::log((1.0 + correlation) / (1.0 - correlation))) *
        std::sqrt(static_cast<double>(samples - 3));
    return 2.0 * (1.0 - normal_cdf(statistic));
}

std::vector<CausalLink> pcmci_lagged_links(const Matrix& series, std::size_t max_lag,
                                           double alpha) {
    require(!series.empty() && series[0].size() > 0 && max_lag > 0 &&
                max_lag < series.size() && alpha > 0.0 && alpha < 1.0,
            "pcmci_lagged_links: inputs");
    const std::size_t n = series.size(), variables = series[0].size();
    std::vector<CausalLink> links;
    for (std::size_t target = 0; target < variables; ++target) {
        for (std::size_t driver = 0; driver < variables; ++driver) {
            for (std::size_t lag = 1; lag <= max_lag; ++lag) {
                const std::size_t rows = n - lag;
                if (rows < 10) continue;
                std::vector<double> driver_series(rows, 0.0), target_series(rows, 0.0);
                for (std::size_t t = 0; t < rows; ++t) {
                    driver_series[t] = series[t][driver];
                    target_series[t] = series[t + lag][target];
                }
                std::vector<std::pair<double, std::size_t>> ranked;
                Matrix candidates(rows, std::vector<double>(0, 0.0));
                std::vector<std::vector<double>> columns;
                for (std::size_t other = 0; other < variables; ++other)
                    for (std::size_t other_lag = 1; other_lag <= max_lag; ++other_lag) {
                        if (other == driver && other_lag == lag) continue;
                        const std::size_t candidate_rows = n - other_lag;
                        if (candidate_rows < rows) continue;
                        std::vector<double> column(rows, 0.0);
                        for (std::size_t t = 0; t < rows; ++t)
                            column[t] = series[t][other];
                        const double value = std::fabs(correlation(column, target_series));
                        ranked.push_back({value, columns.size()});
                        columns.push_back(std::move(column));
                    }
                std::sort(ranked.begin(), ranked.end(),
                          [](const auto& a, const auto& b) { return a.first > b.first; });
                const std::size_t keep = std::min<std::size_t>(4, ranked.size());
                Matrix conditioning(rows, std::vector<double>(keep, 0.0));
                for (std::size_t k = 0; k < keep; ++k)
                    for (std::size_t t = 0; t < rows; ++t)
                        conditioning[t][k] = columns[ranked[k].second][t];
                const double r = partial_correlation(driver_series, target_series, conditioning);
                if (fisher_z_p_value(r, rows) < alpha)
                    links.push_back({driver, target, lag, r});
            }
        }
    }
    return links;
}

double mci_partial_correlation(double coupling, double driver_sd, double target_sd) {
    require(driver_sd > 0.0 && target_sd > 0.0, "mci_partial_correlation: sd");
    return coupling * driver_sd /
           std::sqrt(target_sd * target_sd + coupling * coupling * driver_sd * driver_sd);
}

double network_mann_kendall_sum(const Matrix& field) {
    require(!field.empty() && !field[0].empty(), "network_mann_kendall_sum: empty");
    double sum = 0.0;
    const std::size_t sites = field[0].size();
    for (std::size_t site = 0; site < sites; ++site) {
        std::vector<double> series(field.size(), 0.0);
        for (std::size_t t = 0; t < field.size(); ++t) series[t] = field[t][site];
        sum += mann_kendall_score(series);
    }
    return sum;
}

double rcamk_p_value(const Matrix& field, std::size_t resamples, double ar1,
                     unsigned seed) {
    require(!field.empty() && field.size() > 2 && resamples > 0 && std::fabs(ar1) < 1.0,
            "rcamk_p_value: inputs");
    const double observed = std::fabs(network_mann_kendall_sum(field));
    const Matrix correlation = correlation_matrix(field);
    const std::size_t sites = correlation.size();
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::vector<double> mean(sites, 0.0), innovation(sites, 0.0), residual(sites, 0.0);
    std::size_t exceedances = 0;
    Matrix sample(field.size(), std::vector<double>(sites, 0.0));
    for (std::size_t replicate = 0; replicate < resamples; ++replicate) {
        Matrix lower(sites, std::vector<double>(sites, 0.0));
        for (std::size_t i = 0; i < sites; ++i)
            for (std::size_t j = 0; j <= i; ++j) {
                double sum = correlation[i][j];
                for (std::size_t k = 0; k < j; ++k) sum -= lower[i][k] * lower[j][k];
                if (i == j)
                    lower[i][j] = std::sqrt(std::max(1e-12, sum));
                else
                    lower[i][j] = sum / lower[j][j];
            }
        for (std::size_t t = 0; t < field.size(); ++t) {
            for (double& value : innovation) value = normal(generator);
            for (std::size_t site = 0; site < sites; ++site) {
                double correlated = 0.0;
                for (std::size_t k = 0; k <= site; ++k)
                    correlated += lower[site][k] * innovation[k];
                residual[site] = t == 0 ? correlated : ar1 * residual[site] + correlated;
                sample[t][site] = mean[site] + residual[site];
            }
        }
        if (std::fabs(network_mann_kendall_sum(sample)) >= observed) ++exceedances;
    }
    return static_cast<double>(exceedances + 1) / static_cast<double>(resamples + 1);
}

double gaussian_vine_full_loglik(const Matrix& correlation, std::size_t samples) {
    require(!correlation.empty() && correlation.size() == correlation[0].size(),
            "gaussian_vine_full_loglik: square");
    const double det = determinant(correlation);
    require(det > 0.0, "gaussian_vine_full_loglik: singular correlation");
    return -0.5 * static_cast<double>(samples) * std::log(det);
}

double gaussian_vine_truncated_loglik(const Matrix& correlation, std::size_t samples) {
    double total = 0.0;
    for (std::size_t i = 0; i < correlation.size(); ++i)
        for (std::size_t j = i + 1; j < correlation.size(); ++j)
            total += -0.5 * static_cast<double>(samples) *
                     std::log(std::max(1e-15, 1.0 - correlation[i][j] * correlation[i][j]));
    return total;
}

double higher_tree_diagnostic(const Matrix& correlation, std::size_t samples) {
    return gaussian_vine_full_loglik(correlation, samples) -
           gaussian_vine_truncated_loglik(correlation, samples);
}

}  // namespace climate
