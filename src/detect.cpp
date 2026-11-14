#include "detect.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "linalg.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matrix) {
    const std::size_t n = matrix.size();
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<double> rhs(n, 0.0);
        rhs[column] = 1.0;
        std::vector<double> solution;
        if (!solve_linear_system(matrix, rhs, solution))
            throw std::runtime_error("detect: singular matrix");
        for (std::size_t row = 0; row < n; ++row) result[row][column] = solution[row];
    }
    return result;
}

// Cyclic Jacobi eigen-decomposition of a symmetric matrix.
void jacobi_eigen(const std::vector<std::vector<double>>& input,
                  std::vector<std::vector<double>>& eigenvectors,
                  std::vector<double>& eigenvalues) {
    const std::size_t n = input.size();
    std::vector<std::vector<double>> a = input;
    eigenvectors.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) eigenvectors[i][i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off_diagonal = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j) off_diagonal += a[i][j] * a[i][j];
        if (off_diagonal < 1e-24) break;
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (std::fabs(a[p][q]) < 1e-18) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double cosine = 1.0 / std::sqrt(t * t + 1.0);
                const double sine = t * cosine;
                for (std::size_t k = 0; k < n; ++k) {
                    const double a_pk = a[p][k];
                    const double a_qk = a[q][k];
                    a[p][k] = cosine * a_pk - sine * a_qk;
                    a[q][k] = sine * a_pk + cosine * a_qk;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double a_kp = a[k][p];
                    const double a_kq = a[k][q];
                    a[k][p] = cosine * a_kp - sine * a_kq;
                    a[k][q] = sine * a_kp + cosine * a_kq;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double v_kp = eigenvectors[k][p];
                    const double v_kq = eigenvectors[k][q];
                    eigenvectors[k][p] = cosine * v_kp - sine * v_kq;
                    eigenvectors[k][q] = sine * v_kp + cosine * v_kq;
                }
            }
        }
    }
    eigenvalues.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) eigenvalues[i] = a[i][i];
}

}  // namespace

double gev_cdf(double value, double location, double scale, double shape) {
    require(scale > 0.0, "gev_cdf: scale must be positive");
    const double z = (value - location) / scale;
    if (std::fabs(shape) < 1e-12) return std::exp(-std::exp(-z));
    const double inner = 1.0 + shape * z;
    if (inner <= 0.0) return value < location ? 0.0 : 1.0;
    return std::exp(-std::pow(inner, -1.0 / shape));
}

double gev_return_level(double probability, double location, double scale, double shape) {
    if (probability <= 0.0 || probability >= 1.0)
        throw std::runtime_error("gev_return_level: probability outside (0, 1)");
    const double logarithm = -std::log(1.0 - probability);
    if (std::fabs(shape) < 1e-12) return location - scale * std::log(logarithm);
    return location - scale / shape * (1.0 - std::pow(logarithm, -shape));
}

double gev_log_likelihood(double value, double location, double scale, double shape) {
    require(scale > 0.0, "gev_log_likelihood: scale must be positive");
    const double z = (value - location) / scale;
    if (std::fabs(shape) < 1e-12)
        return -std::log(scale) - z - std::exp(-z);
    const double inner = 1.0 + shape * z;
    require(inner > 0.0, "gev_log_likelihood: observation outside support");
    return -std::log(scale) - (1.0 + 1.0 / shape) * std::log(inner) -
           std::pow(inner, -1.0 / shape);
}

double gp_cdf(double exceedance, double scale, double shape) {
    require(scale > 0.0, "gp_cdf: scale must be positive");
    if (exceedance < 0.0) return 0.0;
    if (std::fabs(shape) < 1e-12) return 1.0 - std::exp(-exceedance / scale);
    const double inner = 1.0 + shape * exceedance / scale;
    if (inner <= 0.0) return 1.0;
    return 1.0 - std::pow(inner, -1.0 / shape);
}

double ou_stationary_variance(double noise_amplitude, double damping_rate) {
    require(damping_rate > 0.0, "ou_stationary_variance: damping must be positive");
    return noise_amplitude * noise_amplitude / (2.0 * damping_rate);
}

double ou_lag1_autocorrelation(double damping_rate, double time_step) {
    require(damping_rate >= 0.0, "ou_lag1_autocorrelation: negative damping");
    return std::exp(-damping_rate * time_step);
}

double ar1_estimate(const std::vector<double>& values) {
    require(values.size() >= 3, "ar1_estimate: need at least three values");
    const std::size_t n = values.size();
    const std::size_t span = n - 1;
    double mean_first = 0.0, mean_second = 0.0;
    for (std::size_t i = 0; i < span; ++i) {
        mean_first += values[i];
        mean_second += values[i + 1];
    }
    mean_first /= static_cast<double>(span);
    mean_second /= static_cast<double>(span);
    double covariance = 0.0, variance_first = 0.0, variance_second = 0.0;
    for (std::size_t i = 0; i < span; ++i) {
        const double a = values[i] - mean_first;
        const double b = values[i + 1] - mean_second;
        covariance += a * b;
        variance_first += a * a;
        variance_second += b * b;
    }
    const double denominator = std::sqrt(variance_first * variance_second);
    require(denominator > 0.0, "ar1_estimate: zero variance");
    return covariance / denominator;
}

double rate_mismatch_criterion(double variance_slope, double autocorrelation_slope,
                               double variance_scale, double autocorrelation_scale) {
    require(variance_scale > 0.0 && autocorrelation_scale > 0.0,
            "rate_mismatch_criterion: scales must be positive");
    return variance_slope / variance_scale + autocorrelation_slope / autocorrelation_scale;
}

double regime_shift_index(const std::vector<double>& values, double new_mean,
                          double regime_sigma) {
    require(!values.empty(), "regime_shift_index: empty window");
    require(regime_sigma > 0.0, "regime_shift_index: sigma must be positive");
    double sum = 0.0;
    for (double value : values) sum += value - new_mean;
    return sum / (static_cast<double>(values.size()) * regime_sigma);
}

double cusum_statistic(const std::vector<double>& values) {
    require(!values.empty(), "cusum_statistic: empty series");
    double mean = 0.0;
    for (double value : values) mean += value;
    mean /= static_cast<double>(values.size());
    double cumulative = 0.0, minimum = 0.0, maximum = 0.0;
    for (double value : values) {
        cumulative += value - mean;
        minimum = std::min(minimum, cumulative);
        maximum = std::max(maximum, cumulative - minimum);
    }
    return maximum;
}

std::vector<double> optimal_fingerprint(const std::vector<std::vector<double>>& responses,
                                        const std::vector<double>& observation,
                                        const std::vector<std::vector<double>>& noise_covariance) {
    require(!responses.empty() && responses.size() == observation.size(),
            "optimal_fingerprint: size mismatch");
    require(noise_covariance.size() == responses.size(),
            "optimal_fingerprint: covariance size mismatch");
    const std::size_t parameters = responses[0].size();
    const auto covariance_inverse = inverse(noise_covariance);

    for (const auto& row : responses)
        require(row.size() == parameters, "optimal_fingerprint: ragged responses");
    std::vector<std::vector<double>> normal(parameters, std::vector<double>(parameters, 0.0));
    std::vector<double> rhs(parameters, 0.0);
    for (std::size_t a = 0; a < parameters; ++a) {
        for (std::size_t b = 0; b < parameters; ++b) {
            double sum = 0.0;
            for (std::size_t i = 0; i < responses.size(); ++i)
                for (std::size_t j = 0; j < responses.size(); ++j)
                    sum += responses[i][a] * covariance_inverse[i][j] * responses[j][b];
            normal[a][b] = sum;
        }
        double sum = 0.0;
        for (std::size_t i = 0; i < responses.size(); ++i)
            for (std::size_t j = 0; j < responses.size(); ++j)
                sum += responses[i][a] * covariance_inverse[i][j] * observation[j];
        rhs[a] = sum;
    }
    std::vector<double> scaling;
    if (!solve_linear_system(normal, rhs, scaling))
        throw std::runtime_error("optimal_fingerprint: singular normal equations");
    return scaling;
}

double marchenko_pastur_bounds_lower(double variance, double aspect_ratio) {
    require(variance > 0.0 && aspect_ratio >= 0.0, "marchenko_pastur_bounds_lower: input");
    const double root = std::sqrt(aspect_ratio);
    return variance * (1.0 - root) * (1.0 - root);
}

double marchenko_pastur_bounds_upper(double variance, double aspect_ratio) {
    require(variance > 0.0 && aspect_ratio >= 0.0, "marchenko_pastur_bounds_upper: input");
    const double root = std::sqrt(aspect_ratio);
    return variance * (1.0 + root) * (1.0 + root);
}

double marchenko_pastur_density(double eigenvalue, double variance, double aspect_ratio) {
    require(variance > 0.0 && aspect_ratio > 0.0, "marchenko_pastur_density: input");
    const double lower = marchenko_pastur_bounds_lower(variance, aspect_ratio);
    const double upper = marchenko_pastur_bounds_upper(variance, aspect_ratio);
    if (eigenvalue <= lower || eigenvalue >= upper) return 0.0;
    return 1.0 / (2.0 * 3.14159265358979323846 * variance * eigenvalue) *
           std::sqrt((upper - eigenvalue) * (eigenvalue - lower));
}

int significant_eigenvalue_count(const std::vector<double>& eigenvalues, double variance,
                                 double aspect_ratio) {
    const double upper = marchenko_pastur_bounds_upper(variance, aspect_ratio);
    int count = 0;
    for (double eigenvalue : eigenvalues)
        if (eigenvalue > upper) ++count;
    return count;
}

Matrix rmt_reconstruction(const std::vector<double>& eigenvalues,
                          const std::vector<std::vector<double>>& eigenvectors) {
    require(eigenvalues.size() == eigenvectors.size(), "rmt_reconstruction: size mismatch");
    if (eigenvalues.empty()) return {};
    const std::size_t dimension = eigenvectors[0].size();
    Matrix result(dimension, std::vector<double>(dimension, 0.0));
    for (std::size_t mode = 0; mode < eigenvalues.size(); ++mode) {
        require(eigenvectors[mode].size() == dimension, "rmt_reconstruction: ragged vector");
        for (std::size_t i = 0; i < dimension; ++i)
            for (std::size_t j = 0; j < dimension; ++j)
                result[i][j] += eigenvalues[mode] * eigenvectors[mode][i] * eigenvectors[mode][j];
    }
    return result;
}

std::vector<double> total_least_squares(const Matrix& predictors,
                                        const std::vector<double>& response) {
    require(!predictors.empty() && predictors.size() == response.size(),
            "total_least_squares: size mismatch");
    const std::size_t parameters = predictors[0].size();
    const std::size_t columns = parameters + 1;
    Matrix augmented(predictors.size(), std::vector<double>(columns, 0.0));
    for (std::size_t row = 0; row < predictors.size(); ++row) {
        require(predictors[row].size() == parameters, "total_least_squares: ragged predictors");
        for (std::size_t column = 0; column < parameters; ++column)
            augmented[row][column] = predictors[row][column];
        augmented[row][parameters] = response[row];
    }
    Matrix gram(columns, std::vector<double>(columns, 0.0));
    for (std::size_t i = 0; i < columns; ++i)
        for (std::size_t j = 0; j < columns; ++j)
            for (const auto& row : augmented) gram[i][j] += row[i] * row[j];

    std::vector<std::vector<double>> eigenvectors;
    std::vector<double> eigenvalues;
    jacobi_eigen(gram, eigenvectors, eigenvalues);

    std::size_t smallest = 0;
    for (std::size_t i = 1; i < eigenvalues.size(); ++i)
        if (eigenvalues[i] < eigenvalues[smallest]) smallest = i;

    const double last = eigenvectors[columns - 1][smallest];
    require(std::fabs(last) > 1e-12, "total_least_squares: degenerate direction");
    std::vector<double> coefficients(parameters, 0.0);
    for (std::size_t i = 0; i < parameters; ++i)
        coefficients[i] = -eigenvectors[i][smallest] / last;
    return coefficients;
}

}  // namespace climate
