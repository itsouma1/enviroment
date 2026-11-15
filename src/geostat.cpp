#include "geostat.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "linalg.hpp"

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Matrix identity(std::size_t n) {
    Matrix result(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) result[i][i] = 1.0;
    return result;
}

Matrix multiply(const Matrix& first, const Matrix& second) {
    const std::size_t rows = first.size();
    const std::size_t inner = second.size();
    const std::size_t columns = second[0].size();
    require(first[0].size() == inner, "multiply: dimension mismatch");
    Matrix result(rows, std::vector<double>(columns, 0.0));
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t k = 0; k < inner; ++k)
            for (std::size_t j = 0; j < columns; ++j)
                result[i][j] += first[i][k] * second[k][j];
    return result;
}

Matrix add_scaled(const Matrix& first, const Matrix& second, double scale) {
    require(first.size() == second.size(), "add_scaled: size mismatch");
    Matrix result = first;
    for (std::size_t i = 0; i < first.size(); ++i)
        for (std::size_t j = 0; j < first[i].size(); ++j)
            result[i][j] += scale * second[i][j];
    return result;
}

Matrix inverse(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    Matrix result(n, std::vector<double>(n, 0.0));
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<double> rhs(n, 0.0);
        rhs[column] = 1.0;
        std::vector<double> solution;
        if (!solve_linear_system(matrix, rhs, solution))
            throw std::runtime_error("inverse: singular matrix");
        for (std::size_t row = 0; row < n; ++row) result[row][column] = solution[row];
    }
    return result;
}

void jacobi_eigen(const Matrix& input, Matrix& eigenvectors, std::vector<double>& eigenvalues) {
    const std::size_t n = input.size();
    Matrix a = input;
    eigenvectors.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) eigenvectors[i][i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j) off += a[i][j] * a[i][j];
        if (off < 1e-24) break;
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (std::fabs(a[p][q]) < 1e-18) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (std::size_t k = 0; k < n; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double vkp = eigenvectors[k][p], vkq = eigenvectors[k][q];
                    eigenvectors[k][p] = c * vkp - s * vkq;
                    eigenvectors[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    eigenvalues.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) eigenvalues[i] = a[i][i];
}

}  // namespace

double spde_matern_variance(double kappa, double tau, int alpha, int dimension) {
    require(kappa > 0.0 && tau > 0.0 && dimension > 0, "spde_matern_variance: parameters");
    require(alpha > dimension / 2.0, "spde_matern_variance: alpha must exceed d/2");
    const double nu = static_cast<double>(alpha) - 0.5 * static_cast<double>(dimension);
    if (std::fabs(nu - 0.5) < 1e-12) return 1.0 / (4.0 * kPi * kappa * tau * tau);
    if (std::fabs(nu - 1.5) < 1e-12)
        return 1.0 / (8.0 * kPi * std::pow(kappa, 3.0) * tau * tau);
    if (std::fabs(nu - 2.5) < 1e-12)
        return 3.0 / (32.0 * kPi * std::pow(kappa, 5.0) * tau * tau);
    throw std::runtime_error("spde_matern_variance: only nu = 1/2, 3/2, 5/2 supported");
}

Matrix spde_precision(const Matrix& mass, const Matrix& stiffness, double kappa, double tau,
                      int alpha) {
    require(alpha >= 1, "spde_precision: alpha must be a positive integer");
    Matrix base = add_scaled(stiffness, mass, kappa * kappa);
    Matrix result = identity(base.size());
    for (int power = 0; power < alpha; ++power) result = multiply(result, base);
    const double scale = tau * tau;
    for (auto& row : result)
        for (double& value : row) value *= scale;
    return result;
}

std::vector<double> gmrf_conditional_mean(const Matrix& precision,
                                         const std::vector<double>& eta) {
    std::vector<double> mean;
    if (!solve_linear_system(precision, eta, mean))
        throw std::runtime_error("gmrf_conditional_mean: singular precision");
    return mean;
}

Matrix gmrf_conditional_covariance(const Matrix& precision) { return inverse(precision); }

double squared_exponential_kernel(double distance, double variance, double length_scale) {
    require(length_scale > 0.0, "squared_exponential_kernel: length scale");
    return variance * std::exp(-0.5 * distance * distance / (length_scale * length_scale));
}

std::vector<double> gp_posterior_mean(const std::vector<double>& train_x,
                                     const std::vector<double>& train_y,
                                     const std::vector<double>& predict_x,
                                     double signal_variance, double length_scale,
                                     double noise_variance) {
    require(!train_x.empty() && train_x.size() == train_y.size(),
            "gp_posterior_mean: size mismatch");
    const std::size_t n = train_x.size();
    Matrix covariance(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j)
            covariance[i][j] = squared_exponential_kernel(std::fabs(train_x[i] - train_x[j]),
                                                          signal_variance, length_scale);
        covariance[i][i] += noise_variance;
    }
    std::vector<double> alpha;
    if (!solve_linear_system(covariance, train_y, alpha))
        throw std::runtime_error("gp_posterior_mean: singular covariance");
    std::vector<double> mean(predict_x.size(), 0.0);
    for (std::size_t p = 0; p < predict_x.size(); ++p)
        for (std::size_t i = 0; i < n; ++i)
            mean[p] += squared_exponential_kernel(std::fabs(predict_x[p] - train_x[i]),
                                                  signal_variance, length_scale) * alpha[i];
    return mean;
}

std::vector<double> gp_posterior_variance(const std::vector<double>& train_x,
                                         const std::vector<double>& predict_x,
                                         double signal_variance, double length_scale,
                                         double noise_variance) {
    const std::size_t n = train_x.size();
    Matrix covariance(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j)
            covariance[i][j] = squared_exponential_kernel(std::fabs(train_x[i] - train_x[j]),
                                                          signal_variance, length_scale);
        covariance[i][i] += noise_variance;
    }
    const Matrix inverse_covariance = inverse(covariance);
    std::vector<double> variance(predict_x.size(), signal_variance);
    for (std::size_t p = 0; p < predict_x.size(); ++p) {
        std::vector<double> kernel(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            kernel[i] = squared_exponential_kernel(std::fabs(predict_x[p] - train_x[i]),
                                                   signal_variance, length_scale);
        double quadratic = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                quadratic += kernel[i] * inverse_covariance[i][j] * kernel[j];
        variance[p] = std::max(signal_variance - quadratic, 0.0);
    }
    return variance;
}

double mad_envelope_statistic(const std::vector<double>& observed,
                              const std::vector<double>& ensemble_member,
                              const std::vector<double>& pointwise_scale) {
    require(observed.size() == ensemble_member.size() &&
                observed.size() == pointwise_scale.size(),
            "mad_envelope_statistic: size mismatch");
    double maximum = 0.0;
    for (std::size_t k = 0; k < observed.size(); ++k) {
        require(pointwise_scale[k] > 0.0, "mad_envelope_statistic: zero scale");
        maximum = std::max(maximum, std::fabs(ensemble_member[k] - observed[k]) /
                                        pointwise_scale[k]);
    }
    return maximum;
}

double mad_envelope_p_value(const std::vector<double>& observed,
                            const std::vector<std::vector<double>>& ensemble,
                            const std::vector<double>& pointwise_scale) {
    require(!ensemble.empty(), "mad_envelope_p_value: empty ensemble");
    const double observed_statistic =
        mad_envelope_statistic(observed, observed, pointwise_scale);
    int exceedances = 0;
    for (const auto& member : ensemble)
        if (mad_envelope_statistic(observed, member, pointwise_scale) >= observed_statistic)
            ++exceedances;
    return static_cast<double>(exceedances) / static_cast<double>(ensemble.size());
}

double gtwr_distance(double u1, double v1, double t1, double u2, double v2, double t2,
                     double spatial_scale, double temporal_scale) {
    require(spatial_scale >= 0.0 && temporal_scale >= 0.0, "gtwr_distance: negative scale");
    return std::sqrt(spatial_scale * ((u1 - u2) * (u1 - u2) + (v1 - v2) * (v1 - v2)) +
                     temporal_scale * (t1 - t2) * (t1 - t2));
}

std::vector<double> gtwr_coefficients(const Matrix& predictors, const std::vector<double>& response,
                                      const std::vector<double>& u, const std::vector<double>& v,
                                      const std::vector<double>& t, std::size_t target,
                                      double spatial_scale, double temporal_scale,
                                      double bandwidth) {
    require(!predictors.empty() && predictors.size() == response.size() &&
                predictors.size() == u.size() && predictors.size() == v.size() &&
                predictors.size() == t.size(),
            "gtwr_coefficients: size mismatch");
    require(target < predictors.size(), "gtwr_coefficients: target out of range");
    require(bandwidth > 0.0, "gtwr_coefficients: bandwidth");
    const std::size_t parameters = predictors[0].size();
    Matrix normal(parameters, std::vector<double>(parameters, 0.0));
    std::vector<double> rhs(parameters, 0.0);
    for (std::size_t i = 0; i < predictors.size(); ++i) {
        require(predictors[i].size() == parameters, "gtwr_coefficients: ragged");
        const double distance =
            gtwr_distance(u[target], v[target], t[target], u[i], v[i], t[i], spatial_scale,
                          temporal_scale);
        const double weight = std::exp(-distance * distance / (bandwidth * bandwidth));
        for (std::size_t a = 0; a < parameters; ++a) {
            rhs[a] += weight * predictors[i][a] * response[i];
            for (std::size_t b = 0; b < parameters; ++b)
                normal[a][b] += weight * predictors[i][a] * predictors[i][b];
        }
    }
    std::vector<double> coefficients;
    if (!solve_linear_system(normal, rhs, coefficients))
        throw std::runtime_error("gtwr_coefficients: singular local design");
    return coefficients;
}

double gtwr_aicc(int sample_count, double residual_sigma, double hat_trace) {
    require(sample_count > 0 && residual_sigma > 0.0, "gtwr_aicc: inputs");
    const double n = static_cast<double>(sample_count);
    require(n - 2.0 - hat_trace > 0.0, "gtwr_aicc: hat trace too large");
    return 2.0 * n * std::log(residual_sigma) + n * std::log(2.0 * kPi) +
           n * (n + hat_trace) / (n - 2.0 - hat_trace);
}

double spherical_variogram(double distance, double nugget, double sill, double range) {
    require(range > 0.0, "spherical_variogram: range");
    if (distance >= range) return nugget + sill;
    const double ratio = distance / range;
    return nugget + sill * (1.5 * ratio - 0.5 * ratio * ratio * ratio);
}

double exponential_variogram(double distance, double nugget, double sill, double range) {
    require(range > 0.0, "exponential_variogram: range");
    return nugget + sill * (1.0 - std::exp(-distance / range));
}

double gaussian_variogram(double distance, double nugget, double sill, double range) {
    require(range > 0.0, "gaussian_variogram: range");
    return nugget + sill * (1.0 - std::exp(-distance * distance / (range * range)));
}

double matern_variogram(double distance, double nugget, double sill, double range,
                        double smoothness) {
    require(range > 0.0, "matern_variogram: range");
    if (distance == 0.0) return nugget;
    if (std::fabs(smoothness - 0.5) < 1e-12)
        return nugget + sill * (1.0 - std::exp(-distance / range));
    if (std::fabs(smoothness - 1.5) < 1e-12) {
        const double scaled = std::sqrt(3.0) * distance / range;
        return nugget + sill * (1.0 - (1.0 + scaled) * std::exp(-scaled));
    }
    if (std::fabs(smoothness - 2.5) < 1e-12) {
        const double scaled = std::sqrt(5.0) * distance / range;
        return nugget + sill * (1.0 - (1.0 + scaled + scaled * scaled / 3.0) *
                                          std::exp(-scaled));
    }
    throw std::runtime_error("matern_variogram: only nu = 1/2, 3/2, 5/2 supported");
}

Matrix weighted_covariance(const Matrix& residuals, const std::vector<double>& weights) {
    require(!residuals.empty() && residuals.size() == weights.size(),
            "weighted_covariance: size mismatch");
    const std::size_t dimensions = residuals[0].size();
    for (const auto& row : residuals)
        require(row.size() == dimensions, "weighted_covariance: ragged");
    double weight_sum = 0.0;
    for (double weight : weights) {
        require(weight >= 0.0, "weighted_covariance: negative weight");
        weight_sum += weight;
    }
    require(weight_sum > 0.0, "weighted_covariance: zero weight");
    Matrix covariance(dimensions, std::vector<double>(dimensions, 0.0));
    for (std::size_t i = 0; i < residuals.size(); ++i)
        for (std::size_t a = 0; a < dimensions; ++a)
            for (std::size_t b = 0; b < dimensions; ++b)
                covariance[a][b] += weights[i] * residuals[i][a] * residuals[i][b];
    for (auto& row : covariance)
        for (double& value : row) value /= weight_sum;
    return covariance;
}

Matrix weighted_pca_residual(const Matrix& residuals, const std::vector<double>& weights,
                             int components) {
    const Matrix covariance = weighted_covariance(residuals, weights);
    const std::size_t dimensions = covariance.size();
    require(components >= 0 && components <= static_cast<int>(dimensions),
            "weighted_pca_residual: component count");
    Matrix eigenvectors;
    std::vector<double> eigenvalues;
    jacobi_eigen(covariance, eigenvectors, eigenvalues);

    std::vector<std::size_t> order(dimensions);
    for (std::size_t i = 0; i < dimensions; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return eigenvalues[a] > eigenvalues[b]; });

    Matrix projected = residuals;
    for (std::size_t i = 0; i < projected.size(); ++i) {
        for (int k = 0; k < components; ++k) {
            const std::size_t mode = order[static_cast<std::size_t>(k)];
            double coefficient = 0.0;
            for (std::size_t a = 0; a < dimensions; ++a)
                coefficient += eigenvectors[a][mode] * residuals[i][a];
            for (std::size_t a = 0; a < dimensions; ++a)
                projected[i][a] -= coefficient * eigenvectors[a][mode];
        }
    }
    return projected;
}

}  // namespace climate
