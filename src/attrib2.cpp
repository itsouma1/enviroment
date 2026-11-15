#include "attrib2.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "linalg.hpp"
#include "regional.hpp"
#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Matrix invert(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    require(!matrix.empty() && matrix[0].size() == n, "invert: square");
    Matrix result(n, std::vector<double>(n, 0.0));
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<double> unit(n, 0.0), solution;
        unit[column] = 1.0;
        require(solve_linear_system(matrix, unit, solution), "invert: singular");
        for (std::size_t row = 0; row < n; ++row) result[row][column] = solution[row];
    }
    return result;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
    require(!a.empty() && !b.empty() && a[0].size() == b.size(), "multiply: size");
    Matrix result(a.size(), std::vector<double>(b[0].size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t k = 0; k < b.size(); ++k)
            for (std::size_t j = 0; j < b[0].size(); ++j)
                result[i][j] += a[i][k] * b[k][j];
    return result;
}

Matrix transpose(const Matrix& a) {
    Matrix result(a[0].size(), std::vector<double>(a.size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[0].size(); ++j) result[j][i] = a[i][j];
    return result;
}

std::vector<double> multiply(const Matrix& a, const std::vector<double>& x) {
    require(!a.empty() && a[0].size() == x.size(), "multiply: vector");
    std::vector<double> result(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < x.size(); ++j) result[i] += a[i][j] * x[j];
    return result;
}

double quadratic(const std::vector<double>& vector, const Matrix& matrix) {
    const std::vector<double> product = multiply(matrix, vector);
    double sum = 0.0;
    for (std::size_t i = 0; i < vector.size(); ++i) sum += vector[i] * product[i];
    return sum;
}

Matrix cholesky(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    Matrix lower(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k) sum -= lower[i][k] * lower[j][k];
            if (i == j) {
                require(sum > 0.0, "cholesky: not positive definite");
                lower[i][j] = std::sqrt(sum);
            } else {
                lower[i][j] = sum / lower[j][j];
            }
        }
    return lower;
}

Matrix lower_inverse(const Matrix& lower) {
    const std::size_t n = lower.size();
    Matrix result(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        result[i][i] = 1.0 / lower[i][i];
        for (std::size_t j = 0; j < i; ++j) {
            double sum = 0.0;
            for (std::size_t k = j; k < i; ++k) sum += lower[i][k] * result[k][j];
            result[i][j] = -sum / lower[i][i];
        }
    }
    return result;
}

std::vector<double> ols_detrend(const std::vector<double>& values, double& intercept,
                                double& slope) {
    const std::size_t n = values.size();
    const double mean_t = 0.5 * (static_cast<double>(n) + 1.0);
    double mean_y = 0.0;
    for (double value : values) mean_y += value;
    mean_y /= static_cast<double>(n);
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dt = static_cast<double>(i + 1) - mean_t;
        numerator += dt * (values[i] - mean_y);
        denominator += dt * dt;
    }
    slope = numerator / denominator;
    intercept = mean_y - slope * mean_t;
    std::vector<double> residual(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        residual[i] = values[i] - intercept - slope * static_cast<double>(i + 1);
    return residual;
}

}  // namespace

FingerprintFit gls_fingerprint(const std::vector<double>& observations,
                               const Matrix& signals, const Matrix& noise_precision) {
    require(!observations.empty() && signals.size() == observations.size() &&
                noise_precision.size() == observations.size(),
            "gls_fingerprint: size");
    const std::size_t n = observations.size(), coefficients = signals[0].size();
    const Matrix weighted_signals = multiply(noise_precision, signals);
    Matrix normal = multiply(transpose(signals), weighted_signals);
    const std::vector<double> rhs = multiply(transpose(signals), multiply(noise_precision, observations));
    FingerprintFit fit;
    fit.scaling.resize(coefficients);
    for (std::size_t j = 0; j < coefficients; ++j) {
        std::vector<double> unit(coefficients, 0.0), solution;
        unit[j] = 1.0;
        require(solve_linear_system(normal, unit, solution),
                "gls_fingerprint: singular normal matrix");
        fit.scaling[j] = 0.0;
        for (std::size_t k = 0; k < coefficients; ++k) fit.scaling[j] += solution[k] * rhs[k];
    }
    fit.covariance = invert(normal);
    fit.standard_error.resize(coefficients);
    for (std::size_t j = 0; j < coefficients; ++j)
        fit.standard_error[j] = std::sqrt(std::max(0.0, fit.covariance[j][j]));
    return fit;
}

double detection_z(double scaling, double standard_error) {
    require(standard_error > 0.0, "detection_z: standard error");
    return scaling / standard_error;
}

bool attribution_consistent(double scaling, double standard_error, double alpha) {
    require(standard_error > 0.0 && alpha > 0.0 && alpha < 1.0,
            "attribution_consistent: inputs");
    return std::fabs(scaling - 1.0) <= normal_quantile(1.0 - 0.5 * alpha) * standard_error;
}

int attribution_decision(const std::vector<double>& confidence_levels, double critical) {
    require(!confidence_levels.empty() && critical > 0.0 && critical < 1.0,
            "attribution_decision: inputs");
    int candidate = -1;
    for (std::size_t i = 0; i < confidence_levels.size(); ++i) {
        if (confidence_levels[i] < critical) continue;
        if (candidate >= 0) return -1;
        candidate = static_cast<int>(i);
    }
    return candidate;
}

double fingerprint_signal_to_noise(const std::vector<double>& fingerprint,
                                   const Matrix& noise_precision) {
    return std::sqrt(std::max(0.0, quadratic(fingerprint, noise_precision)));
}

double joint_signal_to_noise(const std::vector<double>& first,
                             const std::vector<double>& second, const Matrix& first_precision,
                             const Matrix& second_precision, const Matrix& cross_covariance) {
    require(first.size() == second.size() && first_precision.size() == first.size() &&
                second_precision.size() == first.size() &&
                cross_covariance.size() == first.size(),
            "joint_signal_to_noise: size");
    const double direct = quadratic(first, first_precision) +
                          quadratic(second, second_precision) +
                          2.0 * quadratic(first, cross_covariance);
    return std::sqrt(std::max(0.0, direct));
}

double joint_signal_to_noise_schur(const std::vector<double>& first,
                                   const std::vector<double>& second,
                                   const Matrix& first_covariance,
                                   const Matrix& second_covariance,
                                   const Matrix& cross_covariance) {
    require(!first.empty() && first.size() == second.size(), "joint_schur: size");
    const Matrix first_inverse = invert(first_covariance);
    const Matrix cross_transpose = transpose(cross_covariance);
    const Matrix schur = [&] {
        Matrix value = second_covariance;
        const Matrix correction = multiply(cross_transpose,
                                           multiply(first_inverse, cross_covariance));
        for (std::size_t i = 0; i < value.size(); ++i)
            for (std::size_t j = 0; j < value.size(); ++j) value[i][j] -= correction[i][j];
        return value;
    }();
    const std::vector<double> projected =
        [&] {
            std::vector<double> value = multiply(cross_transpose, multiply(first_inverse, first));
            for (std::size_t i = 0; i < value.size(); ++i) value[i] = second[i] - value[i];
            return value;
        }();
    const double first_term = quadratic(first, first_inverse);
    const double second_term = quadratic(projected, invert(schur));
    return std::sqrt(std::max(0.0, first_term + second_term));
}

EeFit estimating_equation_fit(const std::vector<double>& observations,
                              const Matrix& signals, const Matrix& noise_precision,
                              double scale_parameter, double spatial_dimension,
                              const std::vector<double>& model_runs) {
    require(!observations.empty() && signals.size() == observations.size() &&
                model_runs.size() == signals[0].size(),
            "estimating_equation_fit: size");
    const std::size_t n = observations.size(), coefficients = signals[0].size();
    const Matrix weighted = multiply(noise_precision, signals);
    Matrix precision = multiply(transpose(signals), weighted);
    for (std::size_t j = 0; j < coefficients; ++j)
        precision[j][j] += ee_bias_correction(scale_parameter, spatial_dimension, model_runs, j);
    const std::vector<double> rhs =
        multiply(transpose(signals), multiply(noise_precision, observations));
    EeFit fit;
    fit.scaling.resize(coefficients);
    for (std::size_t j = 0; j < coefficients; ++j) {
        std::vector<double> unit(coefficients, 0.0), solution;
        unit[j] = 1.0;
        require(solve_linear_system(precision, unit, solution),
                "estimating_equation_fit: singular precision");
        fit.scaling[j] = 0.0;
        for (std::size_t k = 0; k < coefficients; ++k) fit.scaling[j] += solution[k] * rhs[k];
    }
    fit.precision = precision;
    const EigenSystem eigen = symmetric_eigen(precision);
    const double largest = *std::max_element(eigen.values.begin(), eigen.values.end());
    const double smallest = *std::min_element(eigen.values.begin(), eigen.values.end());
    fit.condition_number = largest / std::max(smallest, 1e-300);
    return fit;
}

double ee_bias_correction(double scale_parameter, double spatial_dimension,
                          const std::vector<double>& model_runs, std::size_t coefficient) {
    require(coefficient < model_runs.size() && model_runs[coefficient] > 0.0,
            "ee_bias_correction: model runs");
    return scale_parameter * spatial_dimension / model_runs[coefficient];
}

double ee_pseudo_residual_scale(const std::vector<double>& scaling,
                                const std::vector<double>& model_runs,
                                double scale_parameter) {
    require(scaling.size() == model_runs.size() && scale_parameter > 0.0,
            "ee_pseudo_residual_scale: size");
    double sum = 0.0;
    for (std::size_t j = 0; j < scaling.size(); ++j) {
        require(model_runs[j] > 0.0, "ee_pseudo_residual_scale: model runs");
        sum += scaling[j] * scaling[j] / model_runs[j];
    }
    return std::sqrt(1.0 + scale_parameter * sum);
}

double confidence_interval_width(double standard_error, double alpha) {
    require(standard_error > 0.0 && alpha > 0.0 && alpha < 1.0,
            "confidence_interval_width: inputs");
    return 2.0 * normal_quantile(1.0 - 0.5 * alpha) * standard_error;
}

double box_cox(double value, double lambda) {
    require(value > 0.0, "box_cox: positive value");
    if (std::fabs(lambda) < 1e-12) return std::log(value);
    return (std::pow(value, lambda) - 1.0) / lambda;
}

double box_cox_lambda_derivative(double value, double lambda) {
    require(value > 0.0, "box_cox_lambda_derivative: positive value");
    if (std::fabs(lambda) < 1e-8)
        return 0.5 * std::log(value) * std::log(value);
    return (lambda * std::pow(value, lambda) * std::log(value) -
            (std::pow(value, lambda) - 1.0)) /
           (lambda * lambda);
}

double box_cox_marginal_effect(double beta, double value, double lambda) {
    require(value > 0.0, "box_cox_marginal_effect: positive value");
    return beta / std::pow(value, lambda - 1.0);
}

double box_cox_marginal_effect_derivative(double beta, double value, double lambda) {
    require(value > 0.0, "box_cox_marginal_effect_derivative: positive value");
    return -beta * (lambda - 1.0) / std::pow(value, lambda);
}

BierensMatrices bierens_matrices(const Matrix& series, double alpha) {
    require(!series.empty() && series[0].size() >= 1 && alpha > 0.0 && alpha < 1.0,
            "bierens_matrices: inputs");
    const std::size_t n = series.size(), variables = series[0].size();
    Matrix residual(n, std::vector<double>(variables, 0.0));
    for (std::size_t variable = 0; variable < variables; ++variable) {
        std::vector<double> column(n, 0.0);
        for (std::size_t t = 0; t < n; ++t) column[t] = series[t][variable];
        double intercept = 0.0, slope = 0.0;
        const std::vector<double> detrended = ols_detrend(column, intercept, slope);
        for (std::size_t t = 0; t < n; ++t) residual[t][variable] = detrended[t];
    }
    Matrix partial_sums(n, std::vector<double>(variables, 0.0));
    std::vector<double> running(variables, 0.0);
    for (std::size_t t = 0; t < n; ++t) {
        for (std::size_t variable = 0; variable < variables; ++variable)
            running[variable] += residual[t][variable];
        for (std::size_t variable = 0; variable < variables; ++variable)
            partial_sums[t][variable] = running[variable] / std::sqrt(static_cast<double>(n));
    }
    BierensMatrices matrices;
    matrices.partial_sums.assign(variables, std::vector<double>(variables, 0.0));
    for (std::size_t t = 0; t < n; ++t)
        for (std::size_t i = 0; i < variables; ++i)
            for (std::size_t j = 0; j < variables; ++j)
                matrices.partial_sums[i][j] +=
                    partial_sums[t][i] * partial_sums[t][j] / static_cast<double>(n);
    const std::size_t window = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::floor(std::pow(static_cast<double>(n), alpha))));
    matrices.moving_average.assign(variables, std::vector<double>(variables, 0.0));
    for (std::size_t t = window - 1; t < n; ++t) {
        std::vector<double> average(variables, 0.0);
        for (std::size_t j = 0; j < window; ++j)
            for (std::size_t variable = 0; variable < variables; ++variable)
                average[variable] += residual[t - j][variable];
        for (double& value : average) value /= static_cast<double>(window);
        for (std::size_t i = 0; i < variables; ++i)
            for (std::size_t j = 0; j < variables; ++j)
                matrices.moving_average[i][j] +=
                    average[i] * average[j] / static_cast<double>(window);
    }
    return matrices;
}

std::vector<double> bierens_eigenvalues(const Matrix& series, double alpha) {
    const BierensMatrices matrices = bierens_matrices(series, alpha);
    Matrix moving = matrices.moving_average;
    double scale = 1.0;
    for (std::size_t i = 0; i < moving.size(); ++i)
        scale = std::max(scale, std::fabs(moving[i][i]));
    for (int attempt = 0; attempt < 8; ++attempt) {
        try {
            const Matrix lower = cholesky(moving);
            const Matrix inverse_lower = lower_inverse(lower);
            const Matrix transformed = multiply(inverse_lower,
                                                multiply(matrices.partial_sums,
                                                         transpose(inverse_lower)));
            const EigenSystem eigen = symmetric_eigen(transformed);
            std::vector<double> values = eigen.values;
            std::sort(values.begin(), values.end());
            return values;
        } catch (const std::runtime_error&) {
            for (std::size_t i = 0; i < moving.size(); ++i)
                moving[i][i] += 1e-10 * scale;
        }
    }
    throw std::runtime_error("bierens_eigenvalues: singular moving-average matrix");
}

double bierens_statistic(const Matrix& series, double alpha, std::size_t r) {
    const std::vector<double> eigenvalues = bierens_eigenvalues(series, alpha);
    require(r < eigenvalues.size(), "bierens_statistic: r");
    return std::pow(static_cast<double>(series.size()), 1.0 - alpha) * eigenvalues[r];
}

double bierens_eigenvalue_sensitivity(const std::vector<double>& right_vector,
                                      const std::vector<double>& left_vector,
                                      const Matrix& perturbed_first,
                                      const Matrix& perturbed_second, double eigenvalue,
                                      const Matrix& second_matrix) {
    require(right_vector.size() == left_vector.size() && !second_matrix.empty(),
            "bierens_eigenvalue_sensitivity: size");
    const std::vector<double> first_product = multiply(perturbed_first, right_vector);
    const std::vector<double> second_product = multiply(perturbed_second, right_vector);
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t i = 0; i < right_vector.size(); ++i) {
        numerator += left_vector[i] * (first_product[i] - eigenvalue * second_product[i]);
        denominator += left_vector[i] * multiply(second_matrix, right_vector)[i];
    }
    require(std::fabs(denominator) > 1e-300,
            "bierens_eigenvalue_sensitivity: zero denominator");
    return numerator / denominator;
}

}  // namespace climate
