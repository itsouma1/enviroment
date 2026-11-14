#include "metrology.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "linalg.hpp"

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

Matrix transpose(const Matrix& matrix) {
    if (matrix.empty()) return {};
    Matrix result(matrix[0].size(), std::vector<double>(matrix.size(), 0.0));
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        require(matrix[i].size() == matrix[0].size(), "transpose: ragged matrix");
        for (std::size_t j = 0; j < matrix[i].size(); ++j) result[j][i] = matrix[i][j];
    }
    return result;
}

Matrix multiply(const Matrix& first, const Matrix& second) {
    require(!first.empty() && !second.empty(), "multiply: empty matrix");
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

Matrix add(const Matrix& first, const Matrix& second) {
    require(first.size() == second.size(), "add: size mismatch");
    Matrix result = first;
    for (std::size_t i = 0; i < first.size(); ++i) {
        require(first[i].size() == second[i].size(), "add: ragged matrix");
        for (std::size_t j = 0; j < first[i].size(); ++j) result[i][j] += second[i][j];
    }
    return result;
}

Matrix inverse(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    require(n > 0, "inverse: empty matrix");
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

}  // namespace

Matrix estimator_covariance(const Matrix& estimator, const Matrix& measurement_covariance) {
    return multiply(multiply(estimator, measurement_covariance), transpose(estimator));
}

Matrix structured_measurement_covariance(const Matrix& functional_map,
                                         const Matrix& functional_covariance,
                                         const Matrix& measurand_covariance,
                                         const Matrix& measurement_error_covariance) {
    return add(add(multiply(multiply(functional_map, functional_covariance),
                            transpose(functional_map)),
                   measurand_covariance),
               measurement_error_covariance);
}

Matrix gls_estimator(const Matrix& design, const Matrix& measurement_covariance) {
    const Matrix design_transpose = transpose(design);
    const Matrix inverse_covariance = inverse(measurement_covariance);
    const Matrix normal = multiply(multiply(design_transpose, inverse_covariance), design);
    return multiply(multiply(inverse(normal), design_transpose), inverse_covariance);
}

Matrix ols_estimator(const Matrix& design) {
    const Matrix design_transpose = transpose(design);
    return multiply(inverse(multiply(design_transpose, design)), design_transpose);
}

Matrix monte_carlo_covariance(const Matrix& realisations) {
    require(realisations.size() >= 2, "monte_carlo_covariance: need two realisations");
    const std::size_t count = realisations.size();
    const std::size_t dimensions = realisations[0].size();
    std::vector<double> mean(dimensions, 0.0);
    for (const auto& row : realisations) {
        require(row.size() == dimensions, "monte_carlo_covariance: ragged");
        for (std::size_t j = 0; j < dimensions; ++j) mean[j] += row[j];
    }
    for (double& value : mean) value /= static_cast<double>(count);
    Matrix covariance(dimensions, std::vector<double>(dimensions, 0.0));
    for (const auto& row : realisations)
        for (std::size_t i = 0; i < dimensions; ++i)
            for (std::size_t j = 0; j < dimensions; ++j)
                covariance[i][j] += (row[i] - mean[i]) * (row[j] - mean[j]);
    for (auto& row : covariance)
        for (double& value : row) value /= static_cast<double>(count - 1);
    return covariance;
}

std::vector<double> ndvi_partials(double near_infrared, double red) {
    const double denominator = near_infrared + red;
    require(std::fabs(denominator) > 0.0, "ndvi_partials: zero denominator");
    const double square = denominator * denominator;
    return {2.0 * red / square, -2.0 * near_infrared / square};
}

double ndvi_uncertainty(double near_infrared, double red, double uncertainty_nir,
                        double uncertainty_red, double covariance) {
    const std::vector<double> partials = ndvi_partials(near_infrared, red);
    const double variance = partials[0] * partials[0] * uncertainty_nir * uncertainty_nir +
                            partials[1] * partials[1] * uncertainty_red * uncertainty_red +
                            2.0 * partials[0] * partials[1] * covariance;
    require(variance >= 0.0, "ndvi_uncertainty: negative variance");
    return std::sqrt(variance);
}

std::vector<double> simple_kriging(const std::vector<double>& observations,
                                   const Matrix& observation_covariance,
                                   const Matrix& cross_covariance, double mean,
                                   const Matrix& observation_error_covariance) {
    const std::size_t n = observations.size();
    require(observation_covariance.size() == n, "simple_kriging: size mismatch");
    Matrix system = observation_covariance;
    if (!observation_error_covariance.empty()) system = add(system, observation_error_covariance);
    std::vector<double> weights;
    if (!solve_linear_system(system, observations, weights))
        throw std::runtime_error("simple_kriging: singular covariance");
    require(cross_covariance.size() == n, "simple_kriging: cross covariance size");
    const std::size_t grid = cross_covariance[0].size();
    std::vector<double> prediction(grid, mean);
    for (std::size_t g = 0; g < grid; ++g)
        for (std::size_t i = 0; i < n; ++i) prediction[g] += cross_covariance[i][g] * weights[i];
    return prediction;
}

Matrix kriging_error_covariance(const Matrix& observation_covariance,
                                const Matrix& cross_covariance,
                                const Matrix& observation_error_covariance) {
    Matrix system = observation_covariance;
    if (!observation_error_covariance.empty())
        system = add(system, observation_error_covariance);
    return multiply(multiply(transpose(cross_covariance), inverse(system)), cross_covariance);
}

double station_homogenisation_uncertainty(double single_station_uncertainty, int station_count) {
    require(station_count > 0, "station_homogenisation_uncertainty: station count");
    require(single_station_uncertainty >= 0.0, "station_homogenisation_uncertainty: uncertainty");
    return single_station_uncertainty / std::sqrt(static_cast<double>(station_count));
}

double sampling_uncertainty(double standard_deviation, int sample_count) {
    require(sample_count > 0, "sampling_uncertainty: sample count");
    require(standard_deviation >= 0.0, "sampling_uncertainty: standard deviation");
    return standard_deviation / std::sqrt(static_cast<double>(sample_count));
}

double bootstrap_mean_standard_error(const std::vector<double>& values, int resamples,
                                     unsigned seed) {
    require(!values.empty(), "bootstrap_mean_standard_error: empty values");
    require(resamples > 1, "bootstrap_mean_standard_error: resamples");
    std::mt19937 generator(seed);
    std::uniform_int_distribution<std::size_t> pick(0, values.size() - 1);
    std::vector<double> means;
    means.reserve(static_cast<std::size_t>(resamples));
    for (int r = 0; r < resamples; ++r) {
        double sum = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i) sum += values[pick(generator)];
        means.push_back(sum / static_cast<double>(values.size()));
    }
    double mean = 0.0;
    for (double value : means) mean += value;
    mean /= static_cast<double>(means.size());
    double variance = 0.0;
    for (double value : means) variance += (value - mean) * (value - mean);
    variance /= static_cast<double>(means.size() - 1);
    return std::sqrt(variance);
}

double simultaneous_nadir_bias(const std::vector<double>& target,
                               const std::vector<double>& reference) {
    require(!target.empty() && target.size() == reference.size(),
            "simultaneous_nadir_bias: size mismatch");
    double sum = 0.0;
    for (std::size_t i = 0; i < target.size(); ++i) sum += target[i] - reference[i];
    return sum / static_cast<double>(target.size());
}

double diurnal_harmonic(double hour, double mean, double cosine1, double sine1,
                        double cosine2, double sine2, double period) {
    require(period > 0.0, "diurnal_harmonic: period must be positive");
    const double angle = 2.0 * kPi * hour / period;
    return mean + cosine1 * std::cos(angle) + sine1 * std::sin(angle) +
           cosine2 * std::cos(2.0 * angle) + sine2 * std::sin(2.0 * angle);
}

}  // namespace climate
