#include "tempac.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "linalg.hpp"
#include "regional.hpp"
#include "ssa.hpp"
#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double ols_slope(const std::vector<double>& values) {
    require(values.size() >= 2, "ols_slope: size");
    const std::size_t n = values.size();
    const double mean_time = 0.5 * (static_cast<double>(n) - 1.0);
    double mean_value = 0.0;
    for (double value : values) mean_value += value;
    mean_value /= static_cast<double>(n);
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t t = 0; t < n; ++t) {
        const double dt = static_cast<double>(t) - mean_time;
        numerator += dt * (values[t] - mean_value);
        denominator += dt * dt;
    }
    return numerator / denominator;
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
    require(!a.empty() && a[0].size() == x.size(), "multiply: vector size");
    std::vector<double> result(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < x.size(); ++j) result[i] += a[i][j] * x[j];
    return result;
}

}  // namespace

double first_harmonic_amplitude(const std::vector<double>& monthly_values) {
    require(monthly_values.size() == 12, "first_harmonic_amplitude: 12 months");
    double real = 0.0, imaginary = 0.0;
    for (std::size_t month = 0; month < 12; ++month) {
        const double angle = 2.0 * 3.14159265358979323846 *
                             static_cast<double>(month + 1) / 12.0;
        real += monthly_values[month] * std::cos(angle);
        imaginary -= monthly_values[month] * std::sin(angle);
    }
    return std::hypot(2.0 * real / 12.0, 2.0 * imaginary / 12.0);
}

std::vector<double> first_harmonic_amplitudes(const Matrix& monthly_by_site) {
    require(!monthly_by_site.empty(), "first_harmonic_amplitudes: empty");
    std::vector<double> amplitudes(monthly_by_site.size(), 0.0);
    for (std::size_t site = 0; site < monthly_by_site.size(); ++site)
        amplitudes[site] = first_harmonic_amplitude(monthly_by_site[site]);
    return amplitudes;
}

std::vector<double> first_eof(const Matrix& field) {
    require(!field.empty() && field[0].size() >= 2, "first_eof: field");
    const std::size_t time = field.size(), sites = field[0].size();
    Matrix covariance(sites, std::vector<double>(sites, 0.0));
    std::vector<double> mean(sites, 0.0);
    for (const auto& row : field)
        for (std::size_t s = 0; s < sites; ++s) mean[s] += row[s];
    for (double& value : mean) value /= static_cast<double>(time);
    for (std::size_t i = 0; i < sites; ++i)
        for (std::size_t j = 0; j < sites; ++j) {
            double sum = 0.0;
            for (std::size_t t = 0; t < time; ++t)
                sum += (field[t][i] - mean[i]) * (field[t][j] - mean[j]);
            covariance[i][j] = sum / static_cast<double>(time - 1);
        }
    const EigenSystem eigen = symmetric_eigen(covariance);
    std::size_t leading = 0;
    for (std::size_t k = 1; k < eigen.values.size(); ++k)
        if (eigen.values[k] > eigen.values[leading]) leading = k;
    return eigen.vectors[leading];
}

std::vector<double> project_fingerprint(const Matrix& field,
                                        const std::vector<double>& fingerprint) {
    require(!field.empty() && field[0].size() == fingerprint.size(),
            "project_fingerprint: size");
    std::vector<double> signal(field.size(), 0.0);
    for (std::size_t t = 0; t < field.size(); ++t)
        for (std::size_t s = 0; s < fingerprint.size(); ++s)
            signal[t] += field[t][s] * fingerprint[s];
    return signal;
}

std::vector<double> fingerprint_trends(const Matrix& field,
                                       const std::vector<double>& fingerprint,
                                       std::size_t window) {
    require(window >= 2 && window <= field.size(), "fingerprint_trends: window");
    const std::vector<double> signal = project_fingerprint(field, fingerprint);
    std::vector<double> trends(field.size() - window + 1, 0.0);
    for (std::size_t start = 0; start < trends.size(); ++start) {
        std::vector<double> slice(signal.begin() + static_cast<long>(start),
                                  signal.begin() + static_cast<long>(start + window));
        trends[start] = ols_slope(slice);
    }
    return trends;
}

double detection_p_value(double observed_trend, std::vector<double> noise_trends) {
    require(!noise_trends.empty(), "detection_p_value: empty noise");
    std::size_t exceedances = 0;
    for (double trend : noise_trends)
        if (trend >= observed_trend) ++exceedances;
    return static_cast<double>(exceedances) / static_cast<double>(noise_trends.size());
}

double pattern_signal_to_noise(double signal_trend, const std::vector<double>& noise_trends) {
    require(noise_trends.size() >= 2, "pattern_signal_to_noise: noise");
    double mean = 0.0;
    for (double trend : noise_trends) mean += trend;
    mean /= static_cast<double>(noise_trends.size());
    double variance = 0.0;
    for (double trend : noise_trends) variance += (trend - mean) * (trend - mean);
    variance /= static_cast<double>(noise_trends.size() - 1);
    require(variance > 0.0, "pattern_signal_to_noise: zero noise variance");
    return signal_trend / std::sqrt(variance);
}

double local_signal_to_noise(double ensemble_trend, double trend_spread) {
    require(trend_spread > 0.0, "local_signal_to_noise: spread");
    return ensemble_trend / trend_spread;
}

double signal_to_noise_length_ratio(double short_spread, double long_spread) {
    require(short_spread > 0.0 && long_spread > 0.0,
            "signal_to_noise_length_ratio: spreads");
    return short_spread / long_spread;
}

std::vector<double> cseof_harmonic_loading(const std::vector<double>& series, double period,
                                           int harmonic, double& amplitude, double& phase) {
    require(harmonic > 0 && period > 0.0, "cseof_harmonic_loading: inputs");
    return sinusoid_fit(series, period / static_cast<double>(harmonic), amplitude, phase);
}

std::vector<double> reduced_space_analysis(const std::vector<double>& forecast,
                                           const std::vector<double>& observations,
                                           const Matrix& observation_operator,
                                           const Matrix& background_covariance,
                                           const Matrix& observation_covariance) {
    const std::size_t state = forecast.size();
    const std::size_t measurements = observations.size();
    require(!observation_operator.empty() && observation_operator.size() == measurements &&
                observation_operator[0].size() == state &&
                background_covariance.size() == state &&
                background_covariance[0].size() == state &&
                observation_covariance.size() == measurements &&
                observation_covariance[0].size() == measurements,
            "reduced_space_analysis: size");
    const Matrix h_transpose = transpose(observation_operator);
    const Matrix bh = multiply(background_covariance, h_transpose);
    const Matrix hbh = multiply(observation_operator, bh);
    Matrix innovation_covariance = hbh;
    for (std::size_t i = 0; i < measurements; ++i)
        for (std::size_t j = 0; j < measurements; ++j)
            innovation_covariance[i][j] += observation_covariance[i][j];
    const std::vector<double> innovation =
        [&] {
            std::vector<double> values = multiply(observation_operator, forecast);
            for (std::size_t i = 0; i < measurements; ++i)
                values[i] = observations[i] - values[i];
            return values;
        }();
    Matrix gain(state, std::vector<double>(measurements, 0.0));
    for (std::size_t column = 0; column < measurements; ++column) {
        std::vector<double> rhs(measurements, 0.0);
        for (std::size_t row = 0; row < measurements; ++row) rhs[row] = bh[column][row];
        std::vector<double> solution;
        require(solve_linear_system(innovation_covariance, rhs, solution),
                "reduced_space_analysis: singular innovation covariance");
        for (std::size_t row = 0; row < state; ++row) gain[row][column] = bh[row][column];
        for (std::size_t row = 0; row < state; ++row) {
            gain[row][column] = 0.0;
            for (std::size_t k = 0; k < measurements; ++k)
                gain[row][column] += bh[row][k] * solution[k];
        }
    }
    std::vector<double> analysis = forecast;
    for (std::size_t row = 0; row < state; ++row)
        for (std::size_t column = 0; column < measurements; ++column)
            analysis[row] += gain[row][column] * innovation[column];
    return analysis;
}

}  // namespace climate
