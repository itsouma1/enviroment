#include "ssa.hpp"

#include <cmath>
#include <stdexcept>

#include "regional.hpp"
#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double ols_slope(const std::vector<double>& series) {
    const std::size_t n = series.size();
    const double mean_time = 0.5 * (static_cast<double>(n) - 1.0);
    double numerator = 0.0, denominator = 0.0, mean_value = 0.0;
    for (double value : series) mean_value += value;
    mean_value /= static_cast<double>(n);
    for (std::size_t t = 0; t < n; ++t) {
        const double dt = static_cast<double>(t) - mean_time;
        numerator += dt * (series[t] - mean_value);
        denominator += dt * dt;
    }
    return numerator / denominator;
}

}  // namespace

Matrix ssa_trajectory(const std::vector<double>& series, std::size_t window) {
    require(window >= 2 && window < series.size(), "ssa_trajectory: window");
    const std::size_t rows = series.size() - window + 1;
    Matrix trajectory(rows, std::vector<double>(window, 0.0));
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < window; ++j) trajectory[i][j] = series[i + j];
    return trajectory;
}

std::vector<double> ssa_trend(const std::vector<double>& series, std::size_t window) {
    const Matrix trajectory = ssa_trajectory(series, window);
    const std::size_t rows = trajectory.size(), columns = window;
    Matrix covariance(columns, std::vector<double>(columns, 0.0));
    for (std::size_t a = 0; a < columns; ++a)
        for (std::size_t b = 0; b < columns; ++b) {
            double sum = 0.0;
            for (std::size_t i = 0; i < rows; ++i) sum += trajectory[i][a] * trajectory[i][b];
            covariance[a][b] = sum;
        }
    const EigenSystem eigen = symmetric_eigen(covariance);
    std::size_t leading = 0;
    for (std::size_t k = 1; k < eigen.values.size(); ++k)
        if (eigen.values[k] > eigen.values[leading]) leading = k;
    const std::vector<double>& loading = eigen.vectors[leading];
    std::vector<double> score(rows, 0.0);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < columns; ++j) score[i] += trajectory[i][j] * loading[j];
    const double norm = std::sqrt(eigen.values[leading]);
    require(norm > 0.0, "ssa_trend: degenerate series");
    std::vector<double> trend(series.size(), 0.0);
    std::vector<int> counts(series.size(), 0);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < columns; ++j) {
            trend[i + j] += score[i] * norm * loading[j];
            ++counts[i + j];
        }
    for (std::size_t t = 0; t < trend.size(); ++t) trend[t] /= static_cast<double>(counts[t]);
    return trend;
}

std::vector<double> sinusoid_fit(const std::vector<double>& series, double period,
                                 double& amplitude, double& phase) {
    require(period > 0.0, "sinusoid_fit: period");
    const double omega = 2.0 * 3.14159265358979323846 / period;
    double ss = 0.0, cc = 0.0, sc = 0.0, ys = 0.0, yc = 0.0;
    for (std::size_t t = 0; t < series.size(); ++t) {
        const double s = std::sin(omega * static_cast<double>(t));
        const double c = std::cos(omega * static_cast<double>(t));
        ss += s * s;
        cc += c * c;
        sc += s * c;
        ys += series[t] * s;
        yc += series[t] * c;
    }
    const double determinant = ss * cc - sc * sc;
    require(std::fabs(determinant) > 1e-18, "sinusoid_fit: degenerate period");
    const double a = (ys * cc - yc * sc) / determinant;
    const double b = (yc * ss - ys * sc) / determinant;
    amplitude = std::hypot(a, b);
    phase = std::atan2(b, a);
    std::vector<double> fitted(series.size(), 0.0);
    for (std::size_t t = 0; t < series.size(); ++t)
        fitted[t] = a * std::sin(omega * static_cast<double>(t)) +
                    b * std::cos(omega * static_cast<double>(t));
    return fitted;
}

SsaFeatures ssa_features(const std::vector<double>& series, std::size_t window,
                         double period) {
    require(series.size() > window + 1, "ssa_features: series");
    SsaFeatures features;
    const std::vector<double> trend = ssa_trend(series, window);
    features.trend = ols_slope(trend);
    std::vector<double> seasonal(series.size(), 0.0);
    for (std::size_t t = 0; t < series.size(); ++t) seasonal[t] = series[t] - trend[t];
    double amplitude = 0.0, phase = 0.0;
    const std::vector<double> fitted = sinusoid_fit(seasonal, period, amplitude, phase);
    features.amplitude = amplitude;
    features.phase = phase;
    std::vector<double> residual(series.size(), 0.0);
    for (std::size_t t = 0; t < series.size(); ++t) residual[t] = seasonal[t] - fitted[t];
    double mean_residual = 0.0;
    for (double value : residual) mean_residual += value;
    mean_residual /= static_cast<double>(residual.size());
    double variance = 0.0;
    for (double value : residual) variance += (value - mean_residual) * (value - mean_residual);
    features.residual = std::sqrt(variance / static_cast<double>(residual.size() - 1));
    return features;
}

}  // namespace climate
