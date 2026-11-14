#include "change.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double log_determinant(const Matrix& matrix) {
    const std::size_t n = matrix.size();
    require(n > 0, "log_determinant: empty matrix");
    Matrix lower(n, std::vector<double>(n, 0.0));
    double logdet = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        require(matrix[i].size() == n, "log_determinant: ragged matrix");
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k) sum -= lower[i][k] * lower[j][k];
            if (i == j) {
                require(sum > 0.0, "log_determinant: matrix not positive definite");
                lower[i][j] = std::sqrt(sum);
                logdet += 2.0 * std::log(lower[i][j]);
            } else {
                lower[i][j] = sum / lower[j][j];
            }
        }
    }
    return logdet;
}

}  // namespace

double piecewise_linear_value(double intercept, double rate, double rate_change, double time,
                              double change_point) {
    return intercept + rate * time + rate_change * std::max(time - change_point, 0.0);
}

double piecewise_rate_before(double rate) { return rate; }

double piecewise_rate_after(double rate, double rate_change) { return rate + rate_change; }

double effective_return_level(double probability, double location, double scale, double shape) {
    require(probability > 0.0 && probability < 1.0,
            "effective_return_level: probability outside (0,1)");
    require(scale > 0.0, "effective_return_level: scale must be positive");
    const double logarithm = -std::log(1.0 - probability);
    if (std::fabs(shape) < 1e-12) return location - scale * std::log(logarithm);
    return location + scale / shape * (std::pow(logarithm, -shape) - 1.0);
}

double return_period_derivative(double observation_interval, double cdf,
                                double cdf_derivative) {
    require(cdf < 1.0, "return_period_derivative: cdf must be below 1");
    const double denominator = (1.0 - cdf) * (1.0 - cdf);
    return observation_interval / denominator * cdf_derivative;
}

double ar2_step(double first_lag, double second_lag, double phi1, double phi2, double shock) {
    return phi1 * first_lag + phi2 * second_lag + shock;
}

double ar2_stability(double phi1, double phi2) { return phi1 + phi2; }

double oht_amplitude_term(double velocity_prime, double temperature_mean) {
    return velocity_prime * temperature_mean;
}

double oht_circulation_term(double velocity_mean, double temperature_prime) {
    return velocity_mean * temperature_prime;
}

double oht_covariance_term(double velocity_prime, double temperature_prime) {
    return velocity_prime * temperature_prime;
}

double fisher_information_gaussian_mean(double variance) {
    require(variance > 0.0, "fisher_information_gaussian_mean: variance must be positive");
    return 1.0 / variance;
}

double information_speed(double fisher_information, double parameter_rate) {
    require(fisher_information >= 0.0, "information_speed: negative information");
    return std::sqrt(fisher_information) * std::fabs(parameter_rate);
}

double clipped_standardised_value(double value, double median, double median_absolute_deviation,
                                  double threshold) {
    require(median_absolute_deviation > 0.0, "clipped_standardised_value: MAD");
    require(threshold > 0.0, "clipped_standardised_value: threshold must be positive");
    const double scaled = (value - median) / median_absolute_deviation;
    const double clipped = std::copysign(std::min(std::fabs(scaled), threshold), scaled);
    return median + median_absolute_deviation * clipped;
}

double local_contrast(double left_mean, double right_mean, double left_variance,
                      double right_variance, double regularisation, double half_width) {
    require(half_width > 0.0, "local_contrast: half width must be positive");
    const double denominator =
        std::sqrt(left_variance + right_variance + regularisation);
    require(denominator > 0.0, "local_contrast: zero denominator");
    return std::sqrt(half_width) * (right_mean - left_mean) / denominator;
}

double max_cusum_statistic(const std::vector<double>& contrasts) {
    require(!contrasts.empty(), "max_cusum_statistic: empty contrasts");
    double maximum = 0.0;
    for (double value : contrasts) maximum = std::max(maximum, std::fabs(value));
    return maximum;
}

double higher_criticism(const std::vector<double>& p_values) {
    require(!p_values.empty(), "higher_criticism: empty p-values");
    std::vector<double> sorted = p_values;
    std::sort(sorted.begin(), sorted.end());
    const double n = static_cast<double>(sorted.size());
    double maximum = 0.0;
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        const double p = std::min(std::max(sorted[i], 1e-12), 1.0 - 1e-12);
        const double statistic =
            std::sqrt(n) * (static_cast<double>(i + 1) / n - p) / std::sqrt(p * (1.0 - p));
        maximum = std::max(maximum, statistic);
    }
    return maximum;
}

double gaussian_covariance_log_likelihood_ratio(int first_count, int second_count,
                                               const Matrix& pooled_covariance,
                                               const Matrix& first_covariance,
                                               const Matrix& second_covariance) {
    require(first_count > 0 && second_count > 0,
            "gaussian_covariance_log_likelihood_ratio: counts");
    const double total = static_cast<double>(first_count + second_count);
    return 0.5 * (static_cast<double>(first_count) * log_determinant(first_covariance) +
                  static_cast<double>(second_count) * log_determinant(second_covariance) -
                  total * log_determinant(pooled_covariance));
}

double optimal_transport_cost_1d(const std::vector<double>& from, const std::vector<double>& to) {
    require(!from.empty() && from.size() == to.size(),
            "optimal_transport_cost_1d: size mismatch");
    std::vector<double> a = from, b = to;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    double cost = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double difference = a[i] - b[i];
        cost += difference * difference;
    }
    return cost / static_cast<double>(a.size());
}

std::vector<double> functional_slope(const Matrix& curves, const std::vector<double>& times) {
    require(!curves.empty() && curves.size() == times.size(),
            "functional_slope: size mismatch");
    const std::size_t positions = curves[0].size();
    require(positions > 0, "functional_slope: empty curves");
    for (const auto& row : curves)
        require(row.size() == positions, "functional_slope: ragged curves");

    const double n = static_cast<double>(times.size());
    double time_mean = 0.0;
    for (double time : times) time_mean += time;
    time_mean /= n;
    double time_scatter = 0.0;
    for (double time : times) time_scatter += (time - time_mean) * (time - time_mean);
    require(time_scatter > 0.0, "functional_slope: zero time scatter");

    std::vector<double> slope(positions, 0.0);
    for (std::size_t position = 0; position < positions; ++position) {
        double curve_mean = 0.0;
        for (const auto& row : curves) curve_mean += row[position];
        curve_mean /= n;
        double covariance = 0.0;
        for (std::size_t i = 0; i < times.size(); ++i)
            covariance += (times[i] - time_mean) * (curves[i][position] - curve_mean);
        slope[position] = covariance / time_scatter;
    }
    return slope;
}

double functional_slope_norm(const std::vector<double>& slope, double spacing) {
    require(spacing > 0.0, "functional_slope_norm: spacing must be positive");
    double total = 0.0;
    for (double value : slope) total += value * value * spacing;
    return total;
}

}  // namespace climate
