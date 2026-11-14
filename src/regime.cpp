#include "regime.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double symmetric_stable_sample(double alpha, double uniform, double exponential) {
    require(alpha > 0.0 && alpha <= 2.0, "symmetric_stable_sample: alpha outside (0, 2]");
    require(exponential > 0.0, "symmetric_stable_sample: exponential draw must be positive");
    const double phi = kPi * (uniform - 0.5);
    if (std::fabs(alpha - 1.0) < 1e-12) return std::tan(phi);
    const double first = std::sin(alpha * phi) / std::pow(std::cos(phi), 1.0 / alpha);
    const double second = std::pow(std::cos((1.0 - alpha) * phi) / exponential,
                                   (1.0 - alpha) / alpha);
    return first * second;
}

double climate_feedback_manifold(double lambda, double kappa, double temperature) {
    require(std::fabs(temperature) > 0.0, "climate_feedback_manifold: zero temperature");
    return lambda + kappa / temperature;
}

double temperature_tendency(double forcing, double lambda, double kappa, double temperature) {
    return forcing - lambda * temperature - kappa;
}

double logistic_transition(double covariate, double gamma, double threshold) {
    require(gamma > 0.0, "logistic_transition: gamma must be positive");
    return 1.0 / (1.0 + std::exp(-gamma * (covariate - threshold)));
}

double exponential_transition(double covariate, double gamma, double threshold) {
    require(gamma >= 0.0, "exponential_transition: gamma must be non-negative");
    return 1.0 - std::exp(-gamma * (covariate - threshold) * (covariate - threshold));
}

double multiple_lstr(double base, double regime1, double regime2, double transition1,
                     double transition2) {
    return base + regime1 * transition1 + regime2 * transition2;
}

double stochastic_volatility_observation(double log_volatility, double shock) {
    return std::exp(0.5 * log_volatility) * shock;
}

double stochastic_volatility_step(double log_volatility, double long_run_mean, double persistence,
                                  double volatility_of_volatility, double shock) {
    return long_run_mean + persistence * (log_volatility - long_run_mean) +
           volatility_of_volatility * shock;
}

double stochastic_volatility_mean_reversion(double autoregressive_coefficient, double time_step) {
    require(autoregressive_coefficient > 0.0 && time_step > 0.0,
            "stochastic_volatility_mean_reversion: inputs must be positive");
    return -std::log(autoregressive_coefficient) / time_step;
}

double stochastic_volatility_intercept(double intercept, double autoregressive_coefficient) {
    require(std::fabs(1.0 - autoregressive_coefficient) > 0.0,
            "stochastic_volatility_intercept: unit root");
    return intercept / (1.0 - autoregressive_coefficient);
}

double gjr_garch_variance(double constant, double shock_coefficient,
                          double asymmetry_coefficient, double persistence,
                          double previous_variance, double previous_shock) {
    const double asymmetry = previous_shock < 0.0 ? 1.0 : 0.0;
    return constant + shock_coefficient * previous_shock * previous_shock +
           asymmetry_coefficient * previous_shock * previous_shock * asymmetry +
           persistence * previous_variance;
}

double beta_weight(int lag, int maximum_lag, double omega1, double omega2) {
    require(maximum_lag > 0 && lag >= 1 && lag <= maximum_lag,
            "beta_weight: lag outside range");
    require(omega1 > 0.0 && omega2 > 0.0, "beta_weight: omega must be positive");
    const double x = static_cast<double>(lag) / static_cast<double>(maximum_lag + 1);
    return std::pow(x, omega1 - 1.0) * std::pow(1.0 - x, omega2 - 1.0);
}

double long_run_log_variance(double intercept, double slope,
                             const std::vector<double>& weights,
                             const std::vector<double>& covariate_lags) {
    require(weights.size() == covariate_lags.size(), "long_run_log_variance: size mismatch");
    double weighted = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i)
        weighted += weights[i] * covariate_lags[i];
    return intercept + slope * weighted;
}

double garch_midas_variance(double previous_component, double previous_return,
                            double previous_long_run, double alpha, double beta) {
    require(previous_long_run > 0.0, "garch_midas_variance: long-run component");
    return (1.0 - alpha - beta) +
           alpha * previous_return * previous_return / previous_long_run +
           beta * previous_component;
}

double asymmetric_laplace_density(double value, double location, double scale,
                                  double quantile) {
    require(scale > 0.0, "asymmetric_laplace_density: scale must be positive");
    require(quantile > 0.0 && quantile < 1.0, "asymmetric_laplace_density: quantile");
    const double z = (value - location) / scale;
    const double weight = z * (quantile - (z < 0.0 ? 1.0 : 0.0));
    return quantile * (1.0 - quantile) / scale * std::exp(-weight);
}

double average_seasonal_trends(const std::vector<double>& seasonal_trends) {
    require(!seasonal_trends.empty(), "average_seasonal_trends: empty");
    double total = 0.0;
    for (double trend : seasonal_trends) total += trend;
    return total / static_cast<double>(seasonal_trends.size());
}

double seasonal_trend_standard_error(const Matrix& covariance) {
    const std::size_t count = covariance.size();
    require(count > 0, "seasonal_trend_standard_error: empty covariance");
    for (const auto& row : covariance)
        require(row.size() == count, "seasonal_trend_standard_error: ragged covariance");
    double total = 0.0;
    for (const auto& row : covariance)
        for (double value : row) total += value;
    require(total >= 0.0, "seasonal_trend_standard_error: negative variance");
    return std::sqrt(total) / static_cast<double>(count);
}

double beast_seasonal(double time, double period, const std::vector<double>& cosines,
                      const std::vector<double>& sines) {
    require(period > 0.0, "beast_seasonal: period must be positive");
    require(cosines.size() == sines.size(), "beast_seasonal: size mismatch");
    double seasonal = 0.0;
    for (std::size_t j = 0; j < cosines.size(); ++j) {
        const double angle = 2.0 * kPi * static_cast<double>(j + 1) * time / period;
        seasonal += cosines[j] * std::cos(angle) + sines[j] * std::sin(angle);
    }
    return seasonal;
}

double beast_steps(double time, const std::vector<double>& change_points,
                   const std::vector<double>& magnitudes) {
    require(change_points.size() == magnitudes.size(), "beast_steps: size mismatch");
    double steps = 0.0;
    for (std::size_t k = 0; k < change_points.size(); ++k)
        if (time >= change_points[k]) steps += magnitudes[k];
    return steps;
}

}  // namespace climate
