#include "magnitude.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {

double quadratic_trend_level(double intercept, double linear, double quadratic, double time) {
    return intercept + linear * time + quadratic * time * time;
}

double quadratic_trend_rate(double linear, double quadratic, double time) {
    return linear + 2.0 * quadratic * time;
}

double quadratic_trend_acceleration(double quadratic) { return 2.0 * quadratic; }

double linearithmic_trend_level(double intercept, double coefficient, double time) {
    if (!(time > 0.0)) throw std::runtime_error("linearithmic_trend_level: time must be positive");
    return intercept + coefficient * std::log(time) * time;
}

double linearithmic_trend_rate(double coefficient, double time) {
    if (!(time > 0.0)) throw std::runtime_error("linearithmic_trend_rate: time must be positive");
    return coefficient * (std::log(time) + 1.0);
}

double linearithmic_trend_acceleration(double coefficient, double time) {
    if (!(time > 0.0))
        throw std::runtime_error("linearithmic_trend_acceleration: time must be positive");
    return coefficient / time;
}

double step_magnitude(double level_after, double level_before) {
    return level_after - level_before;
}

double integrated_change(const std::vector<double>& values, double time_step) {
    if (values.size() < 2) throw std::runtime_error("integrated_change: need two values");
    if (!(time_step > 0.0)) throw std::runtime_error("integrated_change: positive time step");
    double total = 0.0;
    for (std::size_t i = 1; i < values.size(); ++i)
        total += 0.5 * (values[i] + values[i - 1]) * time_step;
    return total;
}

double return_level_shift(double future_return_level, double past_return_level) {
    return future_return_level - past_return_level;
}

double newey_west_variance(const std::vector<double>& residuals, int lag) {
    if (residuals.empty()) throw std::runtime_error("newey_west_variance: empty series");
    if (lag < 0) throw std::runtime_error("newey_west_variance: negative lag");
    const std::size_t n = residuals.size();
    if (lag >= static_cast<int>(n)) lag = static_cast<int>(n) - 1;
    double mean = 0.0;
    for (double residual : residuals) mean += residual;
    mean /= static_cast<double>(n);
    double variance = 0.0;
    for (double residual : residuals) variance += (residual - mean) * (residual - mean);
    variance /= static_cast<double>(n);
    for (int k = 1; k <= lag; ++k) {
        double autocovariance = 0.0;
        for (std::size_t t = static_cast<std::size_t>(k); t < n; ++t)
            autocovariance += (residuals[t] - mean) * (residuals[t - static_cast<std::size_t>(k)] - mean);
        autocovariance /= static_cast<double>(n);
        const double weight = 1.0 - static_cast<double>(k) / static_cast<double>(lag + 1);
        variance += 2.0 * weight * autocovariance;
    }
    return variance;
}

}  // namespace climate
