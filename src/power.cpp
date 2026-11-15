#include "power.hpp"

#include <cmath>
#include <stdexcept>

#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double variance_inflation_factor(double autocorrelation) {
    require(std::fabs(autocorrelation) < 1.0, "variance_inflation_factor: autocorrelation");
    return (1.0 + autocorrelation) / (1.0 - autocorrelation);
}

double stationary_noise_variance(double innovation_sd, double autocorrelation) {
    require(innovation_sd > 0.0 && std::fabs(autocorrelation) < 1.0,
            "stationary_noise_variance: inputs");
    return innovation_sd * innovation_sd / (1.0 - autocorrelation * autocorrelation);
}

double minimum_detectable_trend_power(double noise_sd, double autocorrelation, double years,
                                      double alpha, double power) {
    require(noise_sd > 0.0 && years > 0.0 && alpha > 0.0 && alpha < 1.0 && power > 0.0 &&
                power < 1.0,
            "minimum_detectable_trend_power: inputs");
    const double z_alpha = normal_quantile(1.0 - 0.5 * alpha);
    const double z_power = normal_quantile(power);
    return (z_alpha + z_power) * noise_sd /
           std::pow(years, 1.5) *
           std::sqrt((1.0 + autocorrelation) / (1.0 - autocorrelation));
}

double detection_time_power(double noise_sd, double autocorrelation, double trend_magnitude,
                            double alpha, double power) {
    require(noise_sd > 0.0 && trend_magnitude > 0.0, "detection_time_power: inputs");
    const double z_alpha = normal_quantile(1.0 - 0.5 * alpha);
    const double z_power = normal_quantile(power);
    const double base = (z_alpha + z_power) * noise_sd / trend_magnitude *
                        std::sqrt((1.0 + autocorrelation) / (1.0 - autocorrelation));
    return std::pow(base, 2.0 / 3.0);
}

double trend_power(double trend_magnitude, double standard_error, double alpha) {
    require(standard_error > 0.0, "trend_power: standard error");
    const double ratio = std::fabs(trend_magnitude) / standard_error;
    const double critical = normal_quantile(1.0 - 0.5 * alpha);
    return normal_cdf(ratio - critical) + normal_cdf(-ratio - critical);
}

double signal_to_noise_ratio(double trend, double years, double noise_sd) {
    require(noise_sd > 0.0, "signal_to_noise_ratio: noise");
    return std::fabs(trend) * years / noise_sd;
}

double minimum_detectable_trend_from_snr(double snr50, double noise_sd, double years) {
    require(noise_sd > 0.0 && years > 0.0, "minimum_detectable_trend_from_snr: inputs");
    return snr50 * noise_sd / years;
}

double mdt_autocorrelation_derivative(double noise_sd, double autocorrelation, double years,
                                      double alpha, double power) {
    const double mdt = minimum_detectable_trend_power(noise_sd, autocorrelation, years, alpha,
                                                     power);
    return mdt / ((1.0 - autocorrelation) * (1.0 + autocorrelation));
}

double mdt_noise_derivative(double noise_sd, double autocorrelation, double years, double alpha,
                            double power) {
    return minimum_detectable_trend_power(noise_sd, autocorrelation, years, alpha, power) /
           noise_sd;
}

double detection_time_trend_derivative(double noise_sd, double autocorrelation,
                                       double trend_magnitude, double alpha, double power) {
    const double n = detection_time_power(noise_sd, autocorrelation, trend_magnitude, alpha,
                                          power);
    return -2.0 / 3.0 * n / trend_magnitude;
}

bool variable_is_feasible(double minimum_detectable, double target_trend) {
    return minimum_detectable <= std::fabs(target_trend);
}

double effective_field_significance(double nominal_alpha, double spatial_correlation,
                                    double grid_points) {
    require(nominal_alpha > 0.0 && nominal_alpha < 1.0 && grid_points > 0.0,
            "effective_field_significance: inputs");
    require(spatial_correlation >= 0.0 && spatial_correlation < 1.0,
            "effective_field_significance: correlation");
    const double effective = grid_points * (1.0 - spatial_correlation) /
                             (1.0 + spatial_correlation);
    if (effective <= 1.0) return nominal_alpha;
    return 1.0 - std::pow(1.0 - nominal_alpha, effective);
}

}  // namespace climate
