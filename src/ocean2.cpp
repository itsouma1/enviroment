#include "ocean2.hpp"

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

struct Inverse2D {
    double a = 1.0, b = 0.0, d = 1.0;
};

Inverse2D inverse_2d(const Gaussian2D& distribution) {
    const double determinant = distribution.variance_x * distribution.variance_y -
                               distribution.covariance * distribution.covariance;
    require(determinant > 0.0, "inverse_2d: singular covariance");
    Inverse2D inverse;
    inverse.a = distribution.variance_y / determinant;
    inverse.b = -distribution.covariance / determinant;
    inverse.d = distribution.variance_x / determinant;
    return inverse;
}

}  // namespace

double dynamic_class_mean(double class_effect, double covariate, double trend,
                          double process_effect) {
    return class_effect + covariate * trend + process_effect;
}

double static_region_mixture(double positive_weight, double positive_trend,
                             double negative_weight, double negative_trend) {
    return positive_weight * positive_trend + negative_weight * negative_trend;
}

double class_switch_level_shift(double new_class_effect, double old_class_effect) {
    return new_class_effect - old_class_effect;
}

double depth_warming_gradient(double bottom_trend, double surface_trend, double depth_span) {
    require(depth_span > 0.0, "depth_warming_gradient: depth span");
    return (bottom_trend - surface_trend) / depth_span;
}

std::vector<std::pair<std::size_t, std::size_t>> marine_heatwave_events(
    const std::vector<double>& temperatures, const std::vector<double>& thresholds,
    std::size_t minimum_duration) {
    require(temperatures.size() == thresholds.size() && minimum_duration > 0,
            "marine_heatwave_events: inputs");
    std::vector<std::pair<std::size_t, std::size_t>> events;
    std::size_t start = 0;
    bool active = false;
    for (std::size_t t = 0; t < temperatures.size(); ++t) {
        const bool above = temperatures[t] > thresholds[t];
        if (above && !active) {
            active = true;
            start = t;
        } else if (!above && active) {
            active = false;
            if (t - start >= minimum_duration) events.push_back({start, t - 1});
        }
    }
    if (active && temperatures.size() - start >= minimum_duration)
        events.push_back({start, temperatures.size() - 1});
    return events;
}

double marine_heatwave_intensity(const std::vector<double>& temperatures,
                                 const std::vector<double>& thresholds, std::size_t start,
                                 std::size_t end) {
    require(temperatures.size() == thresholds.size() && start <= end &&
                end < temperatures.size(),
            "marine_heatwave_intensity: range");
    double sum = 0.0;
    for (std::size_t t = start; t <= end; ++t) sum += temperatures[t] - thresholds[t];
    return sum / static_cast<double>(end - start + 1);
}

double marine_heatwave_frequency(std::size_t events, double years) {
    require(years > 0.0, "marine_heatwave_frequency: years");
    return static_cast<double>(events) / years;
}

double mean_event_duration(
    const std::vector<std::pair<std::size_t, std::size_t>>& events) {
    require(!events.empty(), "mean_event_duration: no events");
    double total = 0.0;
    for (const auto& event : events)
        total += static_cast<double>(event.second - event.first + 1);
    return total / static_cast<double>(events.size());
}

double vertically_integrated_anomaly(const std::vector<double>& anomalies,
                                     const std::vector<double>& layer_thicknesses) {
    require(!anomalies.empty() && anomalies.size() == layer_thicknesses.size(),
            "vertically_integrated_anomaly: size");
    double total = 0.0;
    for (std::size_t i = 0; i < anomalies.size(); ++i)
        total += anomalies[i] * layer_thicknesses[i];
    return total;
}

double marine_heatwave_cooccurrence(const std::vector<int>& bottom_flags,
                                    const std::vector<int>& surface_flags) {
    require(!bottom_flags.empty() && bottom_flags.size() == surface_flags.size(),
            "marine_heatwave_cooccurrence: size");
    std::size_t both = 0;
    for (std::size_t i = 0; i < bottom_flags.size(); ++i)
        if (bottom_flags[i] != 0 && surface_flags[i] != 0) ++both;
    return static_cast<double>(both) / static_cast<double>(bottom_flags.size());
}

double independence_delta(const std::vector<int>& bottom_flags,
                          const std::vector<int>& surface_flags) {
    require(!bottom_flags.empty() && bottom_flags.size() == surface_flags.size(),
            "independence_delta: size");
    std::size_t bottom_days = 0, surface_days = 0, both = 0;
    for (std::size_t i = 0; i < bottom_flags.size(); ++i) {
        if (bottom_flags[i] != 0) ++bottom_days;
        if (surface_flags[i] != 0) ++surface_days;
        if (bottom_flags[i] != 0 && surface_flags[i] != 0) ++both;
    }
    require(surface_days > 0, "independence_delta: no surface events");
    const double conditional = static_cast<double>(both) / static_cast<double>(surface_days);
    const double marginal = static_cast<double>(bottom_days) /
                            static_cast<double>(bottom_flags.size());
    return conditional - marginal;
}

bool is_subsurface_intensified(double depth_of_maximum, double mixed_layer_depth) {
    return depth_of_maximum > mixed_layer_depth;
}

std::size_t depth_of_maximum_anomaly(const std::vector<double>& profile) {
    require(!profile.empty(), "depth_of_maximum_anomaly: empty profile");
    std::size_t best = 0;
    for (std::size_t i = 1; i < profile.size(); ++i)
        if (profile[i] > profile[best]) best = i;
    return best;
}

double gaussian_2d_density(const Gaussian2D& distribution, double x, double y) {
    const Inverse2D inverse = inverse_2d(distribution);
    const double dx = x - distribution.mean_x, dy = y - distribution.mean_y;
    const double quadratic =
        inverse.a * dx * dx + 2.0 * inverse.b * dx * dy + inverse.d * dy * dy;
    const double determinant = distribution.variance_x * distribution.variance_y -
                               distribution.covariance * distribution.covariance;
    return std::exp(-0.5 * quadratic) / (2.0 * kPi * std::sqrt(determinant));
}

double pcm_posterior(const Gaussian2D& first, const Gaussian2D& second, double first_weight,
                     double second_weight, double x, double y) {
    const double numerator = first_weight * gaussian_2d_density(first, x, y);
    const double denominator =
        numerator + second_weight * gaussian_2d_density(second, x, y);
    require(denominator > 0.0, "pcm_posterior: zero density");
    return numerator / denominator;
}

std::pair<double, double> pcm_posterior_gradient(const Gaussian2D& first,
                                                 const Gaussian2D& second,
                                                 double first_weight, double second_weight,
                                                 double x, double y) {
    const Inverse2D inverse_first = inverse_2d(first);
    const Inverse2D inverse_second = inverse_2d(second);
    const double posterior =
        pcm_posterior(first, second, first_weight, second_weight, x, y);
    const double dx_first = first.mean_x - x, dy_first = first.mean_y - y;
    const double dx_second = second.mean_x - x, dy_second = second.mean_y - y;
    const double ax = inverse_first.a * dx_first + inverse_first.b * dy_first;
    const double ay = inverse_first.b * dx_first + inverse_first.d * dy_first;
    const double bx = inverse_second.a * dx_second + inverse_second.b * dy_second;
    const double by = inverse_second.b * dx_second + inverse_second.d * dy_second;
    return {posterior * (ax - (posterior * ax + (1.0 - posterior) * bx)),
            posterior * (ay - (posterior * ay + (1.0 - posterior) * by))};
}

double kmeans_objective(const Matrix& points, const std::vector<int>& assignments,
                        const Matrix& centroids) {
    require(!points.empty() && points.size() == assignments.size() && !centroids.empty(),
            "kmeans_objective: size");
    double objective = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const int cluster = assignments[i];
        require(cluster >= 0 && static_cast<std::size_t>(cluster) < centroids.size(),
                "kmeans_objective: assignment outside centroids");
        require(points[i].size() == centroids[static_cast<std::size_t>(cluster)].size(),
                "kmeans_objective: dimension mismatch");
        for (std::size_t j = 0; j < points[i].size(); ++j) {
            const double difference =
                points[i][j] - centroids[static_cast<std::size_t>(cluster)][j];
            objective += difference * difference;
        }
    }
    return objective;
}

double cluster_heat_content(const std::vector<double>& area_weights,
                            const std::vector<double>& responsibilities,
                            const std::vector<double>& layer_integrals, double rho_cp) {
    require(!area_weights.empty() && area_weights.size() == responsibilities.size() &&
                area_weights.size() == layer_integrals.size(),
            "cluster_heat_content: size");
    double total = 0.0;
    for (std::size_t i = 0; i < area_weights.size(); ++i)
        total += area_weights[i] * responsibilities[i] * layer_integrals[i];
    return rho_cp * total;
}

double cluster_weighted_mean(const std::vector<double>& values,
                             const std::vector<double>& weights) {
    require(!values.empty() && values.size() == weights.size(),
            "cluster_weighted_mean: size");
    double weighted = 0.0, total = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        weighted += weights[i] * values[i];
        total += weights[i];
    }
    require(total > 0.0, "cluster_weighted_mean: zero weight");
    return weighted / total;
}

double cluster_area(const std::vector<double>& area_weights, const std::vector<int>& labels,
                    int cluster) {
    require(!area_weights.empty() && area_weights.size() == labels.size(),
            "cluster_area: size");
    double total = 0.0;
    for (std::size_t i = 0; i < labels.size(); ++i)
        if (labels[i] == cluster) total += area_weights[i];
    return total;
}

std::vector<double> ridge_fingerprint(const Matrix& predictors,
                                      const std::vector<double>& target, double lambda) {
    require(!predictors.empty() && predictors.size() == target.size() &&
                !predictors[0].empty(),
            "ridge_fingerprint: size");
    require(lambda >= 0.0, "ridge_fingerprint: lambda");
    const std::size_t parameters = predictors[0].size();
    Matrix normal(parameters, std::vector<double>(parameters, 0.0));
    std::vector<double> rhs(parameters, 0.0);
    for (std::size_t row = 0; row < predictors.size(); ++row) {
        require(predictors[row].size() == parameters, "ridge_fingerprint: ragged");
        for (std::size_t i = 0; i < parameters; ++i) {
            rhs[i] += predictors[row][i] * target[row];
            for (std::size_t j = 0; j < parameters; ++j)
                normal[i][j] += predictors[row][i] * predictors[row][j];
        }
    }
    for (std::size_t i = 0; i < parameters; ++i) normal[i][i] += lambda;
    std::vector<double> coefficients;
    require(solve_linear_system(normal, rhs, coefficients),
            "ridge_fingerprint: singular normal equations");
    return coefficients;
}

double attribution_ambiguity(double first_component, double second_component) {
    const double total = std::fabs(first_component) + std::fabs(second_component);
    if (total == 0.0) return 0.5;
    return std::fabs(first_component) / total;
}

double tracer_anomaly_decomposition(double proportion_change, double source_temperature,
                                    double proportion, double source_temperature_change,
                                    double residual) {
    return proportion_change * source_temperature +
           proportion * source_temperature_change + residual;
}

}  // namespace climate
