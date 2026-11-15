#ifndef CLIMATE_OCEAN2_HPP
#define CLIMATE_OCEAN2_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "ml.hpp"

namespace climate {

double dynamic_class_mean(double class_effect, double covariate, double trend,
                          double process_effect);

double static_region_mixture(double positive_weight, double positive_trend,
                             double negative_weight, double negative_trend);

double class_switch_level_shift(double new_class_effect, double old_class_effect);

double depth_warming_gradient(double bottom_trend, double surface_trend, double depth_span);

// Events as inclusive [start, end] index pairs.
std::vector<std::pair<std::size_t, std::size_t>> marine_heatwave_events(
    const std::vector<double>& temperatures, const std::vector<double>& thresholds,
    std::size_t minimum_duration);

double marine_heatwave_intensity(const std::vector<double>& temperatures,
                                 const std::vector<double>& thresholds, std::size_t start,
                                 std::size_t end);

double marine_heatwave_frequency(std::size_t events, double years);

double mean_event_duration(
    const std::vector<std::pair<std::size_t, std::size_t>>& events);

double vertically_integrated_anomaly(const std::vector<double>& anomalies,
                                     const std::vector<double>& layer_thicknesses);

double marine_heatwave_cooccurrence(const std::vector<int>& bottom_flags,
                                    const std::vector<int>& surface_flags);

double independence_delta(const std::vector<int>& bottom_flags,
                          const std::vector<int>& surface_flags);

bool is_subsurface_intensified(double depth_of_maximum, double mixed_layer_depth);

std::size_t depth_of_maximum_anomaly(const std::vector<double>& profile);

struct Gaussian2D {
    double mean_x = 0.0;
    double mean_y = 0.0;
    double variance_x = 1.0;
    double variance_y = 1.0;
    double covariance = 0.0;
};

double gaussian_2d_density(const Gaussian2D& distribution, double x, double y);

double pcm_posterior(const Gaussian2D& first, const Gaussian2D& second, double first_weight,
                     double second_weight, double x, double y);

std::pair<double, double> pcm_posterior_gradient(const Gaussian2D& first,
                                                 const Gaussian2D& second,
                                                 double first_weight, double second_weight,
                                                 double x, double y);

double kmeans_objective(const Matrix& points, const std::vector<int>& assignments,
                        const Matrix& centroids);

double cluster_heat_content(const std::vector<double>& area_weights,
                            const std::vector<double>& responsibilities,
                            const std::vector<double>& layer_integrals, double rho_cp);

double cluster_weighted_mean(const std::vector<double>& values,
                             const std::vector<double>& weights);

double cluster_area(const std::vector<double>& area_weights, const std::vector<int>& labels,
                    int cluster);

std::vector<double> ridge_fingerprint(const Matrix& predictors,
                                      const std::vector<double>& target, double lambda);

double attribution_ambiguity(double first_component, double second_component);

double tracer_anomaly_decomposition(double proportion_change, double source_temperature,
                                    double proportion, double source_temperature_change,
                                    double residual);

}  // namespace climate

#endif
