#ifndef CLIMATE_FLOW_HPP
#define CLIMATE_FLOW_HPP

#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "ml.hpp"

namespace climate {

using VelocityField =
    std::function<std::pair<double, double>(double, double, double)>;

// Finite-time Lyapunov exponent of the flow map over [-T, 0] or [0, T].
double ftle(double x0, double y0, double t0, double integration_time,
            const VelocityField& velocity, double jacobian_step, int steps);

Matrix ftle_field(double t0, double integration_time, double x_min, double x_max,
                  double y_min, double y_max, std::size_t nx, std::size_t ny,
                  const VelocityField& velocity, double jacobian_step, int steps);

// Arc length of the trajectory through both +tau and -tau.
double lagrangian_descriptor(double x0, double y0, double t0, double tau,
                             const VelocityField& velocity, int steps);

// Grid maximising the Lagrangian descriptor.
// Greedy connectivity-constrained ridge path over local maxima.
std::vector<std::pair<std::size_t, std::size_t>> constrained_ridge_path(
    const Matrix& field, std::size_t max_step, double max_difference);

// True if a full square ring around the candidate stays above field + threshold.
bool has_closed_contour(const Matrix& field, std::size_t row, std::size_t column,
                        double threshold, std::size_t max_radius);

std::pair<double, double> jet_axis(double t0, double tau, double x_min, double x_max,
                                   double y_min, double y_max, std::size_t nx,
                                   std::size_t ny, const VelocityField& velocity, int steps);

double occurrence_percent(std::size_t strong_days, std::size_t sampled_days);

double year_to_year_variability(const std::vector<double>& values);

double regime_shift_difference(double post_mean, double pre_mean);

// Cyclic moving average with a periodic (longitude) domain.
std::vector<double> cyclic_moving_average(const std::vector<double>& values, std::size_t window);

bool can_stitch(double distance, double time_gap, double radius, double max_gap);

// Greedy nearest-neighbour stitching of candidate positions frame by frame.
std::vector<std::vector<std::pair<std::size_t, std::size_t>>> stitch_tracks(
    const std::vector<std::vector<std::pair<double, double>>>& frames, double radius,
    double max_gap);

double jet_trend_correction(double dynamical_trend, double latitude_sensitivity,
                            double isentrope_trend);

// ---- Canny edge detection -------------------------------------------------

Matrix sobel_gradient(const Matrix& field, int axis);

double total_sobel_signal(double spatial_x, double spatial_y, double temporal,
                          double calibration);

double calibration_factor(double spatial_quantile, double temporal_quantile);

std::vector<double> non_maximum_suppression_1d(const std::vector<double>& values);

// 1 = strong edge, 2 = weak edge connected to strong, 0 = suppressed.
std::vector<int> hysteresis_threshold_1d(const std::vector<double>& values, double high,
                                         double low);

double cohens_d(double before_mean, double after_mean, double pooled_sd);

std::size_t connected_component_count(const std::vector<int>& mask, std::size_t nx,
                                      std::size_t ny, std::size_t nz);

}  // namespace climate

#endif
