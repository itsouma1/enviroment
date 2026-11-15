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

}  // namespace climate

#endif
