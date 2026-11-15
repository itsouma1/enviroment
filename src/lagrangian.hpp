#ifndef CLIMATE_LAGRANGIAN_HPP
#define CLIMATE_LAGRANGIAN_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "ml.hpp"

namespace climate {

// Hessian eigenvalues (k1 >= k2) of a 2x2 symmetric matrix.
std::pair<double, double> hessian_eigenvalues(double xx, double xy, double yy);

// SCAFET shape index, bounded to [-1, 1].
double shape_index(double k1, double k2);

double shape_index_k1_derivative(double k1, double k2);

double shape_index_k2_derivative(double k1, double k2);

// Grid-aware Gaussian smoothing with clamped boundaries.
Matrix gaussian_smooth(const Matrix& field, double sigma_cells);

// Eulerian OLS slope of a time series.
double linear_trend(const std::vector<double>& values);

// Lagrangian object statistics.
double object_frequency(const std::vector<int>& present);

double mean_object_intensity(const std::vector<double>& intensities,
                             const std::vector<int>& present);

// Distance to the nearest centroid at the next step, or -1 outside the radius.
double nearest_link_distance(double x, double y, double next_x, double next_y, double radius);

double object_population_trend(const std::vector<double>& counts);

double along_track_trend(const std::vector<double>& intensities);

}  // namespace climate

#endif
