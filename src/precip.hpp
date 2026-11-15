#ifndef CLIMATE_PRECIP_HPP
#define CLIMATE_PRECIP_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

struct PrecipitationField {
    Matrix occurrence;  // days x sites, 0/1
    Matrix amount;      // days x sites, mm on wet days, zero otherwise
};

double wet_day_probability(double p01, double p11);

// Two-state Markov occurrence with spatially correlated latent uniforms.
Matrix markov_occurrence_field(std::size_t rows, std::size_t columns, std::size_t days,
                               double p01, double p11, double neighbour_correlation,
                               unsigned seed);

Matrix markov_occurrence_field_sites(std::size_t rows, std::size_t columns,
                                     std::size_t days, const std::vector<double>& p01,
                                     const std::vector<double>& p11,
                                     double neighbour_correlation, unsigned seed);

PrecipitationField simulate_precipitation(std::size_t rows, std::size_t columns,
                                         std::size_t days, double p01, double p11,
                                         double neighbour_correlation, double shape,
                                         double rate, unsigned seed);

double wilks_correlation_distance_derivative(double distance, double d1, double d2);

double wilks_correlation_shape_derivative(double distance, double d1, double d2);

double temperature_correlation(double distance, double elevation_difference, double b0,
                               double b1, double b2);

// Locally weighted interpolation of one parameter.
double locally_weighted_parameter(const std::vector<double>& site_x,
                                  const std::vector<double>& site_y,
                                  const std::vector<double>& values, double grid_x,
                                  double grid_y, double bandwidth);

double locally_weighted_parameter_elevation(
    const std::vector<double>& site_x, const std::vector<double>& site_y,
    const std::vector<double>& site_elevation, const std::vector<double>& values,
    double grid_x, double grid_y, double grid_elevation, double bandwidth,
    double elevation_bandwidth);

}  // namespace climate

#endif
