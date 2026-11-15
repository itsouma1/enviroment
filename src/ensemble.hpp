#ifndef CLIMATE_ENSEMBLE_HPP
#define CLIMATE_ENSEMBLE_HPP

#include <cstddef>
#include <functional>
#include <vector>

#include "ml.hpp"

namespace climate {

Matrix spatial_correlation_matrix(std::size_t rows, std::size_t columns,
                                  double neighbour_correlation);

bool cholesky_factor(const Matrix& matrix, Matrix& factor);

// Wilks spatial correlation: exp(d1 * distance^d2).
double wilks_correlation(double distance, double d1, double d2);

Matrix simulate_occurrence_field(std::size_t rows, std::size_t columns, std::size_t days,
                                 double neighbour_correlation, double occurrence_probability,
                                 unsigned seed);

Matrix simulate_gamma_field(std::size_t rows, std::size_t columns, double shape,
                            double rate, unsigned seed);

void detrend_field(Matrix& field);

// Trend-free AR(1) field with C0^{d(i,j)} spatial correlation; rows are times.
Matrix simulate_zero_trend_field(std::size_t rows, std::size_t columns, std::size_t years,
                                 double neighbour_correlation, double temporal_correlation,
                                 unsigned seed);

std::size_t positive_trend_count(const Matrix& field);

std::vector<std::size_t> sign_test_null_ensemble(std::size_t rows, std::size_t columns,
                                                 std::size_t years,
                                                 double neighbour_correlation,
                                                 double temporal_correlation, int replicates,
                                                 unsigned seed);

void inject_uniform_trend(Matrix& field, double trend);

// Monte Carlo power of the Huth-Dubrovsky sign test.
double sign_test_power(std::size_t rows, std::size_t columns, std::size_t years,
                       double neighbour_correlation, double temporal_correlation,
                       double trend, double alpha, int null_replicates,
                       int signal_replicates, unsigned seed);

double central_difference(const std::function<double(double)>& function, double point,
                          double step);

}  // namespace climate

#endif
