#ifndef CLIMATE_METROLOGY_HPP
#define CLIMATE_METROLOGY_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// ---------------------------------------------------------------------------
// Covariance-aware propagation
// ---------------------------------------------------------------------------
// Q_xhat = G Q_y G^T.
Matrix estimator_covariance(const Matrix& estimator, const Matrix& measurement_covariance);

// Q_y = A_bar Q_xE A_bar^T + Q_Y + Q_E.
Matrix structured_measurement_covariance(const Matrix& functional_map,
                                         const Matrix& functional_covariance,
                                         const Matrix& measurand_covariance,
                                         const Matrix& measurement_error_covariance);

// G_GLS = (A^T Q_y^-1 A)^-1 A^T Q_y^-1.
Matrix gls_estimator(const Matrix& design, const Matrix& measurement_covariance);

// G_OLS = (A^T A)^-1 A^T.
Matrix ols_estimator(const Matrix& design);

// Sample covariance from Monte Carlo realisations (one realisation per row).
Matrix monte_carlo_covariance(const Matrix& realisations);

// ---------------------------------------------------------------------------
// NDVI propagation
// ---------------------------------------------------------------------------
// ndvi(double, double) is declared in retrieval.hpp.
std::vector<double> ndvi_partials(double near_infrared, double red);

double ndvi_uncertainty(double near_infrared, double red, double uncertainty_nir,
                        double uncertainty_red, double covariance);

// ---------------------------------------------------------------------------
// Simple kriging
// ---------------------------------------------------------------------------
// Prediction at every grid point: C_cross^T (C_obs + E)^-1 y + mu.
std::vector<double> simple_kriging(const std::vector<double>& observations,
                                   const Matrix& observation_covariance,
                                   const Matrix& cross_covariance, double mean,
                                   const Matrix& observation_error_covariance = {});

// C_cross^T (C_obs + E)^-1 C_cross.
Matrix kriging_error_covariance(const Matrix& observation_covariance,
                                const Matrix& cross_covariance,
                                const Matrix& observation_error_covariance = {});

// ---------------------------------------------------------------------------
// Component error budgets
// ---------------------------------------------------------------------------
double station_homogenisation_uncertainty(double single_station_uncertainty, int station_count);

double sampling_uncertainty(double standard_deviation, int sample_count);

double bootstrap_mean_standard_error(const std::vector<double>& values, int resamples,
                                     unsigned seed = 0);

// ---------------------------------------------------------------------------
// Intercalibration helpers
// ---------------------------------------------------------------------------
double simultaneous_nadir_bias(const std::vector<double>& target,
                               const std::vector<double>& reference);

// Diurnal cycle harmonic reconstruction.
double diurnal_harmonic(double hour, double mean, double cosine1, double sine1,
                        double cosine2, double sine2, double period = 24.0);

}  // namespace climate

#endif
