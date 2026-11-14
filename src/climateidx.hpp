#ifndef CLIMATE_CLIMATEIDX_HPP
#define CLIMATE_CLIMATEIDX_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// ---------------------------------------------------------------------------
// Modal indices
// ---------------------------------------------------------------------------
double area_mean_anomaly(const std::vector<double>& anomalies,
                         const std::vector<double>& weights);

double relative_nino34(double nino34, double nino34_climatology, double tropical_mean,
                       double tropical_climatology, double scale = 1.0);

double amo_relative(double north_atlantic_anomaly, double global_anomaly);

double pdo_projection(const std::vector<double>& leading_eof,
                      const std::vector<double>& anomaly_field);

// Maximum overturning streamfunction over the supplied latitude/depth box.
double amoc_index(const std::vector<std::vector<double>>& streamfunction);

// ---------------------------------------------------------------------------
// Ocean heat content and Earth energy imbalance
// ---------------------------------------------------------------------------
double ocean_heat_content(double density, double specific_heat,
                          const std::vector<double>& temperatures,
                          const std::vector<double>& climatology,
                          const std::vector<double>& layer_thickness);

double layer_heat_content_anomaly(double density, double specific_heat,
                                  double temperature_anomaly, double layer_thickness);

double earth_energy_imbalance(double incoming, double outgoing_longwave,
                              double outgoing_shortwave);

double energy_imbalance_decomposition(double forcing, double feedback_parameter,
                                      double temperature_anomaly);

// ---------------------------------------------------------------------------
// Sea level
// ---------------------------------------------------------------------------
double sea_level_budget(double sterodynamic, double gravity_rotation_deformation,
                        double glacial_isostatic_adjustment, double residual);

double barystatic_sea_level(double greenland, double antarctic, double glaciers,
                            double terrestrial_water, double atmospheric_water);

double quadratic_sea_level(double offset, double rate, double acceleration, double time);

double quadratic_sea_level_acceleration(double acceleration);

// Broken linear model with a slope change at the change point.
double broken_linear_sea_level(double offset, double rate, double slope_change,
                               double time, double change_point);

// Local linear dynamic linear model state.
struct DynamicLinearState {
    double level = 0.0;
    double rate = 0.0;
    double covariance_00 = 1.0;
    double covariance_01 = 0.0;
    double covariance_11 = 1.0;

    // One Kalman update with observation variance V and state noise W.
    void update(double observation, double observation_variance, double level_noise,
                double rate_noise);
};

// ---------------------------------------------------------------------------
// Groundwater
// ---------------------------------------------------------------------------
double groundwater_storage_anomaly(double total_storage, double snow_water_equivalent,
                                   double surface_water, double soil_moisture);

// h(t) = h0 + sum_tau R(t - tau) exp(-tau / characteristic_time).
double recharge_convolution(double initial_head, const std::vector<double>& recharge,
                            double characteristic_time);

double storage_depletion(const std::vector<double>& drawdown,
                         const std::vector<double>& cell_area,
                         double storage_coefficient);

// ---------------------------------------------------------------------------
// Diagnostic arrays
// ---------------------------------------------------------------------------
// H(x, t) = (1/N_y) sum_y F(x, y, t); input is one field per time step.
Matrix hovmoller_average(const std::vector<Matrix>& fields);

// Central-difference time derivative, one-sided at the ends.
std::vector<double> phase_space_derivative(const std::vector<double>& values, double time_step);

// ---------------------------------------------------------------------------
// Stochastic volatility and regimes
// ---------------------------------------------------------------------------
double heston_variance_step(double variance, double mean_reversion, double long_run_variance,
                            double volatility_of_volatility, double time_step, double shock);

std::vector<double> markov_stationary_distribution(
    const std::vector<std::vector<double>>& transition);

double regime_expected_value(const std::vector<double>& regime_means,
                             const std::vector<double>& regime_probabilities);

}  // namespace climate

#endif
