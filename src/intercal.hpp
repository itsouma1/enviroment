#ifndef CLIMATE_INTERCAL_HPP
#define CLIMATE_INTERCAL_HPP

#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// GSICS linear intercalibration
// ---------------------------------------------------------------------------
double gsics_leo_radiance(double geo_radiance, double offset, double slope);

double gsics_uncertainty(double geo_radiance, double offset, double slope,
                         double offset_uncertainty, double slope_uncertainty,
                         double offset_slope_covariance);

// ---------------------------------------------------------------------------
// MSU/AMSU calibration and drift
// ---------------------------------------------------------------------------
double msu_linear_radiance(double cold_space_radiance, double slope, double earth_counts,
                           double cold_counts);

double msu_drift_term(double kappa, double mechanism_one, double mechanism_two);

double msu_tlt_synthesis(double tmt, double tut, double tls);

double msu_delta_radiance(double alpha, double beta, double time, double reference_time);

double msu_nonlinearity(double mu0, double k, double time, double reference_time);

double msu_z_factor(double scene_parameter, double earth_emissivity, double cold_emissivity,
                    double warm_emissivity);

double msu_calibrated_radiance(double reference_radiance, double delta_radiance,
                               double nonlinearity, double z_factor);

double diurnal_drift_correction(double observed, double factor, double drift_anomaly);

// ---------------------------------------------------------------------------
// Double difference
// ---------------------------------------------------------------------------
double double_difference(double target_observation, double reference_observation,
                         double target_model, double reference_model);

double predicted_brightness(double reference_observation, double target_model,
                            double reference_model);

// ---------------------------------------------------------------------------
// Diurnal and seasonal harmonics
// ---------------------------------------------------------------------------
double diurnal_anomaly(double alpha, double beta, double gamma, double local_solar_time);

double seasonal_amplitude(double base, double sine_coefficient, double cosine_coefficient,
                          int month);

// ---------------------------------------------------------------------------
// Microwave retrievals
// ---------------------------------------------------------------------------
double cwv_retrieval(const std::vector<double>& brightness_temperatures,
                     const std::vector<double>& coefficients);

double cwv_channel_sensitivity(double brightness_temperature, double coefficient);

double sst_first_step(const std::vector<double>& brightness_temperatures,
                      const std::vector<double>& linear_coefficients,
                      const std::vector<double>& quadratic_coefficients);

double sst_channel_sensitivity(double brightness_temperature, double linear_coefficient,
                               double quadratic_coefficient);

// ---------------------------------------------------------------------------
// SSA hybrid recalibration
// ---------------------------------------------------------------------------
double ssa_recalibration(double years_since_launch, double a0, double a1, double a2, double A1,
                         double A2, double A3, double A4);

double ssa_derivative(double years_since_launch, double a1, double a2, double A1, double A2,
                      double A3);

// ---------------------------------------------------------------------------
// CLARREO inverse-variance weighting
// ---------------------------------------------------------------------------
double clarreo_weight(double spatial_variance, double temporal_variance,
                      double noise_variance);

double clarreo_bias(const std::vector<double>& differences,
                    const std::vector<double>& weights);

double clarreo_uncertainty(const std::vector<double>& weights);

double clarreo_temporal_variance(double overlap_fraction, double sounder_variance,
                                 double footprint_count);

double clarreo_inflate(double variance, double sounder_samples);

// ---------------------------------------------------------------------------
// TRUTHS OBCS measurement function
// ---------------------------------------------------------------------------
double obcs_radiance(double calibration_coefficient, double nonlinearity_corrected_signal,
                     double offset, double background);

double angular_standard_deviation(const std::vector<double>& radiance_differences);

double bias_regression_slope(const std::vector<double>& target,
                             const std::vector<double>& reference);

// ---------------------------------------------------------------------------
// SNO variance budget
// ---------------------------------------------------------------------------
double sno_scene_variance(double instrument_variance, double scene_variance);

double clarreo_bias_variance(const std::vector<double>& weights);

double truths_total_uncertainty(double spectral, double geometric, double noise);

// One random-walk Kalman update; returns the updated mean and writes the
// updated variance, gain, innovation and innovation variance.
double kalman_random_walk(double mean, double variance, double process_variance, double dt,
                          double observation, double observation_variance, double& gain,
                          double& innovation, double& innovation_variance);

// CycleGAN losses.
double adversarial_loss(const std::vector<double>& real_scores,
                        const std::vector<double>& fake_scores);

double cycle_consistency_loss(const std::vector<double>& forward_cycle,
                              const std::vector<double>& original);

double cyclegan_objective(double forward_adversarial, double backward_adversarial,
                          double cycle_loss, double lambda);

double msu_z_counts(double slope, double earth_counts, double cold_counts, double warm_counts);

// Aerosol-cloud interaction and susceptibility.
double aci_susceptibility(double droplet_slope, double extinction_slope);

double nd_relative_bias(double cloud_fraction_slope, double cloud_fraction_intercept,
                        double cloud_fraction);

double lwp_error_variance(double effective_radius, double optical_depth,
                          double optical_depth_variance, double effective_radius_variance,
                          double cross_covariance);

double dml_treatment_effect(const std::vector<double>& residual_outcome,
                            const std::vector<double>& residual_treatment);

double g_computation_effect(const std::vector<double>& stage_effects,
                            const std::vector<double>& stage_probabilities);

double natural_effect_total(double direct_effect, double indirect_effect);

double mean_squared_error(const std::vector<double>& predicted,
                          const std::vector<double>& observed);

double coefficient_of_determination(const std::vector<double>& predicted,
                                    const std::vector<double>& observed);

double sno_variance_decomposition(double spatial_variance, double temporal_variance,
                                  double noise_variance);

double sno_noise_fraction(double noise_variance, double total_variance);

double sno_t_statistic(double mean_bias, double ensemble_standard_deviation, int event_count);

}  // namespace climate

#endif
