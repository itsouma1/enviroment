#ifndef CLIMATE_POWER_HPP
#define CLIMATE_POWER_HPP

namespace climate {

double variance_inflation_factor(double autocorrelation);

double stationary_noise_variance(double innovation_sd, double autocorrelation);

// MDT at a specified significance and power.
double minimum_detectable_trend_power(double noise_sd, double autocorrelation, double years,
                                      double alpha, double power);

// Detection time with an explicit power quantile.
double detection_time_power(double noise_sd, double autocorrelation, double trend_magnitude,
                            double alpha, double power);

double trend_power(double trend_magnitude, double standard_error, double alpha);

double signal_to_noise_ratio(double trend, double years, double noise_sd);

double minimum_detectable_trend_from_snr(double snr50, double noise_sd, double years);

double mdt_autocorrelation_derivative(double noise_sd, double autocorrelation, double years,
                                      double alpha, double power);

double mdt_noise_derivative(double noise_sd, double autocorrelation, double years, double alpha,
                            double power);

double detection_time_trend_derivative(double noise_sd, double autocorrelation,
                                       double trend_magnitude, double alpha, double power);

bool variable_is_feasible(double minimum_detectable, double target_trend);

// Effective field significance inflation under spatial correlation.
double effective_field_significance(double nominal_alpha, double spatial_correlation,
                                    double grid_points);

}  // namespace climate

#endif
