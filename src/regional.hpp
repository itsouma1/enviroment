#ifndef CLIMATE_REGIONAL_HPP
#define CLIMATE_REGIONAL_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

// ---- Huth-Dubrovsky sign test -------------------------------------------

double sign_difference_statistic(std::size_t positive, std::size_t negative);

// Huth-Dubrovsky percentage form: S = (n_+/N - 1) * 100.
double sign_test_percentage(std::size_t positive, std::size_t total);

// Exact one-sided binomial tail P(M >= observed) with p = 0.5.
double sign_test_binomial_p_value(std::size_t observed, std::size_t trials);

double sign_test_empirical_p_value(std::size_t observed,
                                   const std::vector<std::size_t>& ensemble_counts);

double effective_sample_size(double samples, double spatial_autocorrelation);

// S_crit ~ z_{1-alpha} sqrt((1+C0)/(N(1-C0))).
double sign_test_critical_value(double samples, double spatial_autocorrelation,
                                double alpha);

double sign_test_critical_value_derivative(double samples, double spatial_autocorrelation,
                                           double alpha);

// ---- Cai zonal decomposition --------------------------------------------

double zonal_mean_at_latitude(const std::vector<double>& trends,
                              const std::vector<int>& latitude_index, int latitude);

std::vector<double> zonally_asymmetric_trend(const std::vector<double>& trends,
                                             const std::vector<int>& latitude_index);

double regional_mean(const std::vector<double>& field, const std::vector<int>& region_mask);

double sign_agreement(const std::vector<double>& observed,
                      const std::vector<double>& estimated,
                      const std::vector<int>& land_mask);

// ---- Hierarchical GEV partial pooling -----------------------------------

double gev_location_score(double value, double location, double scale, double shape);

double gev_trend_score(double value, double location, double scale, double shape,
                       double time);

double shrinkage_log_prior(double site_trend, double cluster_mean,
                           double between_site_variance);

double shrinkage_derivative(double site_trend, double cluster_mean,
                            double between_site_variance);

double pooled_posterior_variance(double at_site_variance, double between_site_variance);

// SE_pooled / SE_no-pool = (1 + tau^2 * data precision)^(-1/2).
double pooling_standard_error_ratio(double between_site_variance, double data_precision);

double pooled_cluster_mean(const std::vector<double>& site_trends,
                           const std::vector<double>& precisions);

// ---- REOF regionalisation ----------------------------------------------

struct EigenSystem {
    std::vector<double> values;
    Matrix vectors;
};

EigenSystem symmetric_eigen(const Matrix& matrix);

// North et al. eigenvalue standard error: lambda sqrt(2/n*).
double north_test_standard_error(double eigenvalue, double effective_dof);

Matrix temporal_correlation_matrix(const Matrix& data);

// Returns the rotated loadings; `rotation` receives the K x K rotation.
Matrix varimax_rotation(const Matrix& loadings, int max_iterations, double tolerance,
                        Matrix& rotation);

Matrix rotate_loadings(const Matrix& loadings, const Matrix& rotation);

std::vector<double> reof_rpc_trends(const Matrix& data, const Matrix& loadings,
                                    const Matrix& rotation);

double reof_back_projection(const std::vector<double>& trends, const Matrix& loadings,
                            const Matrix& rotation, std::size_t site);

// ---- Weatherhead detectability -----------------------------------------

double weatherhead_trend_standard_error(double noise_sd, double autocorrelation, double years);

double weatherhead_se_autocorrelation_derivative(double noise_sd, double autocorrelation,
                                                 double years);

double weatherhead_se_noise_derivative(double noise_sd, double autocorrelation,
                                       double years);

double minimum_detectable_trend(double noise_sd, double autocorrelation, double years,
                                double alpha);

double detection_time(double noise_sd, double autocorrelation, double trend_magnitude,
                      double alpha);

double detection_time_autocorrelation_derivative(double noise_sd, double autocorrelation,
                                                 double trend_magnitude, double alpha);

double detection_time_noise_derivative(double noise_sd, double autocorrelation,
                                       double trend_magnitude, double alpha);

// ---- Vertical coherence -------------------------------------------------

double trend_pattern_correlation(const std::vector<double>& first,
                                 const std::vector<double>& second);

double effective_spatial_dof(double points, double spatial_autocorrelation);

double correlation_t_statistic(double correlation, double effective_dof);

double coherence_height(const std::vector<double>& pressures,
                        const std::vector<double>& correlations, double critical_value);

double segmented_trend_difference(double second_period, double first_period);

// Point-wise Bayesian membership inference.
double bayesian_posterior(double likelihood_inside, double likelihood_outside,
                          double prior_probability);

double bayesian_log_odds(double likelihood_inside, double likelihood_outside,
                         double prior_probability);

double posterior_prior_derivative(double likelihood_inside, double likelihood_outside,
                                  double prior_probability);

double probability_weighted_centroid(const std::vector<double>& coordinates,
                                     const std::vector<double>& probabilities);

double centroid_probability_sensitivity(double coordinate, double centroid,
                                        double total_probability);

// One-sided Kolmogorov-Smirnov statistic sup(F_null - F_observed).
double kolmogorov_smirnov_d_plus(const std::vector<double>& observed,
                                 const std::vector<double>& null_sample);

// Chi-square optical-class trend precision.
double optical_trend_precision(double prior_precision, const std::vector<double>& covariate,
                               double measurement_variance, double process_variance);

// ---- Elevation-dependent warming ---------------------------------------

double elevation_trend_slope(const std::vector<double>& elevations,
                             const std::vector<double>& trends);

double elevation_leverage(double elevation, const std::vector<double>& elevations);

double latitude_corrected_edw(double raw_edw, double latitude_sensitivity,
                              double latitude_elevation_slope);

double lapse_rate_trend(double later_lapse_rate, double earlier_lapse_rate, double elapsed);

}  // namespace climate

#endif
