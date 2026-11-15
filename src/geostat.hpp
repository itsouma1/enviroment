#ifndef CLIMATE_GEOSTAT_HPP
#define CLIMATE_GEOSTAT_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// ---------------------------------------------------------------------------
// SPDE / GMRF
// ---------------------------------------------------------------------------
// Matern variance sigma^2 = Gamma(nu) / (Gamma(alpha) (4 pi)^(d/2) kappa^(2 nu) tau^2),
// with nu = alpha - d/2.
double spde_matern_variance(double kappa, double tau, int alpha, int dimension);

// Q = tau^2 (kappa^2 C + G)^alpha for integer alpha.
Matrix spde_precision(const Matrix& mass, const Matrix& stiffness, double kappa, double tau,
                      int alpha);

// Conditional mean Q^-1 eta and covariance Q^-1.
std::vector<double> gmrf_conditional_mean(const Matrix& precision,
                                         const std::vector<double>& eta);

Matrix gmrf_conditional_covariance(const Matrix& precision);

// ---------------------------------------------------------------------------
// Gaussian process regression and MAD envelope
// ---------------------------------------------------------------------------
double squared_exponential_kernel(double distance, double variance, double length_scale);

std::vector<double> gp_posterior_mean(const std::vector<double>& train_x,
                                     const std::vector<double>& train_y,
                                     const std::vector<double>& predict_x,
                                     double signal_variance, double length_scale,
                                     double noise_variance);

std::vector<double> gp_posterior_variance(const std::vector<double>& train_x,
                                         const std::vector<double>& predict_x,
                                         double signal_variance, double length_scale,
                                         double noise_variance);

// max_k |T_ik - T_0k| / s_k for one ensemble member.
double mad_envelope_statistic(const std::vector<double>& observed,
                              const std::vector<double>& ensemble_member,
                              const std::vector<double>& pointwise_scale);

// Scaled MAD statistic of one function against an ensemble mean and sd.
double mad_global_statistic(const std::vector<double>& function,
                            const std::vector<double>& ensemble_mean,
                            const std::vector<double>& ensemble_scale);

// Proper scaled MAD global envelope p-value: standardise the observed
// function and every ensemble member with the ensemble mean and sd.
double mad_envelope_p_value(const std::vector<double>& observed,
                            const std::vector<std::vector<double>>& ensemble);

// Overlap of two equal circles with centre separation d and radius r.
double circle_overlap_fraction(double separation, double radius);

// ---------------------------------------------------------------------------
// GTWR
// ---------------------------------------------------------------------------
double gtwr_distance(double u1, double v1, double t1, double u2, double v2, double t2,
                     double spatial_scale, double temporal_scale);

std::vector<double> gtwr_coefficients(const Matrix& predictors, const std::vector<double>& response,
                                      const std::vector<double>& u, const std::vector<double>& v,
                                      const std::vector<double>& t, std::size_t target,
                                      double spatial_scale, double temporal_scale,
                                      double bandwidth);

double gtwr_aicc(int sample_count, double residual_sigma, double hat_trace);

// ---------------------------------------------------------------------------
// Variograms
// ---------------------------------------------------------------------------
double spherical_variogram(double distance, double nugget, double sill, double range);

double exponential_variogram(double distance, double nugget, double sill, double range);

double gaussian_variogram(double distance, double nugget, double sill, double range);

// Matern variogram for nu = 1/2, 3/2 and 5/2.
double matern_variogram(double distance, double nugget, double sill, double range,
                        double smoothness);

// ---------------------------------------------------------------------------
// Weighted PCA
// ---------------------------------------------------------------------------
Matrix weighted_covariance(const Matrix& residuals, const std::vector<double>& weights);

// Project out the leading components of the weighted covariance.
Matrix weighted_pca_residual(const Matrix& residuals, const std::vector<double>& weights,
                             int components);

}  // namespace climate

#endif
