#ifndef CLIMATE_TEMPAC_HPP
#define CLIMATE_TEMPAC_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

// First-harmonic amplitude of a 12-month annual cycle.
double first_harmonic_amplitude(const std::vector<double>& monthly_values);

std::vector<double> first_harmonic_amplitudes(const Matrix& monthly_by_site);

// First EOF of a time x site anomaly field (unweighted covariance).
std::vector<double> first_eof(const Matrix& field);

std::vector<double> project_fingerprint(const Matrix& field,
                                        const std::vector<double>& fingerprint);

std::vector<double> fingerprint_trends(const Matrix& field,
                                       const std::vector<double>& fingerprint,
                                       std::size_t window);

double detection_p_value(double observed_trend, std::vector<double> noise_trends);

double pattern_signal_to_noise(double signal_trend, const std::vector<double>& noise_trends);

double local_signal_to_noise(double ensemble_trend, double trend_spread);

// Ratio of S/N at two analysis lengths; the noise spread must contract with length.
double signal_to_noise_length_ratio(double short_spread, double long_spread);

// Cyclostationary harmonic loading of a series at a given period and harmonic.
std::vector<double> cseof_harmonic_loading(const std::vector<double>& series, double period,
                                           int harmonic, double& amplitude, double& phase);

// Reduced-space optimal interpolation: x_a = x_f + K (y - H x_f).
std::vector<double> reduced_space_analysis(const std::vector<double>& forecast,
                                           const std::vector<double>& observations,
                                           const Matrix& observation_operator,
                                           const Matrix& background_covariance,
                                           const Matrix& observation_covariance);

}  // namespace climate

#endif
