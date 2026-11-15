#include "intercal.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double gsics_leo_radiance(double geo_radiance, double offset, double slope) {
    require(std::fabs(slope) > 0.0, "gsics_leo_radiance: zero slope");
    return -offset / slope + geo_radiance / slope;
}

double gsics_uncertainty(double geo_radiance, double offset, double slope,
                         double offset_uncertainty, double slope_uncertainty,
                         double offset_slope_covariance) {
    require(std::fabs(slope) > 0.0, "gsics_uncertainty: zero slope");
    const double first = offset_uncertainty / slope;
    const double second = (geo_radiance - offset) * slope_uncertainty / (slope * slope);
    const double variance = first * first + second * second -
                            2.0 * (geo_radiance - offset) / (slope * slope * slope) *
                                offset_slope_covariance;
    require(variance >= 0.0, "gsics_uncertainty: negative variance");
    return std::sqrt(variance);
}

double msu_linear_radiance(double cold_space_radiance, double slope, double earth_counts,
                           double cold_counts) {
    return cold_space_radiance + slope * (earth_counts - cold_counts);
}

double msu_drift_term(double kappa, double mechanism_one, double mechanism_two) {
    return kappa * (mechanism_one + mechanism_two);
}

double msu_tlt_synthesis(double tmt, double tut, double tls) {
    return 1.430 * tmt - 0.462 * tut + 0.032 * tls;
}

double msu_delta_radiance(double alpha, double beta, double time, double reference_time) {
    return alpha + beta * (time - reference_time);
}

double msu_nonlinearity(double mu0, double k, double time, double reference_time) {
    return mu0 + k * (time - reference_time);
}

double msu_z_factor(double scene_parameter, double earth_emissivity, double cold_emissivity,
                    double warm_emissivity) {
    return scene_parameter * scene_parameter * (earth_emissivity - cold_emissivity) *
           (earth_emissivity - warm_emissivity);
}

double msu_calibrated_radiance(double reference_radiance, double delta_radiance,
                               double nonlinearity, double z_factor) {
    return reference_radiance - delta_radiance + nonlinearity * z_factor;
}

double diurnal_drift_correction(double observed, double factor, double drift_anomaly) {
    return observed + factor * drift_anomaly;
}

double double_difference(double target_observation, double reference_observation,
                         double target_model, double reference_model) {
    return (target_observation - target_model) - (reference_observation - reference_model);
}

double predicted_brightness(double reference_observation, double target_model,
                            double reference_model) {
    return reference_observation + target_model - reference_model;
}

double diurnal_anomaly(double alpha, double beta, double gamma, double local_solar_time) {
    const double angle = 2.0 * kPi * local_solar_time / 24.0;
    return alpha + beta * std::sin(angle) + gamma * std::cos(angle);
}

double seasonal_amplitude(double base, double sine_coefficient, double cosine_coefficient,
                          int month) {
    const double angle = 2.0 * kPi * static_cast<double>(month) / 12.0;
    return base + sine_coefficient * std::sin(angle) + cosine_coefficient * std::cos(angle);
}

double cwv_retrieval(const std::vector<double>& brightness_temperatures,
                     const std::vector<double>& coefficients) {
    require(coefficients.size() == brightness_temperatures.size() + 1,
            "cwv_retrieval: coefficient count");
    double value = coefficients[0];
    for (std::size_t i = 0; i < brightness_temperatures.size(); ++i) {
        require(brightness_temperatures[i] < 288.0, "cwv_retrieval: brightness above 288 K");
        value += coefficients[i + 1] * std::log((288.0 - brightness_temperatures[i]) / 288.0);
    }
    return value;
}

double cwv_channel_sensitivity(double brightness_temperature, double coefficient) {
    require(brightness_temperature < 288.0, "cwv_channel_sensitivity: brightness above 288 K");
    return -coefficient / (288.0 - brightness_temperature);
}

double sst_first_step(const std::vector<double>& brightness_temperatures,
                      const std::vector<double>& linear_coefficients,
                      const std::vector<double>& quadratic_coefficients) {
    require(linear_coefficients.size() == brightness_temperatures.size() + 1 &&
                quadratic_coefficients.size() == brightness_temperatures.size(),
            "sst_first_step: coefficient count");
    double value = linear_coefficients[0];
    for (std::size_t i = 0; i < brightness_temperatures.size(); ++i) {
        require(brightness_temperatures[i] < 288.0, "sst_first_step: brightness above 288 K");
        const double logarithm = std::log((288.0 - brightness_temperatures[i]) / 288.0);
        value += linear_coefficients[i + 1] * logarithm +
                 quadratic_coefficients[i] * logarithm * logarithm;
    }
    return value;
}

double sst_channel_sensitivity(double brightness_temperature, double linear_coefficient,
                               double quadratic_coefficient) {
    require(brightness_temperature < 288.0, "sst_channel_sensitivity: brightness above 288 K");
    const double logarithm = std::log((288.0 - brightness_temperature) / 288.0);
    return -(linear_coefficient + 2.0 * quadratic_coefficient * logarithm) /
           (288.0 - brightness_temperature);
}

double ssa_recalibration(double years_since_launch, double a0, double a1, double a2, double A1,
                         double A2, double A3, double A4) {
    return a0 + a1 * years_since_launch + a2 * years_since_launch * years_since_launch +
           A1 * std::sin(A2 * years_since_launch + A3) + A4;
}

double ssa_derivative(double years_since_launch, double a1, double a2, double A1, double A2,
                      double A3) {
    return a1 + 2.0 * a2 * years_since_launch +
           A1 * A2 * std::cos(A2 * years_since_launch + A3);
}

double clarreo_weight(double spatial_variance, double temporal_variance,
                      double noise_variance) {
    const double total = spatial_variance + temporal_variance + noise_variance;
    require(total > 0.0, "clarero_weight: zero total variance");
    return 1.0 / total;
}

double clarreo_bias(const std::vector<double>& differences,
                    const std::vector<double>& weights) {
    require(!differences.empty() && differences.size() == weights.size(),
            "clarero_bias: size mismatch");
    double weight_sum = 0.0, weighted = 0.0;
    for (std::size_t i = 0; i < differences.size(); ++i) {
        weighted += weights[i] * differences[i];
        weight_sum += weights[i];
    }
    require(weight_sum > 0.0, "clarero_bias: zero weight");
    return weighted / weight_sum;
}

double clarreo_uncertainty(const std::vector<double>& weights) {
    double total = 0.0;
    for (double weight : weights) {
        require(weight > 0.0, "clarero_uncertainty: non-positive weight");
        total += weight;
    }
    require(total > 0.0, "clarero_uncertainty: zero total weight");
    return 1.0 / std::sqrt(total);
}

double clarreo_temporal_variance(double overlap_fraction, double sounder_variance,
                                 double footprint_count) {
    require(footprint_count > 0.0, "clarero_temporal_variance: footprint count");
    require(overlap_fraction >= 0.0 && overlap_fraction <= 1.0,
            "clarero_temporal_variance: overlap fraction");
    return 2.0 * (1.0 - overlap_fraction) * sounder_variance / footprint_count;
}

double clarreo_inflate(double variance, double sounder_samples) {
    require(sounder_samples > 1.0, "clarero_inflate: need more than one sample");
    return variance * (1.0 + std::sqrt(2.0 / (sounder_samples - 1.0)));
}

double obcs_radiance(double calibration_coefficient, double nonlinearity_corrected_signal,
                     double offset, double background) {
    return calibration_coefficient * (nonlinearity_corrected_signal - offset) - background;
}

double angular_standard_deviation(const std::vector<double>& radiance_differences) {
    require(radiance_differences.size() >= 2, "angular_standard_deviation: need two values");
    double mean = 0.0;
    for (double value : radiance_differences) mean += value;
    mean /= static_cast<double>(radiance_differences.size());
    double variance = 0.0;
    for (double value : radiance_differences) variance += (value - mean) * (value - mean);
    variance /= static_cast<double>(radiance_differences.size());
    return std::sqrt(variance);
}

double bias_regression_slope(const std::vector<double>& target,
                             const std::vector<double>& reference) {
    require(!target.empty() && target.size() == reference.size(),
            "bias_regression_slope: size mismatch");
    double target_mean = 0.0, reference_mean = 0.0;
    for (std::size_t i = 0; i < target.size(); ++i) {
        target_mean += target[i];
        reference_mean += reference[i];
    }
    target_mean /= static_cast<double>(target.size());
    reference_mean /= static_cast<double>(reference.size());
    double covariance = 0.0, variance = 0.0;
    for (std::size_t i = 0; i < target.size(); ++i) {
        covariance += (target[i] - target_mean) * (reference[i] - reference_mean);
        variance += (reference[i] - reference_mean) * (reference[i] - reference_mean);
    }
    require(variance > 0.0, "bias_regression_slope: zero reference variance");
    return covariance / variance;
}

double msu_z_counts(double slope, double earth_counts, double cold_counts, double warm_counts) {
    return slope * slope * (earth_counts - cold_counts) * (earth_counts - warm_counts);
}

double aci_susceptibility(double droplet_slope, double extinction_slope) {
    return droplet_slope * extinction_slope;
}

double nd_relative_bias(double cloud_fraction_slope, double cloud_fraction_intercept,
                        double cloud_fraction) {
    return cloud_fraction_slope * cloud_fraction + cloud_fraction_intercept;
}

double lwp_error_variance(double effective_radius, double optical_depth,
                          double optical_depth_variance, double effective_radius_variance,
                          double cross_covariance) {
    return effective_radius * effective_radius * optical_depth_variance +
           optical_depth * optical_depth * effective_radius_variance +
           2.0 * effective_radius * optical_depth * cross_covariance;
}

double dml_treatment_effect(const std::vector<double>& residual_outcome,
                            const std::vector<double>& residual_treatment) {
    require(!residual_outcome.empty() && residual_outcome.size() == residual_treatment.size(),
            "dml_treatment_effect: size mismatch");
    double covariance = 0.0, variance = 0.0;
    for (std::size_t i = 0; i < residual_outcome.size(); ++i) {
        covariance += residual_treatment[i] * residual_outcome[i];
        variance += residual_treatment[i] * residual_treatment[i];
    }
    require(variance > 0.0, "dml_treatment_effect: zero treatment variance");
    return covariance / variance;
}

double g_computation_effect(const std::vector<double>& stage_effects,
                            const std::vector<double>& stage_probabilities) {
    require(stage_effects.size() == stage_probabilities.size(),
            "g_computation_effect: size mismatch");
    double total = 0.0;
    for (std::size_t i = 0; i < stage_effects.size(); ++i)
        total += stage_effects[i] * stage_probabilities[i];
    return total;
}

double natural_effect_total(double direct_effect, double indirect_effect) {
    return direct_effect + indirect_effect;
}

double mean_squared_error(const std::vector<double>& predicted,
                          const std::vector<double>& observed) {
    require(!predicted.empty() && predicted.size() == observed.size(),
            "mean_squared_error: size mismatch");
    double total = 0.0;
    for (std::size_t i = 0; i < predicted.size(); ++i) {
        const double difference = predicted[i] - observed[i];
        total += difference * difference;
    }
    return total / static_cast<double>(predicted.size());
}

double coefficient_of_determination(const std::vector<double>& predicted,
                                    const std::vector<double>& observed) {
    require(!predicted.empty() && predicted.size() == observed.size(),
            "coefficient_of_determination: size mismatch");
    double mean = 0.0;
    for (double value : observed) mean += value;
    mean /= static_cast<double>(observed.size());
    double residual = 0.0, total = 0.0;
    for (std::size_t i = 0; i < predicted.size(); ++i) {
        residual += (predicted[i] - observed[i]) * (predicted[i] - observed[i]);
        total += (observed[i] - mean) * (observed[i] - mean);
    }
    require(total > 0.0, "coefficient_of_determination: zero variance");
    return 1.0 - residual / total;
}

double sno_scene_variance(double instrument_variance, double scene_variance) {
    return instrument_variance + scene_variance;
}

double clarreo_bias_variance(const std::vector<double>& weights) {
    double total = 0.0;
    for (double weight : weights) {
        require(weight > 0.0, "clarero_bias_variance: non-positive weight");
        total += weight;
    }
    require(total > 0.0, "clarero_bias_variance: zero total weight");
    return 1.0 / total;
}

double truths_total_uncertainty(double spectral, double geometric, double noise) {
    return std::sqrt(spectral * spectral + geometric * geometric + noise * noise);
}

double kalman_random_walk(double mean, double variance, double process_variance, double dt,
                          double observation, double observation_variance, double& gain,
                          double& innovation, double& innovation_variance) {
    require(variance >= 0.0 && process_variance >= 0.0 && dt >= 0.0,
            "kalman_random_walk: negative variance");
    require(observation_variance > 0.0, "kalman_random_walk: observation variance");
    const double predicted_variance = variance + process_variance * dt;
    innovation_variance = predicted_variance + observation_variance;
    gain = predicted_variance / innovation_variance;
    innovation = observation - mean;
    const double updated = mean + gain * innovation;
    variance = (1.0 - gain) * predicted_variance;
    return updated;
}

double adversarial_loss(const std::vector<double>& real_scores,
                        const std::vector<double>& fake_scores) {
    require(!real_scores.empty() && !fake_scores.empty(), "adversarial_loss: empty scores");
    double total = 0.0;
    for (double score : real_scores) {
        require(score > 0.0 && score <= 1.0, "adversarial_loss: score outside (0,1]");
        total += std::log(score);
    }
    for (double score : fake_scores) {
        require(score >= 0.0 && score < 1.0, "adversarial_loss: score outside [0,1)");
        total += std::log(1.0 - score);
    }
    return total / static_cast<double>(real_scores.size() + fake_scores.size());
}

double cycle_consistency_loss(const std::vector<double>& forward_cycle,
                              const std::vector<double>& original) {
    require(!forward_cycle.empty() && forward_cycle.size() == original.size(),
            "cycle_consistency_loss: size mismatch");
    double total = 0.0;
    for (std::size_t i = 0; i < original.size(); ++i)
        total += std::fabs(forward_cycle[i] - original[i]);
    return total / static_cast<double>(original.size());
}

double cyclegan_objective(double forward_adversarial, double backward_adversarial,
                          double cycle_loss, double lambda) {
    return forward_adversarial + backward_adversarial + lambda * cycle_loss;
}

double sno_variance_decomposition(double spatial_variance, double temporal_variance,
                                  double noise_variance) {
    return spatial_variance + temporal_variance + noise_variance;
}

double sno_noise_fraction(double noise_variance, double total_variance) {
    require(total_variance > 0.0, "sno_noise_fraction: zero total variance");
    return 100.0 * noise_variance / total_variance;
}

double sno_t_statistic(double mean_bias, double ensemble_standard_deviation, int event_count) {
    require(event_count > 0 && ensemble_standard_deviation > 0.0,
            "sno_t_statistic: inputs");
    return mean_bias / (ensemble_standard_deviation / std::sqrt(static_cast<double>(event_count)));
}

}  // namespace climate
