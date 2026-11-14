#include "attrib.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double mean(const std::vector<double>& values) {
    require(!values.empty(), "mean: empty series");
    double total = 0.0;
    for (double value : values) total += value;
    return total / static_cast<double>(values.size());
}

}  // namespace

double prp_feedback(double flux_perturbed, double flux_control,
                    double surface_temperature_perturbed, double surface_temperature_control) {
    const double denominator = surface_temperature_perturbed - surface_temperature_control;
    require(std::fabs(denominator) > 0.0, "prp_feedback: zero temperature difference");
    return (flux_perturbed - flux_control) / denominator;
}

double prp_feedback_two_sided(double flux_forward_perturbed, double flux_forward_control,
                              double flux_reverse_perturbed, double flux_reverse_control,
                              double delta_temperature_forward, double delta_temperature_reverse) {
    require(std::fabs(delta_temperature_forward) > 0.0 &&
                std::fabs(delta_temperature_reverse) > 0.0,
            "prp_feedback_two_sided: zero temperature difference");
    const double forward = (flux_forward_perturbed - flux_forward_control) /
                           delta_temperature_forward;
    const double reverse = (flux_reverse_perturbed - flux_reverse_control) /
                           delta_temperature_reverse;
    return 0.5 * (forward + reverse);
}

double fSST_erf(const std::vector<double>& perturbed_imbalance,
                const std::vector<double>& control_imbalance) {
    require(!perturbed_imbalance.empty() && perturbed_imbalance.size() == control_imbalance.size(),
            "fSST_erf: size mismatch");
    return mean(perturbed_imbalance) - mean(control_imbalance);
}

double erf_decomposition(double instantaneous_forcing,
                         const std::vector<double>& rapid_adjustments) {
    double total = instantaneous_forcing;
    for (double adjustment : rapid_adjustments) total += adjustment;
    return total;
}

double cloud_adjustment_factor(double land_cloud_fraction, double land_area_fraction,
                               double cloud_fraction_change, double surface_albedo_change,
                               double cloud_albedo, double surface_albedo,
                               double latitude_area_fraction, double reference_state_term) {
    require(land_area_fraction > 0.0, "cloud_adjustment_factor: zero global area");
    require(reference_state_term > 0.0,
            "cloud_adjustment_factor: reference-state term must be positive");
    const double cloud_fraction_derivative =
        std::fabs(surface_albedo_change) > 0.0
            ? cloud_fraction_change / surface_albedo_change
            : 0.0;
    const double bracket =
        (1.0 - land_cloud_fraction) * land_area_fraction +
        cloud_fraction_derivative * (cloud_albedo - surface_albedo) *
            latitude_area_fraction;
    return bracket / reference_state_term;
}

double surface_temperature_eight_term(
    double incoming_shortwave, double albedo_change, double albedo,
    double incoming_shortwave_change, double surface_emissivity, double air_temperature,
    double air_emissivity_change, double latent_heat_of_vaporisation,
    double evapotranspiration_change, double air_density, double drag_coefficient,
    double surface_temperature, double air_temperature_change, double aerodynamic_resistance,
    double aerodynamic_resistance_change, double ground_heat_change, double residual,
    double stefan_boltzmann) {
    require(aerodynamic_resistance > 0.0,
            "surface_temperature_eight_term: aerodynamic resistance");
    const double turbulent = air_density * drag_coefficient / aerodynamic_resistance;
    const double surface_radiative =
        4.0 * surface_emissivity * stefan_boltzmann * surface_temperature * surface_temperature *
        surface_temperature;
    const double factor = turbulent + surface_radiative;
    require(std::fabs(factor) > 0.0,
            "surface_temperature_eight_term: zero energy redistribution factor");
    const double bracket =
        -incoming_shortwave * albedo_change + (1.0 - albedo) * incoming_shortwave_change +
        surface_emissivity * stefan_boltzmann * air_temperature * air_temperature *
            air_temperature * air_temperature * air_emissivity_change -
        latent_heat_of_vaporisation * evapotranspiration_change +
        air_density * drag_coefficient * (surface_temperature - air_temperature) /
            (aerodynamic_resistance * aerodynamic_resistance) * aerodynamic_resistance_change -
        ground_heat_change + residual;
    const double air_term =
        (turbulent +
         4.0 * surface_emissivity * stefan_boltzmann * air_temperature * air_temperature *
             air_temperature) /
        factor * air_temperature_change;
    return bracket / factor + air_term;
}

FeedbackTarget stratospheric_swv_feedback_target() { return {0.25, 0.06}; }

FeedbackTarget full_column_swv_feedback_target() { return {0.020, 0.007}; }

}  // namespace climate
