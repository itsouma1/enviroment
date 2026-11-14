#ifndef CLIMATE_ATTRIB_HPP
#define CLIMATE_ATTRIB_HPP

#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// PRP: partial radiative perturbation
// ---------------------------------------------------------------------------
// One-sided feedback parameter lambda_x = dR/dx * dx/dTs.
double prp_feedback(double flux_perturbed, double flux_control, double surface_temperature_perturbed,
                    double surface_temperature_control);

// Two-sided average of the forward and reverse substitutions.
double prp_feedback_two_sided(double flux_forward_perturbed, double flux_forward_control,
                              double flux_reverse_perturbed, double flux_reverse_control,
                              double delta_temperature_forward, double delta_temperature_reverse);

// ---------------------------------------------------------------------------
// fSST: fixed-SST effective radiative forcing
// ---------------------------------------------------------------------------
double fSST_erf(const std::vector<double>& perturbed_imbalance,
                const std::vector<double>& control_imbalance);

double erf_decomposition(double instantaneous_forcing,
                         const std::vector<double>& rapid_adjustments);

// ---------------------------------------------------------------------------
// Cloud adjustment factor
// ---------------------------------------------------------------------------
// Bracket of the published decomposition; the reference-state isomorphic term
// is supplied by the caller because it is not closed in the literature.
double cloud_adjustment_factor(double land_cloud_fraction, double land_area_fraction,
                               double cloud_fraction_change, double surface_albedo_change,
                               double cloud_albedo, double surface_albedo,
                               double latitude_area_fraction, double reference_state_term);

// ---------------------------------------------------------------------------
// Extended eight-term surface temperature decomposition
// ---------------------------------------------------------------------------
// Adds ground heat flux and a closure/residual term to the six-term core.
double surface_temperature_eight_term(
    double incoming_shortwave, double albedo_change, double albedo,
    double incoming_shortwave_change, double surface_emissivity, double air_temperature,
    double air_emissivity_change, double latent_heat_of_vaporisation,
    double evapotranspiration_change, double air_density, double drag_coefficient,
    double surface_temperature, double air_temperature_change, double aerodynamic_resistance,
    double aerodynamic_resistance_change, double ground_heat_change, double residual,
    double stefan_boltzmann = 5.670374419e-8);

// ---------------------------------------------------------------------------
// CMIP6 SWV feedback validation targets
// ---------------------------------------------------------------------------
struct FeedbackTarget {
    double median;
    double spread;
};

FeedbackTarget stratospheric_swv_feedback_target();
FeedbackTarget full_column_swv_feedback_target();

}  // namespace climate

#endif
