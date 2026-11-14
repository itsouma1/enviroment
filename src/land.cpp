#include "land.hpp"

#include <cmath>
#include <stdexcept>

#include "common.hpp"

namespace climate {
namespace {

constexpr double GAS_CONSTANT = 8.314462618;  // J mol^-1 K^-1
constexpr double KELVIN_OFFSET = 273.15;

}  // namespace

double land_net_radiation(double albedo, double shortwave_down, double emissivity_surface,
                          double emissivity_air, double air_temperature_k,
                          double surface_temperature_k) {
    if (albedo < 0.0 || albedo > 1.0)
        throw std::runtime_error("land_net_radiation: albedo must be in [0,1]");
    if (air_temperature_k <= 0.0 || surface_temperature_k <= 0.0)
        throw std::runtime_error("land_net_radiation: temperatures must be positive");
    const double ta4 = std::pow(air_temperature_k, 4.0);
    const double ts4 = std::pow(surface_temperature_k, 4.0);
    return (1.0 - albedo) * shortwave_down +
           emissivity_surface * STEFAN_BOLTZMANN * (emissivity_air * ta4 - ts4);
}

double sensible_heat_flux(double air_density, double specific_heat, double surface_temperature_k,
                          double air_temperature_k, double aerodynamic_resistance) {
    if (aerodynamic_resistance <= 0.0)
        throw std::runtime_error("sensible_heat_flux: aerodynamic resistance must be positive");
    return air_density * specific_heat * (surface_temperature_k - air_temperature_k) /
           aerodynamic_resistance;
}

double aerodynamic_resistance(double von_karman, double friction_velocity, double measurement_height,
                              double displacement_height, double scalar_roughness_length,
                              double stability_correction) {
    if (von_karman <= 0.0 || friction_velocity <= 0.0)
        throw std::runtime_error("aerodynamic_resistance: k and u* must be positive");
    if (scalar_roughness_length <= 0.0)
        throw std::runtime_error("aerodynamic_resistance: z_oh must be positive");
    const double height = measurement_height - displacement_height;
    if (height <= scalar_roughness_length)
        throw std::runtime_error("aerodynamic_resistance: measurement height must exceed d + z_oh");
    return (std::log(height / scalar_roughness_length) - stability_correction) /
           (von_karman * friction_velocity);
}

namespace {

double penman_monteith_denominator(double delta, double psychrometric_constant,
                                   double surface_resistance, double aerodynamic_resistance) {
    if (aerodynamic_resistance <= 0.0)
        throw std::runtime_error("penman_monteith: aerodynamic resistance must be positive");
    return delta + psychrometric_constant * (1.0 + surface_resistance / aerodynamic_resistance);
}

}  // namespace

double penman_monteith_latent_heat(double delta, double net_radiation_minus_ground,
                                   double air_density, double specific_heat,
                                   double humidity_deficit, double aerodynamic_resistance,
                                   double psychrometric_constant, double surface_resistance) {
    const double denominator =
        penman_monteith_denominator(delta, psychrometric_constant, surface_resistance,
                                    aerodynamic_resistance);
    if (denominator == 0.0)
        throw std::runtime_error("penman_monteith_latent_heat: zero denominator");
    return (delta * net_radiation_minus_ground +
            air_density * specific_heat * humidity_deficit / aerodynamic_resistance) /
           denominator;
}

double penman_monteith_sensible_heat(double delta, double net_radiation_minus_ground,
                                     double air_density, double specific_heat,
                                     double humidity_deficit, double aerodynamic_resistance,
                                     double psychrometric_constant, double surface_resistance) {
    if (aerodynamic_resistance <= 0.0)
        throw std::runtime_error("penman_monteith_sensible_heat: r_a must be positive");
    const double denominator =
        penman_monteith_denominator(delta, psychrometric_constant, surface_resistance,
                                    aerodynamic_resistance);
    if (denominator == 0.0)
        throw std::runtime_error("penman_monteith_sensible_heat: zero denominator");
    const double resistance_ratio = 1.0 + surface_resistance / aerodynamic_resistance;
    return (net_radiation_minus_ground * psychrometric_constant * resistance_ratio -
            air_density * specific_heat * humidity_deficit / aerodynamic_resistance) /
           denominator;
}

double surface_resistance_soil(double minimum_resistance, double saturation_moisture,
                               double residual_moisture, double top_layer_moisture) {
    if (minimum_resistance < 0.0)
        throw std::runtime_error("surface_resistance_soil: resistance must be >= 0");
    const double denominator = top_layer_moisture - residual_moisture;
    if (denominator <= 0.0)
        throw std::runtime_error("surface_resistance_soil: top-layer moisture must exceed residual");
    const double ratio = (saturation_moisture - residual_moisture) / denominator;
    return minimum_resistance * ratio * ratio * ratio;
}

double farquhar_rubisco(double vcmax, double intercellular_co2, double compensation_point,
                        double kc, double intercellular_o2, double ko) {
    if (ko == 0.0) throw std::runtime_error("farquhar_rubisco: Ko must be non-zero");
    const double denominator = intercellular_co2 + kc * (1.0 + intercellular_o2 / ko);
    if (denominator == 0.0) throw std::runtime_error("farquhar_rubisco: zero denominator");
    return vcmax * (intercellular_co2 - compensation_point) / denominator;
}

double farquhar_electron_transport(double electron_transport, double intercellular_co2,
                                   double compensation_point) {
    const double denominator = 4.0 * intercellular_co2 + 8.0 * compensation_point;
    if (denominator == 0.0)
        throw std::runtime_error("farquhar_electron_transport: zero denominator");
    return electron_transport * (intercellular_co2 - compensation_point) / denominator;
}

double farquhar_assimilation(double rubisco_rate, double electron_transport_rate,
                             double dark_respiration) {
    return (rubisco_rate < electron_transport_rate ? rubisco_rate : electron_transport_rate) -
           dark_respiration;
}

double ball_berry_conductance(double g1, double assimilation, double atmospheric_pressure,
                              double surface_co2, double compensation_point, double vpd,
                              double vpd_scale) {
    if (vpd_scale <= 0.0) throw std::runtime_error("ball_berry: D0 must be positive");
    const double denominator = (surface_co2 - compensation_point) * (1.0 + vpd / vpd_scale);
    if (denominator == 0.0) throw std::runtime_error("ball_berry: zero denominator");
    return g1 * assimilation * atmospheric_pressure / denominator;
}

double medlyn_conductance(double g1, double assimilation, double surface_co2, double vpd) {
    if (vpd <= 0.0) throw std::runtime_error("medlyn: VPD must be positive");
    if (surface_co2 == 0.0) throw std::runtime_error("medlyn: Cs must be non-zero");
    return (1.0 + g1 / std::sqrt(vpd)) * assimilation / surface_co2;
}

double leuning_conductance(double cuticular_conductance, double a1, double assimilation,
                           double leaf_internal_co2, double compensation_point,
                           double vpd_response, double pressure_scaling) {
    const double denominator = leaf_internal_co2 - compensation_point;
    if (denominator == 0.0) throw std::runtime_error("leuning: Ct must differ from Gamma*");
    return cuticular_conductance +
           a1 * assimilation / denominator * vpd_response * pressure_scaling;
}

double npp_lpj(double gpp, double maintenance_respiration) {
    return 0.75 * (gpp - maintenance_respiration);
}

double maintenance_respiration(double respiration_coefficient, double q10, double temperature_c,
                               double leaf_nitrogen, double leaf_carbon, double root_carbon,
                               double sapwood_coefficient, double wood_leaf_carbon,
                               double wood_root_carbon) {
    if (q10 <= 0.0) throw std::runtime_error("maintenance_respiration: Q10 must be positive");
    const double temperature_factor = std::pow(q10, (temperature_c - 20.0) / 10.0);
    const double pools = (leaf_carbon + root_carbon) +
                         sapwood_coefficient * (wood_leaf_carbon + wood_root_carbon);
    return respiration_coefficient * temperature_factor * leaf_nitrogen * pools;
}

double lai_from_leaf_carbon(double leaf_carbon, double specific_leaf_area) {
    if (leaf_carbon < 0.0 || specific_leaf_area < 0.0)
        throw std::runtime_error("lai_from_leaf_carbon: inputs must be non-negative");
    return leaf_carbon * specific_leaf_area;
}

double xylem_conductance(double soil_water_potential, double p50, double shape) {
    if (shape <= 0.0) throw std::runtime_error("xylem_conductance: shape must be positive");
    const double exponent = shape * (soil_water_potential - p50);
    return 1.0 - 1.0 / (1.0 + std::exp(exponent));
}

double soil_water_step(double storage, double precipitation, double evaporation, double runoff,
                       double drainage, double dt) {
    if (dt < 0.0) throw std::runtime_error("soil_water_step: dt must be non-negative");
    return storage + (precipitation - evaporation - runoff - drainage) * dt;
}

double soil_thermal_inertia(double dry_soil_inertia, double volumetric_moisture,
                            double residual_moisture, double thermal_scaling) {
    return dry_soil_inertia + (volumetric_moisture - residual_moisture) * thermal_scaling;
}

double vcmax_temperature(double vcmax_25, double temperature_c, double activation_energy) {
    const double temperature_k = temperature_c + KELVIN_OFFSET;
    if (temperature_k <= 0.0)
        throw std::runtime_error("vcmax_temperature: temperature must be above absolute zero");
    const double reference = 298.15;  // 25 degC
    return vcmax_25 *
           std::exp(activation_energy * (temperature_k - reference) /
                    (reference * GAS_CONSTANT * temperature_k));
}

double leaf_turnover(double specific_leaf_area, double leaf_biomass) {
    if (specific_leaf_area <= 0.0)
        throw std::runtime_error("leaf_turnover: specific leaf area must be positive");
    if (leaf_biomass < 0.0) throw std::runtime_error("leaf_turnover: biomass must be >= 0");
    const double longevity = 29664.0 * std::pow(specific_leaf_area * 10.0, -0.909);
    return leaf_biomass / longevity;
}

double root_uptake(double top_fraction, double surface_potential, double root_fraction,
                   double root_potential) {
    if (top_fraction < 0.0 || root_fraction < 0.0)
        throw std::runtime_error("root_uptake: fractions must be >= 0");
    return top_fraction * surface_potential + root_fraction * root_potential;
}

double transpiration_moisture(double moisture, double saturation, double exponent,
                              double canopy_fraction, double potential) {
    if (!(saturation > 0.0)) throw std::runtime_error("transpiration_moisture: saturation");
    if (moisture < 0.0) throw std::runtime_error("transpiration_moisture: moisture >= 0");
    if (canopy_fraction < 0.0) throw std::runtime_error("transpiration_moisture: canopy >= 0");
    return std::pow(moisture / saturation, exponent) * canopy_fraction * potential;
}

double optimal_stomatal_objective(double assimilation, double transpiration,
                                  double critical_transpiration) {
    if (!(critical_transpiration > 0.0))
        throw std::runtime_error("optimal_stomatal_objective: critical transpiration");
    return assimilation * (1.0 - transpiration / critical_transpiration);
}

double hydraulic_profit(double assimilation, double jmax, double carbon_cost,
                        double hydraulic_cost, double water_potential_drop) {
    return assimilation - carbon_cost * jmax -
           hydraulic_cost * water_potential_drop * water_potential_drop;
}

double growth_temperature_acclimation(double intercept, double slope,
                                      double growth_temperature_c) {
    return intercept + slope * growth_temperature_c;
}

double mesophyll_limitation(double assimilation_ambient, double assimilation_intercellular,
                            double mesophyll_conductance) {
    const double denominator = assimilation_intercellular - mesophyll_conductance;
    if (!(std::fabs(denominator) > 0.0))
        throw std::runtime_error("mesophyll_limitation: zero denominator");
    return (1.0 - (assimilation_ambient - mesophyll_conductance) / denominator) * 100.0;
}

double nutrient_colimitation(double nitrogen_limitation, double phosphorus_limitation) {
    return std::min(nitrogen_limitation, phosphorus_limitation);
}

double nutrient_demand(double net_carbon, double stoichiometric_ratio, double retranslocation) {
    if (!(stoichiometric_ratio > 0.0))
        throw std::runtime_error("nutrient_demand: stoichiometric ratio must be positive");
    return std::max(net_carbon / stoichiometric_ratio - retranslocation, 0.0);
}

double nitrogen_fixation(double rate, double density, double mass, double dt) {
    if (density < 0.0 || mass < 0.0 || dt < 0.0)
        throw std::runtime_error("nitrogen_fixation: negative input");
    return rate * density * mass * dt;
}

double canopy_air_exchange(double air_density, double canopy_specific_humidity,
                           double air_specific_humidity, double resistance) {
    if (!(resistance > 0.0)) throw std::runtime_error("canopy_air_exchange: resistance");
    return -air_density * (canopy_specific_humidity - air_specific_humidity) / resistance;
}

double multilayer_scalar_tendency(double flux_up, double flux_down, double source,
                                  double layer_thickness) {
    if (!(layer_thickness > 0.0))
        throw std::runtime_error("multilayer_scalar_tendency: thickness");
    return (-(flux_up - flux_down) + source * layer_thickness) / layer_thickness;
}

double soil_production(double reference_rate, double regolith_thickness, double decay_depth) {
    if (!(decay_depth > 0.0)) throw std::runtime_error("soil_production: decay depth");
    if (regolith_thickness < 0.0)
        throw std::runtime_error("soil_production: thickness must be >= 0");
    return reference_rate * std::exp(-regolith_thickness / decay_depth);
}

double soil_production_reference(double coefficient, double runoff, double activation_energy,
                                 double temperature_k, double reference_temperature) {
    if (!(temperature_k > 0.0))
        throw std::runtime_error("soil_production_reference: temperature must be positive");
    return coefficient * runoff *
           std::exp(-activation_energy / GAS_CONSTANT *
                    (1.0 / temperature_k - 1.0 / reference_temperature));
}

}  // namespace climate
