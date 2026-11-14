#include "cryo.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {
constexpr double kGasConstant = 8.314462618;      // J mol^-1 K^-1
constexpr double kStefanBoltzmann = 5.670374419e-8;

void require_positive(double value, const char* name) {
    if (!(value > 0.0)) throw std::runtime_error(std::string(name) + " must be positive");
}
}  // namespace

double glen_rate_factor(double reference_rate, double activation_energy,
                        double temperature_k) {
    require_positive(temperature_k, "glen_rate_factor: temperature");
    require_positive(reference_rate, "glen_rate_factor: reference rate");
    return reference_rate * std::exp(-activation_energy / (kGasConstant * temperature_k));
}

double glen_strain_rate(double rate_factor, double effective_stress, double n) {
    require_positive(rate_factor, "glen_strain_rate: rate factor");
    if (effective_stress < 0.0)
        throw std::runtime_error("glen_strain_rate: effective stress must be >= 0");
    return rate_factor * std::pow(effective_stress, n);
}

double effective_strain(double e_xx, double e_yy, double e_zz, double e_xy, double e_xz,
                        double e_yz) {
    const double invariant = e_xx * e_xx + e_yy * e_yy + e_zz * e_zz +
                             2.0 * (e_xy * e_xy + e_xz * e_xz + e_yz * e_yz);
    return std::sqrt(0.5 * invariant);
}

double effective_stress(double t_xx, double t_yy, double t_zz, double t_xy, double t_xz,
                        double t_yz) {
    const double invariant = t_xx * t_xx + t_yy * t_yy + t_zz * t_zz +
                             2.0 * (t_xy * t_xy + t_xz * t_xz + t_yz * t_yz);
    return std::sqrt(0.5 * invariant);
}

double ice_viscosity(double hardness, double effective_strain, double n) {
    require_positive(hardness, "ice_viscosity: hardness");
    require_positive(effective_strain, "ice_viscosity: effective strain");
    return hardness / (2.0 * std::pow(effective_strain, n - 1.0));
}

double sia_velocity(double coefficient, double surface_slope, double thickness, double n) {
    require_positive(thickness, "sia_velocity: thickness");
    return coefficient * surface_slope * std::pow(thickness, n + 1.0);
}

double weertman_sliding(double coefficient, double basal_shear, double effective_pressure,
                        double p, double q) {
    require_positive(coefficient, "weertman_sliding: coefficient");
    require_positive(effective_pressure, "weertman_sliding: effective pressure");
    if (basal_shear < 0.0)
        throw std::runtime_error("weertman_sliding: basal shear must be >= 0");
    return coefficient * std::pow(basal_shear, p) / std::pow(effective_pressure, q);
}

double coulomb_sliding(double friction, double effective_pressure) {
    return friction * effective_pressure;
}

HorizontalStress bulk_ice_stress(double density, double drag, double fluid_u,
                                 double fluid_v, double ice_u, double ice_v,
                                 double rotation) {
    require_positive(density, "bulk_ice_stress: density");
    if (drag < 0.0) throw std::runtime_error("bulk_ice_stress: drag must be >= 0");
    const double relative_u = fluid_u - ice_u;
    const double relative_v = fluid_v - ice_v;
    const double speed = std::hypot(relative_u, relative_v);
    const double magnitude = density * drag * speed;
    const double cosine = std::cos(rotation);
    const double sine = std::sin(rotation);
    HorizontalStress stress;
    stress.x = magnitude * (cosine * relative_u - sine * relative_v);
    stress.y = magnitude * (sine * relative_u + cosine * relative_v);
    return stress;
}

void sea_ice_momentum_tendency(double& du_dt, double& dv_dt, double u, double v,
                               double mass, double coriolis,
                               const HorizontalStress& air,
                               const HorizontalStress& ocean,
                               double potential_gradient_x, double potential_gradient_y,
                               const HorizontalStress& internal) {
    require_positive(mass, "sea_ice_momentum_tendency: mass");
    // m Du/Dt = -m f k x u + tau_air + tau_ocean - m grad(phi) + F
    // k x u = (-v, u), so -m f k x u = (m f v, -m f u).
    du_dt = (mass * coriolis * v + air.x + ocean.x - mass * potential_gradient_x + internal.x) /
            mass;
    dv_dt = (-mass * coriolis * u + air.y + ocean.y - mass * potential_gradient_y + internal.y) /
            mass;
}

double surface_melt_energy(double albedo, double shortwave, double longwave,
                           double emissivity, double surface_temperature) {
    if (surface_temperature < 0.0)
        throw std::runtime_error("surface_melt_energy: temperature must be in kelvin");
    return (1.0 - albedo) * shortwave + longwave -
           emissivity * kStefanBoltzmann * std::pow(surface_temperature, 4.0);
}

double stefan_front_velocity(double frozen_conductivity, double frozen_gradient,
                             double thawed_conductivity, double thawed_gradient,
                             double density, double latent_heat) {
    require_positive(density, "stefan_front_velocity: density");
    require_positive(latent_heat, "stefan_front_velocity: latent heat");
    return (frozen_conductivity * frozen_gradient - thawed_conductivity * thawed_gradient) /
           (density * latent_heat);
}

double stefan_number(double bulk_heat_capacity, double bulk_density,
                     double surface_temperature, double freezing_temperature,
                     double water_saturation, double freeze_fraction, double water_density,
                     double porosity, double latent_heat) {
    const double denominator =
        water_saturation * freeze_fraction * water_density * porosity * latent_heat;
    require_positive(denominator, "stefan_number: latent denominator");
    return bulk_heat_capacity * bulk_density * (surface_temperature - freezing_temperature) /
           denominator;
}

double active_layer_depth(double scaling, double conductivity, double thawing_degree_days,
                          double latent_heat) {
    require_positive(scaling, "active_layer_depth: scaling");
    require_positive(conductivity, "active_layer_depth: conductivity");
    require_positive(latent_heat, "active_layer_depth: latent heat");
    if (thawing_degree_days < 0.0)
        throw std::runtime_error("active_layer_depth: thawing degree days must be >= 0");
    return std::sqrt(2.0 * scaling * conductivity * thawing_degree_days / latent_heat);
}

double basal_melt_three_equation(double ocean_density, double ice_density, double drag,
                                 double transfer_coefficient, double latent_heat,
                                 double ice_heat_capacity, double ice_temperature,
                                 double freezing_temperature, double ocean_heat_capacity,
                                 double current_speed, double boundary_temperature) {
    require_positive(ocean_density, "basal_melt_three_equation: ocean density");
    require_positive(ice_density, "basal_melt_three_equation: ice density");
    require_positive(drag, "basal_melt_three_equation: drag");
    require_positive(latent_heat, "basal_melt_three_equation: latent heat");
    require_positive(ocean_heat_capacity, "basal_melt_three_equation: ocean heat capacity");
    const double effective_latent =
        (latent_heat - ice_heat_capacity * (ice_temperature - freezing_temperature)) /
        ocean_heat_capacity;
    if (!(std::fabs(effective_latent) > 0.0))
        throw std::runtime_error("basal_melt_three_equation: zero effective latent heat");
    return (ocean_density / ice_density) * (std::sqrt(drag) * transfer_coefficient /
                                            effective_latent) *
           current_speed * (boundary_temperature - freezing_temperature);
}

double basal_melt_simplified(double ocean_heat_capacity, double latent_heat,
                             double ocean_density, double ice_density,
                             double transfer_coefficient, double mixing_efficiency,
                             double local_temperature, double local_freezing) {
    require_positive(latent_heat, "basal_melt_simplified: latent heat");
    require_positive(ice_density, "basal_melt_simplified: ice density");
    return (ocean_heat_capacity / latent_heat) * (ocean_density / ice_density) *
           transfer_coefficient * mixing_efficiency * (local_temperature - local_freezing);
}

double plume_melt(double ambient_temperature, double grounding_line_freezing,
                  double melt_function, double haline_contraction, double ambient_salinity,
                  double gravity, double exchange_velocity, double lambda3, double plume_length,
                  double density_ratio, double gamma, double drag, double entrainment,
                  double sin_slope, double tau_coefficient) {
    require_positive(gravity, "plume_melt: gravity");
    require_positive(lambda3, "plume_melt: lambda3");
    require_positive(plume_length, "plume_melt: plume length");
    require_positive(exchange_velocity, "plume_melt: exchange velocity");
    const double buoyancy = haline_contraction * ambient_salinity * gravity *
                            std::pow(exchange_velocity, 3.0) /
                            (lambda3 * std::pow(plume_length, 3.0));
    const double first_root =
        std::sqrt((1.0 - density_ratio * gamma) / (drag + entrainment * sin_slope));
    const double second_root = std::pow(gamma * entrainment * sin_slope /
                                            (gamma + tau_coefficient + entrainment * sin_slope),
                                        1.5);
    const double delta = ambient_temperature - grounding_line_freezing;
    return delta * delta * melt_function * std::sqrt(buoyancy) * first_root * second_root;
}

double subglacial_flux(double coefficient, double gap_height, double alpha, double beta,
                       double head_gradient) {
    require_positive(coefficient, "subglacial_flux: coefficient");
    if (gap_height < 0.0) throw std::runtime_error("subglacial_flux: gap must be >= 0");
    const double magnitude = std::pow(std::fabs(head_gradient), beta - 1.0);
    return -coefficient * std::pow(gap_height, alpha) * magnitude * head_gradient;
}

double hydraulic_potential(double water_density, double gravity, double bed_elevation,
                           double water_pressure) {
    return water_density * gravity * bed_elevation + water_pressure;
}

double subglacial_melt(double geothermal_flux, double basal_shear, double basal_speed,
                       double discharge, double head_gradient, double water_heat_capacity,
                       double water_density, double pressure_gradient, double ice_density,
                       double latent_heat) {
    require_positive(ice_density, "subglacial_melt: ice density");
    require_positive(latent_heat, "subglacial_melt: latent heat");
    const double energy = geothermal_flux + basal_shear * basal_speed +
                          discharge * head_gradient -
                          water_heat_capacity * water_density * discharge * pressure_gradient;
    return energy / (ice_density * latent_heat);
}

double calving_flux(double calving_speed, double cross_section) {
    if (cross_section < 0.0)
        throw std::runtime_error("calving_flux: cross section must be >= 0");
    return calving_speed * cross_section;
}

double phase_field_notch(double thickness, double tensile_strength) {
    if (thickness < 0.0)
        throw std::runtime_error("phase_field_notch: thickness must be >= 0");
    return 0.33 * std::sqrt(3.75 * thickness + tensile_strength * tensile_strength);
}

double firn_densification_rate(double overburden, double grain_radius, double ice_density,
                               double firn_density, double activation_energy,
                               double temperature_k) {
    require_positive(grain_radius, "firn_densification_rate: grain radius");
    require_positive(temperature_k, "firn_densification_rate: temperature");
    require_positive(ice_density, "firn_densification_rate: ice density");
    if (firn_density < 0.0 || firn_density > ice_density)
        throw std::runtime_error("firn_densification_rate: density outside [0, rho_ice]");
    const double kc = firn_density < 550.0 ? 9.2e-9 : 3.7e-9;
    return kc * overburden * std::pow(grain_radius, -2.0) * (ice_density - firn_density) *
           std::exp(-activation_energy / (kGasConstant * temperature_k));
}

double overburden_pressure(double density, double depth) {
    if (depth < 0.0) throw std::runtime_error("overburden_pressure: depth must be >= 0");
    return density * 9.80665 * depth;
}

double herron_langway_gradient(double coefficient, double density, double ice_density) {
    if (ice_density <= 0.0)
        throw std::runtime_error("herron_langway_gradient: ice density must be positive");
    return coefficient / ice_density * density * (ice_density - density);
}

double vertical_ice_velocity(double accumulation, double thickness_change, double basal_melt,
                             double z, double thickness, double exponent) {
    require_positive(thickness, "vertical_ice_velocity: thickness");
    if (z < 0.0 || z > thickness)
        throw std::runtime_error("vertical_ice_velocity: z outside the column");
    return -(accumulation - thickness_change - basal_melt) *
               std::pow(z / thickness, exponent + 1.0) -
           basal_melt;
}

double tipping_equilibrium(double thickness, double exponent, double gamma, double Gamma,
                           double offset) {
    if (thickness < 0.0)
        throw std::runtime_error("tipping_equilibrium: thickness must be >= 0");
    return std::pow(thickness, exponent) - gamma * Gamma * thickness - offset;
}

double tipping_critical_thickness(double exponent, double offset) {
    if (std::fabs(exponent - 1.0) < 1e-12 || exponent <= 0.0)
        throw std::runtime_error("tipping_critical_thickness: exponent must differ from 1");
    const double base = offset / (1.0 - exponent);
    if (base < 0.0)
        throw std::runtime_error("tipping_critical_thickness: no positive fold thickness");
    return std::pow(base, 1.0 / exponent);
}

double scalar_kalman_update(double background, double observation,
                            double background_variance, double observation_variance) {
    const double total = background_variance + observation_variance;
    require_positive(total, "scalar_kalman_update: variances");
    const double gain = background_variance / total;
    return background + gain * (observation - background);
}

}  // namespace climate
