#include "ocean.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "common.hpp"

namespace climate {
namespace {

double trapezoid(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2)
        throw std::runtime_error("trapezoid: need two or more paired samples");
    double total = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        const double dx = x[i + 1] - x[i];
        if (dx <= 0.0) throw std::runtime_error("trapezoid: x must increase");
        total += 0.5 * (y[i] + y[i + 1]) * dx;
    }
    return total;
}

}  // namespace

double buoyancy_frequency_squared(double gravity, double density, double density_gradient) {
    if (density <= 0.0) throw std::runtime_error("buoyancy_frequency_squared: density must be positive");
    return -gravity * density_gradient / density;
}

double mixed_layer_depth(const std::vector<double>& depth_m, const std::vector<double>& n2) {
    const double weight = trapezoid(depth_m, n2);
    if (weight == 0.0)
        throw std::runtime_error("mixed_layer_depth: integral of N^2 is zero");
    std::vector<double> moment(n2.size());
    for (std::size_t i = 0; i < n2.size(); ++i) moment[i] = depth_m[i] * n2[i];
    return trapezoid(depth_m, moment) / weight;
}

double ekman_velocity(double wind_stress, double coriolis, double ekman_depth) {
    if (coriolis == 0.0) throw std::runtime_error("ekman_velocity: Coriolis parameter is zero");
    if (ekman_depth <= 0.0) throw std::runtime_error("ekman_velocity: Ekman depth must be positive");
    return wind_stress / (coriolis * ekman_depth);
}

double barotropic_compensation(double ekman_velocity, double ekman_depth, double depth) {
    if (depth == 0.0) throw std::runtime_error("barotropic_compensation: depth must be non-zero");
    return -ekman_velocity * ekman_depth / depth;
}

double eddy_diffusivity_pv(double pv_flux, double beta, double slope_gradient) {
    const double denominator = -beta + slope_gradient;
    if (denominator == 0.0) throw std::runtime_error("eddy_diffusivity_pv: zero denominator");
    return pv_flux / denominator;
}

double eddy_diffusivity_gm(double buoyancy_flux, double velocity_gradient) {
    if (velocity_gradient == 0.0)
        throw std::runtime_error("eddy_diffusivity_gm: zero velocity gradient");
    return buoyancy_flux / velocity_gradient;
}

double internal_wave_frequency_squared(double f, double n_squared, double k, double m) {
    if (n_squared <= 0.0) throw std::runtime_error("internal_wave: N^2 must be positive");
    if (m == 0.0) throw std::runtime_error("internal_wave: vertical wavenumber must be non-zero");
    return f * f + n_squared * (k * k) / (m * m);
}

double gradient_richardson(double gravity, double reference_density, double density_gradient,
                           double shear_u, double shear_v) {
    const double shear_squared = shear_u * shear_u + shear_v * shear_v;
    if (shear_squared == 0.0)
        throw std::runtime_error("gradient_richardson: zero shear");
    if (reference_density <= 0.0)
        throw std::runtime_error("gradient_richardson: density must be positive");
    return -gravity * density_gradient / (reference_density * shear_squared);
}

double osborn_diffusivity(double mixing_efficiency, double dissipation_rate,
                          double n_squared) {
    if (n_squared <= 0.0) throw std::runtime_error("osborn_diffusivity: N^2 must be positive");
    return mixing_efficiency * dissipation_rate / n_squared;
}

double dissipation_from_shear(double kinematic_viscosity, double shear_variance,
                              double constant) {
    if (kinematic_viscosity < 0.0 || shear_variance < 0.0)
        throw std::runtime_error("dissipation_from_shear: inputs must be >= 0");
    return constant * kinematic_viscosity * shear_variance;
}

double nori_viscosity(double richardson_number, double shear_viscosity, double convective_viscosity,
                      double neutral_viscosity, double critical_richardson,
                      double transition_width) {
    if (critical_richardson <= 0.0 || transition_width <= 0.0)
        throw std::runtime_error("nori_viscosity: Ri_c and transition width must be positive");
    if (richardson_number < 0.0)
        return (shear_viscosity - convective_viscosity) *
                   std::tanh(richardson_number / transition_width) +
               shear_viscosity;
    if (richardson_number < critical_richardson) {
        return (neutral_viscosity - shear_viscosity) * richardson_number / critical_richardson +
               shear_viscosity;
    }
    return neutral_viscosity;
}

double martin_export_flux(double flux_at_100m, double depth_m, double exponent) {
    if (depth_m <= 0.0) throw std::runtime_error("martin_export_flux: depth must be positive");
    return flux_at_100m * std::pow(depth_m / 100.0, exponent);
}

double air_sea_co2_flux(double gas_transfer_velocity, double solubility, double pco2_sea,
                        double pco2_atmosphere) {
    return gas_transfer_velocity * solubility * (pco2_sea - pco2_atmosphere);
}

double dissolved_inorganic_carbon(double co2, double bicarbonate, double carbonate) {
    if (co2 < 0.0 || bicarbonate < 0.0 || carbonate < 0.0)
        throw std::runtime_error("dissolved_inorganic_carbon: concentrations must be >= 0");
    return co2 + bicarbonate + carbonate;
}

double latent_heat_flux_bulk(double air_density, double latent_heat_vaporization,
                             double transfer_coefficient, double wind_speed_10m,
                             double saturation_specific_humidity,
                             double actual_specific_humidity) {
    return air_density * latent_heat_vaporization * transfer_coefficient * wind_speed_10m *
           (saturation_specific_humidity - actual_specific_humidity);
}

double sensible_heat_flux_bulk(double air_density, double specific_heat,
                               double transfer_coefficient, double wind_speed_10m,
                               double surface_temperature, double air_temperature) {
    return air_density * specific_heat * transfer_coefficient * wind_speed_10m *
           (surface_temperature - air_temperature);
}

double net_air_sea_heat_flux(double shortwave, double longwave, double latent, double sensible) {
    return shortwave + longwave + latent + sensible;
}

double iron_limitation(double nitrogen, double kn, double phosphorus, double kp, double par,
                       double kpar, double iron, double kfe) {
    const double terms[4] = {nitrogen / (kn + nitrogen), phosphorus / (kp + phosphorus),
                             par / (kpar + par), iron / (kfe + iron)};
    double limit = terms[0];
    for (int i = 1; i < 4; ++i) limit = std::min(limit, terms[i]);
    return limit;
}

double free_iron_ligand_speciation(double total_ligand, double total_iron,
                                   double equilibrium_constant) {
    if (total_ligand < 0.0 || total_iron < 0.0)
        throw std::runtime_error("free_iron_ligand_speciation: totals must be >= 0");
    if (equilibrium_constant <= 0.0)
        throw std::runtime_error("free_iron_ligand_speciation: K_eq must be positive");
    if (total_iron == 0.0) return 0.0;

    // FeL = K (L_T - FeL) Fe', with FeL = Fe_T - Fe'.
    // K Fe'^2 + (K (L_T - Fe_T) + 1) Fe' - Fe_T = 0
    const double a = equilibrium_constant;
    const double b = equilibrium_constant * (total_ligand - total_iron) + 1.0;
    const double c = -total_iron;
    const double discriminant = b * b - 4.0 * a * c;
    return (-b + std::sqrt(discriminant)) / (2.0 * a);
}

double mcdr_potential_ratio(double dic, double alkalinity, double gamma_alkalinity,
                            double gamma_dic) {
    if (alkalinity == 0.0) throw std::runtime_error("mcdr_potential_ratio: alkalinity is zero");
    if (gamma_dic == 0.0) throw std::runtime_error("mcdr_potential_ratio: gamma_DIC is zero");
    return -(dic / alkalinity) * (gamma_alkalinity / gamma_dic);
}

double pierson_moskowitz_spectrum(double angular_frequency, double wind_speed_19_5m,
                                  double alpha, double beta, double gravity) {
    if (angular_frequency <= 0.0)
        throw std::runtime_error("pierson_moskowitz_spectrum: frequency must be positive");
    if (wind_speed_19_5m <= 0.0)
        throw std::runtime_error("pierson_moskowitz_spectrum: wind speed must be positive");
    const double omega0 = gravity / wind_speed_19_5m;
    return alpha * gravity * gravity / std::pow(angular_frequency, 5.0) *
           std::exp(-beta * std::pow(omega0 / angular_frequency, 4.0));
}

double doppler_shifted_frequency(double intrinsic_frequency, double wavenumber, double current) {
    return intrinsic_frequency + wavenumber * current;
}

double hibler_ice_pressure(double p_star, double thickness_m, double concentration, double c) {
    if (thickness_m < 0.0) throw std::runtime_error("hibler_ice_pressure: thickness must be >= 0");
    if (concentration < 0.0 || concentration > 1.0)
        throw std::runtime_error("hibler_ice_pressure: concentration must be in [0,1]");
    return p_star * thickness_m * std::exp(-c * (1.0 - concentration));
}

}  // namespace climate
