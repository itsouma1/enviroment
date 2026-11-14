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

double idemix_mu0(double mu, double buoyancy_frequency, double coriolis, double n0, double f0) {
    if (buoyancy_frequency <= 0.0 || coriolis <= 0.0)
        throw std::runtime_error("idemix_mu0: N and f must be positive");
    if (n0 <= 0.0 || f0 <= 0.0)
        throw std::runtime_error("idemix_mu0: reference N0 and f0 must be positive");
    if (buoyancy_frequency <= coriolis)
        throw std::runtime_error("idemix_mu0: N must exceed f");
    return mu / std::acosh(n0 / f0);
}

double idemix_dissipation(double mu0, double coriolis, double buoyancy_frequency,
                          double gm_bandwidth, double wave_energy) {
    if (buoyancy_frequency <= coriolis)
        throw std::runtime_error("idemix_dissipation: N must exceed f");
    if (buoyancy_frequency <= 0.0)
        throw std::runtime_error("idemix_dissipation: N must be positive");
    const double f_e = coriolis * std::acosh(buoyancy_frequency / coriolis);
    const double ratio = gm_bandwidth / buoyancy_frequency;
    return mu0 * f_e * ratio * ratio * wave_energy * wave_energy;
}

double idemix_dissipation_derivative(double mu0, double coriolis, double buoyancy_frequency,
                                     double gm_bandwidth, double wave_energy) {
    if (buoyancy_frequency <= coriolis || buoyancy_frequency <= 0.0)
        throw std::runtime_error("idemix_dissipation_derivative: N must exceed positive f");
    const double f_e = coriolis * std::acosh(buoyancy_frequency / coriolis);
    const double ratio = gm_bandwidth / buoyancy_frequency;
    return 2.0 * mu0 * f_e * ratio * ratio * wave_energy;
}

double idemix_vertical_term(double vertical_group_velocity, double relaxation_time,
                            double second_derivative) {
    if (relaxation_time < 0.0)
        throw std::runtime_error("idemix_vertical_term: relaxation time must be >= 0");
    return -vertical_group_velocity * vertical_group_velocity * relaxation_time *
           second_derivative;
}

double kpp_bulk_richardson(double reference_buoyancy, double buoyancy_at_depth, double depth,
                           double reference_velocity, double velocity_at_depth,
                           double turbulent_velocity_scale) {
    if (depth < 0.0) throw std::runtime_error("kpp_bulk_richardson: depth must be >= 0");
    const double shear = reference_velocity - velocity_at_depth;
    const double denominator = shear * shear + turbulent_velocity_scale * turbulent_velocity_scale;
    if (denominator == 0.0)
        throw std::runtime_error("kpp_bulk_richardson: zero shear and turbulent scale");
    return (reference_buoyancy - buoyancy_at_depth) * depth / denominator;
}

double kpp_shape_function(double sigma, double a1, double a2, double a3) {
    if (sigma < 0.0 || sigma > 1.0)
        throw std::runtime_error("kpp_shape_function: sigma outside [0, 1]");
    return a1 * sigma + a2 * sigma * sigma + a3 * sigma * sigma * sigma;
}

double kpp_velocity_scale(double w_lambda, double sigma, double a1, double a2, double a3) {
    return w_lambda * kpp_shape_function(sigma, a1, a2, a3);
}

double kpp_bulk_richardson_sensitivity(double reference_buoyancy, double buoyancy_at_depth,
                                       double depth, double reference_velocity,
                                       double velocity_at_depth,
                                       double turbulent_velocity_scale) {
    const double shear = reference_velocity - velocity_at_depth;
    const double denominator = shear * shear + turbulent_velocity_scale * turbulent_velocity_scale;
    if (denominator == 0.0)
        throw std::runtime_error("kpp_bulk_richardson_sensitivity: zero denominator");
    return -2.0 * turbulent_velocity_scale * (reference_buoyancy - buoyancy_at_depth) * depth /
           (denominator * denominator);
}

HorizontalVector stochastic_pressure_acceleration(double dp_dx, double dp_dy,
                                                  double reference_density) {
    if (reference_density <= 0.0)
        throw std::runtime_error("stochastic_pressure_acceleration: density must be positive");
    HorizontalVector acceleration;
    acceleration.x = -dp_dx / reference_density;
    acceleration.y = -dp_dy / reference_density;
    return acceleration;
}

double bbm_effective_viscosity(double undamaged_viscosity, double damage, double damage_exponent,
                               double concentration_parameter, double concentration) {
    if (undamaged_viscosity < 0.0)
        throw std::runtime_error("bbm_effective_viscosity: viscosity must be >= 0");
    if (damage < 0.0 || damage > 1.0)
        throw std::runtime_error("bbm_effective_viscosity: damage must be in [0,1]");
    if (concentration < 0.0 || concentration > 1.0)
        throw std::runtime_error("bbm_effective_viscosity: concentration must be in [0,1]");
    return undamaged_viscosity * std::pow(1.0 - damage, damage_exponent) *
           std::exp(-damage_exponent * concentration_parameter * (1.0 - concentration));
}

double bbm_viscosity_damage_derivative(double undamaged_viscosity, double damage,
                                       double damage_exponent, double concentration_parameter,
                                       double concentration) {
    if (damage <= 0.0 || damage >= 1.0)
        throw std::runtime_error("bbm_viscosity_damage_derivative: damage must be in (0,1)");
    const double common = undamaged_viscosity * std::pow(1.0 - damage, damage_exponent - 1.0) *
                          std::exp(-damage_exponent * concentration_parameter *
                                   (1.0 - concentration));
    return -damage_exponent * common;
}

double mohr_coulomb_yield(double normal_stress, double cohesion, double friction_angle_deg) {
    if (cohesion < 0.0) throw std::runtime_error("mohr_coulomb_yield: cohesion must be >= 0");
    return cohesion + normal_stress * std::tan(deg_to_rad(friction_angle_deg));
}

double bbm_elastic_stress(double damage, double stiffness, double elastic_strain) {
    if (damage < 0.0 || damage > 1.0)
        throw std::runtime_error("bbm_elastic_stress: damage must be in [0,1]");
    return (1.0 - damage) * stiffness * elastic_strain;
}

double bbm_friction_stress(double normal_stress, double max_pressure, double elastic_stress,
                           double viscous_stress) {
    if (max_pressure < 0.0)
        throw std::runtime_error("bbm_friction_stress: P_max must be >= 0");
    if (normal_stress > 0.0) return viscous_stress;
    if (normal_stress > -max_pressure) return elastic_stress;
    return viscous_stress;
}

double teos10_gibbs_water(double t_celsius, double pressure_pa,
                          const std::vector<double>& g_jk, double gu, double tu, double pu) {
    if (g_jk.size() != 8 * 7)
        throw std::runtime_error("teos10_gibbs_water: expected 56 (8 x 7) coefficients");
    if (tu <= 0.0 || pu <= 0.0)
        throw std::runtime_error("teos10_gibbs_water: scaling constants must be positive");
    const double t = t_celsius / tu;
    const double p = pressure_pa / pu;
    double total = 0.0;
    double t_power = 1.0;
    for (int j = 0; j <= 7; ++j) {
        double p_power = 1.0;
        for (int k = 0; k <= 6; ++k) {
            total += g_jk[static_cast<std::size_t>(j * 7 + k)] * t_power * p_power;
            p_power *= p;
        }
        t_power *= t;
    }
    return gu * total;
}

double teos10_gibbs_saline(double salinity, double t_celsius, double pressure_pa,
                           const std::vector<SalineTerm>& terms, double gu, double su, double tu,
                           double pu) {
    if (salinity < 0.0) throw std::runtime_error("teos10_gibbs_saline: salinity must be >= 0");
    if (su <= 0.0 || tu <= 0.0 || pu <= 0.0)
        throw std::runtime_error("teos10_gibbs_saline: scaling constants must be positive");
    const double s = salinity / su;
    const double t = t_celsius / tu;
    const double p = pressure_pa / pu;
    double total = 0.0;
    for (const auto& term : terms) {
        if (term.i < 1) throw std::runtime_error("teos10_gibbs_saline: i must be >= 1");
        if (term.j < 0 || term.k < 0)
            throw std::runtime_error("teos10_gibbs_saline: j and k must be >= 0");
        double salinity_part = 0.0;
        if (term.i == 1) {
            if (s <= 0.0)
                throw std::runtime_error("teos10_gibbs_saline: i = 1 term needs S > 0");
            salinity_part = 0.5 * s * std::log(s);
        } else {
            salinity_part = std::pow(s, term.i / 2.0);
        }
        total += term.coefficient * salinity_part * std::pow(t, term.j) * std::pow(p, term.k);
    }
    return gu * total;
}

const std::vector<double>& pure_water_gibbs_coefficients() {
    // Row-major (j, k); coefficients not listed in the table are zero.
    static const std::vector<double> coefficients = {
        // j = 0
        0.101342743139674e3, 0.100015695367145e6, -0.254457654203630e4,
        0.284517778446287e3, -0.333146754253611e2, 0.420263108803084e1,
        -0.546428511471039,
        // j = 1
        0.590578347909402e1, -0.270983805184062e3, 0.776153611613101e3,
        -0.196512550881220e3, 0.289796526294175e2, -0.213290083518327e1, 0.0,
        // j = 2
        -0.123577859330390e5, 0.145503645404680e4, -0.756558385769359e3,
        0.273479662323528e3, -0.555604063817218e2, 0.434420671917197e1, 0.0,
        // j = 3
        0.736741204151612e3, -0.672507783145070e3, 0.499360390819152e3,
        -0.239545330654412e3, 0.488012518593872e2, -0.166307106208905e1, 0.0,
        // j = 4
        -0.148185936433658e3, 0.397968445406972e3, -0.301815380621876e3,
        0.152196371733841e3, -0.263748377232802e2, 0.0, 0.0,
        // j = 5
        0.580259125842571e2, -0.194618310617595e3, 0.120520654902025e3,
        -0.552723052340152e2, 0.648190668077221e1, 0.0, 0.0,
        // j = 6
        -0.189843846514172e2, 0.635113936641785e2, -0.222897317140459e2,
        0.817060541818112e1, 0.0, 0.0, 0.0,
        // j = 7
        0.305081646487967e1, -0.963108119393062e1, 0.0, 0.0, 0.0, 0.0, 0.0,
    };
    return coefficients;
}

GibbsDerivatives teos10_pure_water_derivatives(double t_celsius, double pressure_pa, double gu,
                                               double tu, double pu) {
    if (tu <= 0.0 || pu <= 0.0)
        throw std::runtime_error("teos10_pure_water_derivatives: scales must be positive");
    const auto& coefficients = pure_water_gibbs_coefficients();
    const double t = t_celsius / tu;
    const double p = pressure_pa / pu;

    double g = 0.0, gt = 0.0, gp = 0.0, gtt = 0.0, gtp = 0.0, gpp = 0.0;
    double t_power = 1.0;  // t^j
    for (int j = 0; j <= 7; ++j) {
        double p_power = 1.0;  // p^k
        for (int k = 0; k <= 6; ++k) {
            const double c = coefficients[static_cast<std::size_t>(j * 7 + k)];
            if (c != 0.0) {
                g += c * t_power * p_power;
                if (j >= 1) gt += c * j * std::pow(t, j - 1) * p_power;
                if (k >= 1) gp += c * k * t_power * std::pow(p, k - 1);
                if (j >= 2) gtt += c * j * (j - 1) * std::pow(t, j - 2) * p_power;
                if (j >= 1 && k >= 1)
                    gtp += c * j * k * std::pow(t, j - 1) * std::pow(p, k - 1);
                if (k >= 2) gpp += c * k * (k - 1) * t_power * std::pow(p, k - 2);
            }
            p_power *= p;
        }
        t_power *= t;
    }
    GibbsDerivatives derivatives;
    derivatives.g = gu * g;
    derivatives.g_t = gu / tu * gt;
    derivatives.g_p = gu / pu * gp;
    derivatives.g_tt = gu / (tu * tu) * gtt;
    derivatives.g_tp = gu / (tu * pu) * gtp;
    derivatives.g_pp = gu / (pu * pu) * gpp;
    return derivatives;
}

GibbsDerivatives teos10_saline_derivatives(double salinity, double t_celsius, double pressure_pa,
                                           const std::vector<SalineTerm>& terms, double gu,
                                           double su, double tu, double pu) {
    if (salinity < 0.0) throw std::runtime_error("teos10_saline_derivatives: salinity must be >= 0");
    if (su <= 0.0 || tu <= 0.0 || pu <= 0.0)
        throw std::runtime_error("teos10_saline_derivatives: scales must be positive");
    const double s = salinity / su;
    const double t = t_celsius / tu;
    const double p = pressure_pa / pu;

    double g = 0.0, gt = 0.0, gp = 0.0, gtt = 0.0, gtp = 0.0, gpp = 0.0;
    for (const auto& term : terms) {
        if (term.i < 1 || term.j < 0 || term.k < 0)
            throw std::runtime_error("teos10_saline_derivatives: invalid term indices");
        double salinity_part = 0.0;
        if (term.i == 1) {
            if (s <= 0.0)
                throw std::runtime_error("teos10_saline_derivatives: i = 1 term needs S > 0");
            salinity_part = 0.5 * s * std::log(s);
        } else {
            salinity_part = std::pow(s, term.i / 2.0);
        }
        const double tp = std::pow(t, term.j) * std::pow(p, term.k);
        g += term.coefficient * salinity_part * tp;
        if (term.j >= 1)
            gt += term.coefficient * salinity_part * term.j * std::pow(t, term.j - 1) *
                  std::pow(p, term.k);
        if (term.k >= 1)
            gp += term.coefficient * salinity_part * std::pow(t, term.j) * term.k *
                  std::pow(p, term.k - 1);
        if (term.j >= 2)
            gtt += term.coefficient * salinity_part * term.j * (term.j - 1) *
                   std::pow(t, term.j - 2) * std::pow(p, term.k);
        if (term.j >= 1 && term.k >= 1)
            gtp += term.coefficient * salinity_part * term.j * term.k *
                   std::pow(t, term.j - 1) * std::pow(p, term.k - 1);
        if (term.k >= 2)
            gpp += term.coefficient * salinity_part * std::pow(t, term.j) * term.k *
                   (term.k - 1) * std::pow(p, term.k - 2);
    }
    GibbsDerivatives derivatives;
    derivatives.g = gu * g;
    derivatives.g_t = gu / tu * gt;
    derivatives.g_p = gu / pu * gp;
    derivatives.g_tt = gu / (tu * tu) * gtt;
    derivatives.g_tp = gu / (tu * pu) * gtp;
    derivatives.g_pp = gu / (pu * pu) * gpp;
    return derivatives;
}

double sound_speed_from_derivatives(const GibbsDerivatives& derivatives) {
    const double denominator =
        derivatives.g_tp * derivatives.g_tp - derivatives.g_tt * derivatives.g_pp;
    if (denominator == 0.0 || derivatives.g_tt / denominator < 0.0)
        throw std::runtime_error("sound_speed_from_derivatives: non-physical curvature");
    return derivatives.g_p * std::sqrt(derivatives.g_tt / denominator);
}

double gibbs_sound_speed(const std::function<double(double, double, double)>& gibbs,
                         double salinity, double t_celsius, double pressure_pa,
                         double temperature_step, double pressure_step) {
    if (!gibbs) throw std::runtime_error("gibbs_sound_speed: empty Gibbs function");
    if (temperature_step <= 0.0 || pressure_step <= 0.0)
        throw std::runtime_error("gibbs_sound_speed: steps must be positive");
    const double h_t = temperature_step;
    const double h_p = pressure_step;
    const auto g = [&](double t, double p) { return gibbs(salinity, t, p); };

    const double g_t = (g(t_celsius + h_t, pressure_pa) - g(t_celsius - h_t, pressure_pa)) /
                       (2.0 * h_t);
    const double g_p = (g(t_celsius, pressure_pa + h_p) - g(t_celsius, pressure_pa - h_p)) /
                       (2.0 * h_p);
    const double g_tt = (g(t_celsius + h_t, pressure_pa) - 2.0 * g(t_celsius, pressure_pa) +
                         g(t_celsius - h_t, pressure_pa)) /
                        (h_t * h_t);
    const double g_pp = (g(t_celsius, pressure_pa + h_p) - 2.0 * g(t_celsius, pressure_pa) +
                         g(t_celsius, pressure_pa - h_p)) /
                        (h_p * h_p);
    const double g_tp = (g(t_celsius + h_t, pressure_pa + h_p) -
                         g(t_celsius + h_t, pressure_pa - h_p) -
                         g(t_celsius - h_t, pressure_pa + h_p) +
                         g(t_celsius - h_t, pressure_pa - h_p)) /
                        (4.0 * h_t * h_p);
    (void)g_t;
    const double denominator = g_tp * g_tp - g_tt * g_pp;
    if (denominator <= 0.0)
        throw std::runtime_error("gibbs_sound_speed: non-physical Gibbs curvature");
    return g_p * std::sqrt(g_tt / denominator);
}

}  // namespace climate
