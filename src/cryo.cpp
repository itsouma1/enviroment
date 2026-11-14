#include "cryo.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

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

double damaged_viscosity(double intact_viscosity, double damage) {
    require_positive(intact_viscosity, "damaged_viscosity: viscosity");
    if (damage < 0.0 || damage >= 1.0)
        throw std::runtime_error("damaged_viscosity: damage must lie in [0, 1)");
    return (1.0 - damage) * intact_viscosity;
}

double damaged_glen_strain(double rate_factor, double damage, double effective_stress,
                           double n) {
    require_positive(rate_factor, "damaged_glen_strain: rate factor");
    if (damage < 0.0 || damage >= 1.0)
        throw std::runtime_error("damaged_glen_strain: damage must lie in [0, 1)");
    if (effective_stress < 0.0)
        throw std::runtime_error("damaged_glen_strain: stress must be >= 0");
    return rate_factor / std::pow(1.0 - damage, n) * std::pow(effective_stress, n);
}

double biot_effective_stress(double intact_stress, double damage, double pore_pressure) {
    if (damage < 0.0 || damage > 1.0)
        throw std::runtime_error("biot_effective_stress: damage outside [0, 1]");
    return (1.0 - damage) * intact_stress - damage * pore_pressure;
}

double maxwell_elasto_brittle_stress(double damage, double stiffness, double elastic_strain) {
    if (damage < 0.0 || damage > 1.0)
        throw std::runtime_error("maxwell_elasto_brittle_stress: damage outside [0, 1]");
    return (1.0 - damage) * stiffness * elastic_strain;
}

HorizontalStress floe_acceleration(const HorizontalStress& contact_force,
                                   const HorizontalStress& air_stress,
                                   const HorizontalStress& water_stress, double mass) {
    require_positive(mass, "floe_acceleration: mass");
    HorizontalStress acceleration;
    acceleration.x = (contact_force.x + air_stress.x + water_stress.x) / mass;
    acceleration.y = (contact_force.y + air_stress.y + water_stress.y) / mass;
    return acceleration;
}

double roethlisberger_channel_rate(double melt_input, double channel_area,
                                   double effective_pressure, double n,
                                   double closure_coefficient, double water_density,
                                   double ice_density) {
    require_positive(ice_density, "roethlisberger_channel_rate: ice density");
    if (channel_area < 0.0)
        throw std::runtime_error("roethlisberger_channel_rate: area must be >= 0");
    if (effective_pressure < 0.0)
        throw std::runtime_error("roethlisberger_channel_rate: pressure must be >= 0");
    return water_density / ice_density * melt_input -
           closure_coefficient * channel_area * std::pow(effective_pressure, n);
}

double channel_energy_melt(double discharge, double head_gradient, double channel_area,
                           double temperature_change, double temperature_gradient,
                           double water_density, double water_heat_capacity,
                           double latent_heat) {
    require_positive(water_density, "channel_energy_melt: density");
    require_positive(latent_heat, "channel_energy_melt: latent heat");
    if (channel_area < 0.0)
        throw std::runtime_error("channel_energy_melt: area must be >= 0");
    const double energy = -discharge * head_gradient -
                          water_density * water_heat_capacity * channel_area *
                              temperature_change -
                          water_density * water_heat_capacity * discharge *
                              temperature_gradient;
    return energy / (water_density * latent_heat);
}

double phase_field_energy_density(double phi, double gradient_squared,
                                  double tensile_energy, double compressive_energy,
                                  double fracture_toughness, double length_scale) {
    require_positive(length_scale, "phase_field_energy_density: length scale");
    if (phi < 0.0 || phi > 1.0)
        throw std::runtime_error("phase_field_energy_density: phi outside [0, 1]");
    return (1.0 - phi) * (1.0 - phi) * tensile_energy + compressive_energy +
           fracture_toughness / (2.0 * length_scale) *
               (phi * phi + length_scale * length_scale * gradient_squared);
}

double phase_field_evolution(double phi, double laplacian, double tensile_energy,
                             double fracture_toughness, double length_scale,
                             double mobility) {
    require_positive(length_scale, "phase_field_evolution: length scale");
    if (phi < 0.0 || phi > 1.0)
        throw std::runtime_error("phase_field_evolution: phi outside [0, 1]");
    const double variational =
        -2.0 * (1.0 - phi) * tensile_energy +
        fracture_toughness / length_scale *
            (phi - length_scale * length_scale * laplacian);
    return -mobility * variational;
}

double plume_flux_derivative(double density, double melt_rate) {
    return density * melt_rate;
}

double plume_momentum_derivative(double density_difference, double gravity, double thickness,
                                 double sin_slope, double drag_stress) {
    if (thickness < 0.0)
        throw std::runtime_error("plume_momentum_derivative: thickness must be >= 0");
    return density_difference * gravity * thickness * sin_slope - drag_stress;
}

double probabilistic_age_increment(double correction, double accumulation, double thinning,
                                   double dz) {
    require_positive(accumulation, "probabilistic_age_increment: accumulation");
    require_positive(thinning, "probabilistic_age_increment: thinning");
    return correction / (accumulation * thinning) * dz;
}

double hayhurst_stress(double alpha, double beta, double max_principal,
                       double von_mises, double trace) {
    if (alpha < 0.0 || beta < 0.0 || alpha + beta > 1.0)
        throw std::runtime_error("hayhurst_stress: weights outside the simplex");
    return alpha * max_principal + beta * von_mises + (1.0 - alpha - beta) * trace;
}

double hayhurst_principal_sensitivity(double alpha) {
    if (alpha < 0.0 || alpha > 1.0)
        throw std::runtime_error("hayhurst_principal_sensitivity: alpha outside [0, 1]");
    return alpha;
}

double hayhurst_von_mises_sensitivity(double beta) {
    if (beta < 0.0 || beta > 1.0)
        throw std::runtime_error("hayhurst_von_mises_sensitivity: beta outside [0, 1]");
    return beta;
}

double hayhurst_damage_rate(double hayhurst, double threshold, double exponent) {
    require_positive(threshold, "hayhurst_damage_rate: threshold");
    require_positive(exponent, "hayhurst_damage_rate: exponent");
    if (hayhurst <= threshold) return 0.0;
    return std::pow(hayhurst / threshold, exponent);
}

double hayhurst_damage_rate_derivative(double hayhurst, double threshold, double exponent) {
    require_positive(threshold, "hayhurst_damage_rate_derivative: threshold");
    require_positive(exponent, "hayhurst_damage_rate_derivative: exponent");
    if (hayhurst <= threshold) return 0.0;
    return exponent / threshold * std::pow(hayhurst / threshold, exponent - 1.0);
}

ContactForce hertz_mindlin_contact(double youngs_modulus, double poisson_ratio,
                                   double damping_ratio, double packing_parameter,
                                   double contact_radius, double centroid_distance,
                                   double mass_i, double mass_j, double normal_displacement,
                                   double normal_velocity, double tangential_displacement,
                                   double tangential_velocity, double friction) {
    require_positive(youngs_modulus, "hertz_mindlin_contact: Young's modulus");
    require_positive(centroid_distance, "hertz_mindlin_contact: centroid distance");
    require_positive(contact_radius, "hertz_mindlin_contact: contact radius");
    require_positive(mass_i + mass_j, "hertz_mindlin_contact: masses");
    if (poisson_ratio <= -1.0 || poisson_ratio >= 0.5)
        throw std::runtime_error("hertz_mindlin_contact: Poisson ratio outside (-1, 0.5)");
    const double area = 3.14159265358979323846 * contact_radius * contact_radius;
    const double equivalent_mass = mass_i * mass_j / (mass_i + mass_j);
    const double normal_stiffness =
        packing_parameter * area / centroid_distance * youngs_modulus;
    // C_n = 2 xi sqrt(m_eq K_n); the geometric factor belongs to K_n, not to
    // the damping coefficient (otherwise C_n would carry units of area).
    const double normal_damping =
        2.0 * damping_ratio * std::sqrt(equivalent_mass * normal_stiffness);
    const double shear_stiffness = normal_stiffness / (2.0 * (1.0 + poisson_ratio));
    const double shear_damping = 2.0 * damping_ratio * std::sqrt(equivalent_mass * shear_stiffness);

    ContactForce force;
    force.normal = normal_stiffness * normal_displacement + normal_damping * normal_velocity;
    double tangential = shear_stiffness * tangential_displacement + shear_damping * tangential_velocity;
    const double limit = friction * std::fabs(force.normal);
    if (std::fabs(tangential) > limit) tangential = std::copysign(limit, tangential);
    force.tangential = tangential;
    force.normal_stiffness = normal_stiffness;
    force.shear_stiffness = shear_stiffness;
    return force;
}

double firn_strain_softening(double e_xx, double e_yy, double e_zz, double e_xy, double e_xz,
                             double e_yz, double classical_vertical_strain) {
    if (!(std::fabs(classical_vertical_strain) > 0.0))
        throw std::runtime_error("firn_strain_softening: classical strain must be non-zero");
    const double norm = std::sqrt(e_xx * e_xx + e_yy * e_yy + e_zz * e_zz + 2.0 * e_xy * e_xy +
                                  2.0 * e_xz * e_xz + 2.0 * e_yz * e_yz);
    return norm / std::fabs(classical_vertical_strain);
}

double softened_firn_rate(double classical_rate, double softening) {
    if (softening < 0.0)
        throw std::runtime_error("softened_firn_rate: softening must be >= 0");
    return classical_rate * softening;
}

namespace {
std::vector<std::vector<double>> cholesky_lower(const std::vector<std::vector<double>>& matrix) {
    const std::size_t n = matrix.size();
    std::vector<std::vector<double>> lower(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        if (matrix[i].size() != n)
            throw std::runtime_error("icechrono: covariance must be square");
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k) sum -= lower[i][k] * lower[j][k];
            if (i == j) {
                if (!(sum > 0.0))
                    throw std::runtime_error("icechrono: covariance is not positive definite");
                lower[i][j] = std::sqrt(sum);
            } else {
                lower[i][j] = sum / lower[j][j];
            }
        }
    }
    return lower;
}

std::vector<double> forward_substitute(const std::vector<std::vector<double>>& lower,
                                       const std::vector<double>& rhs) {
    const std::size_t n = lower.size();
    std::vector<double> solution(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = rhs[i];
        for (std::size_t k = 0; k < i; ++k) sum -= lower[i][k] * solution[k];
        solution[i] = sum / lower[i][i];
    }
    return solution;
}

std::vector<double> back_substitute(const std::vector<std::vector<double>>& lower,
                                    const std::vector<double>& rhs) {
    const std::size_t n = lower.size();
    std::vector<double> solution(n, 0.0);
    for (std::size_t offset = 0; offset < n; ++offset) {
        const std::size_t i = n - 1 - offset;
        double sum = rhs[i];
        for (std::size_t k = i + 1; k < n; ++k) sum -= lower[k][i] * solution[k];
        solution[i] = sum / lower[i][i];
    }
    return solution;
}
}  // namespace

double icechrono_cost(const std::vector<std::vector<double>>& covariance,
                      const std::vector<double>& residuals) {
    if (covariance.size() != residuals.size())
        throw std::runtime_error("icechrono_cost: size mismatch");
    const auto lower = cholesky_lower(covariance);
    const auto whitened = forward_substitute(lower, residuals);
    double cost = 0.0;
    for (double value : whitened) cost += value * value;
    return cost;
}

std::vector<double> icechrono_gradient(const std::vector<std::vector<double>>& covariance,
                                       const std::vector<double>& residuals,
                                       const std::vector<std::vector<double>>& jacobian) {
    if (jacobian.size() != residuals.size())
        throw std::runtime_error("icechrono_gradient: Jacobian row mismatch");
    const auto lower = cholesky_lower(covariance);
    const auto whitened = forward_substitute(lower, residuals);
    const auto weighted = back_substitute(lower, whitened);  // P^-1 r
    const std::size_t parameters = jacobian.empty() ? 0 : jacobian[0].size();
    std::vector<double> gradient(parameters, 0.0);
    for (std::size_t i = 0; i < jacobian.size(); ++i) {
        if (jacobian[i].size() != parameters)
            throw std::runtime_error("icechrono_gradient: ragged Jacobian");
        for (std::size_t j = 0; j < parameters; ++j)
            gradient[j] += 2.0 * jacobian[i][j] * weighted[i];
    }
    return gradient;
}

}  // namespace climate
