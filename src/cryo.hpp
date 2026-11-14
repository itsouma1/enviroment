#ifndef CLIMATE_CRYO_HPP
#define CLIMATE_CRYO_HPP

namespace climate {

// Horizontal vector used by the sea-ice momentum balance.
struct HorizontalStress {
    double x = 0.0;
    double y = 0.0;
};

// Glen's flow law -----------------------------------------------------------

// A(T) = A0 exp(-Q / (R T)), T in kelvin.
double glen_rate_factor(double reference_rate, double activation_energy,
                        double temperature_k);

// eps = A tau^n for the scalar effective quantities (n = 3 conventionally).
double glen_strain_rate(double rate_factor, double effective_stress, double n = 3.0);

// Effective strain from the second invariant: sqrt(0.5 eps_ij eps_ij).
double effective_strain(double e_xx, double e_yy, double e_zz, double e_xy, double e_xz,
                        double e_yz);

// Effective deviatoric stress: sqrt(0.5 tau_ij tau_ij).
double effective_stress(double t_xx, double t_yy, double t_zz, double t_xy, double t_xz,
                        double t_yz);

// mu = B / (2 eps_e^(n-1)).
double ice_viscosity(double hardness, double effective_strain, double n = 3.0);

// Shallow-ice velocity u = C (ds/dx) h^(n+1).
double sia_velocity(double coefficient, double surface_slope, double thickness,
                    double n = 3.0);

// Weertman-Budd sliding v_b = C_b tau_b^p / N_b^q.
double weertman_sliding(double coefficient, double basal_shear, double effective_pressure,
                        double p = 3.0, double q = 1.0);

// Coulomb-plastic sliding tau_b = mu N.
double coulomb_sliding(double friction, double effective_pressure);

// Sea-ice momentum balance ------------------------------------------------

// Bulk stress rho C |U| R U with a rotation matrix through `rotation` radians.
HorizontalStress bulk_ice_stress(double density, double drag, double fluid_u,
                                 double fluid_v, double ice_u, double ice_v,
                                 double rotation);

// du/dt and dv/dt from m Du/Dt = -m f k x u + tau_air + tau_ocean
// - m grad(phi) + F, with F the internal-stress divergence.
void sea_ice_momentum_tendency(double& du_dt, double& dv_dt, double u, double v,
                               double mass, double coriolis,
                               const HorizontalStress& air,
                               const HorizontalStress& ocean,
                               double potential_gradient_x, double potential_gradient_y,
                               const HorizontalStress& internal);

// Snow and surface energy balance -----------------------------------------

// Melt energy (1 - alpha) Qs + QL - eps sigma Ts^4.
double surface_melt_energy(double albedo, double shortwave, double longwave,
                           double emissivity, double surface_temperature);

// Permafrost ---------------------------------------------------------------

// Stefan condition: rho L dz_f/dt = k_f dT_f/dz - k_u dT_u/dz.
double stefan_front_velocity(double frozen_conductivity, double frozen_gradient,
                             double thawed_conductivity, double thawed_gradient,
                             double density, double latent_heat);

// ST = c_b rho_b (Ts - Tf) / (S_w f rho_w phi L_f).
double stefan_number(double bulk_heat_capacity, double bulk_density,
                     double surface_temperature, double freezing_temperature,
                     double water_saturation, double freeze_fraction, double water_density,
                     double porosity, double latent_heat);

// Z = sqrt(2 s lambda TDD / L).
double active_layer_depth(double scaling, double conductivity, double thawing_degree_days,
                          double latent_heat);

// Ice-ocean interaction ---------------------------------------------------

// Three-equation basal melt (metres of ice per second).
double basal_melt_three_equation(double ocean_density, double ice_density, double drag,
                                 double transfer_coefficient, double latent_heat,
                                 double ice_heat_capacity, double ice_temperature,
                                 double freezing_temperature, double ocean_heat_capacity,
                                 double current_speed, double boundary_temperature);

// m = (c_oc / L_i) (rho_oc / rho_i) gamma_Ts eps (T_loc - T_t,loc).
double basal_melt_simplified(double ocean_heat_capacity, double latent_heat,
                             double ocean_density, double ice_density,
                             double transfer_coefficient, double mixing_efficiency,
                             double local_temperature, double local_freezing);

// Depth-dependent plume melt.
double plume_melt(double ambient_temperature, double grounding_line_freezing,
                  double melt_function, double haline_contraction, double ambient_salinity,
                  double gravity, double exchange_velocity, double lambda3, double plume_length,
                  double density_ratio, double gamma, double drag, double entrainment,
                  double sin_slope, double tau_coefficient);

// Subglacial hydrology ----------------------------------------------------

// q = -k_q b^alpha |grad phi|^(beta-1) grad phi.
double subglacial_flux(double coefficient, double gap_height, double alpha, double beta,
                       double head_gradient);

// phi = rho_w g b + p_w.
double hydraulic_potential(double water_density, double gravity, double bed_elevation,
                           double water_pressure);

// Bed melt including geothermal, frictional and dissipation terms.
double subglacial_melt(double geothermal_flux, double basal_shear, double basal_speed,
                       double discharge, double head_gradient, double water_heat_capacity,
                       double water_density, double pressure_gradient, double ice_density,
                       double latent_heat);

// Calving -----------------------------------------------------------------

double calving_flux(double calving_speed, double cross_section);

// L_r = 0.33 sqrt(3.75 H + S_t^2).
double phase_field_notch(double thickness, double tensile_strength);

// Firn --------------------------------------------------------------------

// Arthern et al. (2010) dry-firn densification rate; k_c switches at 550 kg/m^3.
double firn_densification_rate(double overburden, double grain_radius, double ice_density,
                               double firn_density, double activation_energy,
                               double temperature_k);

// sigma = g int_z^zs rho dz, for a uniform density column.
double overburden_pressure(double density, double depth);

// d rho / dz = (c / rho_ice) rho (rho_ice - rho).
double herron_langway_gradient(double coefficient, double density, double ice_density);

// Ice-core age and flow ---------------------------------------------------

double vertical_ice_velocity(double accumulation, double thickness_change, double basal_melt,
                             double z, double thickness, double exponent);

// Hysteresis --------------------------------------------------------------

double tipping_equilibrium(double thickness, double exponent, double gamma, double Gamma,
                           double offset);

// Thickness at the fold where F(h) = F'(h) = 0.
double tipping_critical_thickness(double exponent, double offset);

// Scalar Kalman update X_a = X_b + K (X_o - X_b).
double scalar_kalman_update(double background, double observation,
                            double background_variance, double observation_variance);

// Advanced cryosphere (notes 112-128 KB) ----------------------------------

// eta_tilde = (1 - D) eta, D in [0, 1].
double damaged_viscosity(double intact_viscosity, double damage);

// Glen law with the (1 - D)^-n strain enhancement.
double damaged_glen_strain(double rate_factor, double damage, double effective_stress,
                           double n = 3.0);

// Biot poroelastic Cauchy stress: (1 - D) sigma_bar - D P_h.
double biot_effective_stress(double intact_stress, double damage, double pore_pressure);

// Maxwell elasto-brittle stress (1 - D) C : eps, isotropic stiffness.
double maxwell_elasto_brittle_stress(double damage, double stiffness, double elastic_strain);

// Lagrangian floe acceleration: (sum F + tau_air + tau_water) / m.
HorizontalStress floe_acceleration(const HorizontalStress& contact_force,
                                   const HorizontalStress& air_stress,
                                   const HorizontalStress& water_stress, double mass);

// Rothlisberger channel: dS/dt = (rho_w/rho_i) M - A_tilde S N^n.
double roethlisberger_channel_rate(double melt_input, double channel_area,
                                   double effective_pressure, double n,
                                   double closure_coefficient, double water_density,
                                   double ice_density);

// Channel melt from rho_w M L + rho_w c_w S dT/dt + rho_w c_w Q dT/ds = -Q dphi/ds.
double channel_energy_melt(double discharge, double head_gradient, double channel_area,
                           double temperature_change, double temperature_gradient,
                           double water_density, double water_heat_capacity,
                           double latent_heat);

// Phase-field free-energy density (tensile + compressive + fracture terms).
double phase_field_energy_density(double phi, double gradient_squared,
                                  double tensile_energy, double compressive_energy,
                                  double fracture_toughness, double length_scale);

// dphi/dt = -M delta F / delta phi, with delta F/delta phi = -2(1-phi) W+ 
// + (Gc/ell)(phi - ell^2 laplacian(phi)).
double phase_field_evolution(double phi, double laplacian, double tensile_energy,
                             double fracture_toughness, double length_scale,
                             double mobility);

// d(rho Q U)/dxi = rho m_dot.
double plume_flux_derivative(double density, double melt_rate);

// d(rho Q U^2)/dxi = delta_rho g h sin(theta) - tau.
double plume_momentum_derivative(double density_difference, double gravity, double thickness,
                                 double sin_slope, double drag_stress);

// Probabilistic age increment D(z) / (A(z) tau(z)) dz.
double probabilistic_age_increment(double correction, double accumulation, double thinning,
                                   double dz);

}  // namespace climate

#endif
