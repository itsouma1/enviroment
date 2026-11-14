#pragma once

#include <functional>
#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Stratification, mixed layer and overturning
// ---------------------------------------------------------------------------
// N^2 = -(g/rho) d rho/dz
double buoyancy_frequency_squared(double gravity, double density, double density_gradient);

// Bulk mixed-layer depth as the N^2-weighted centroid, h = int z N^2 dz / int N^2 dz.
// The notes give h = int 2 N^2 dz / int N^2 dz, which is dimensionless and
// therefore cannot be a depth; the centroid is used instead.
double mixed_layer_depth(const std::vector<double>& depth_m, const std::vector<double>& n2);

// Ekman surface velocity v_e = tau / (f h) and its barotropic compensation.
double ekman_velocity(double wind_stress, double coriolis, double ekman_depth);
double barotropic_compensation(double ekman_velocity, double ekman_depth, double depth);

// Eddy diffusivities: K_R = vq / (-beta + S'(z)), K_GM = vb / U'(z).
double eddy_diffusivity_pv(double pv_flux, double beta, double slope_gradient);
double eddy_diffusivity_gm(double buoyancy_flux, double velocity_gradient);

// Internal-wave dispersion: omega^2 = f^2 + N^2 (k^2 / m^2).
// The notes carry an extra g^2/(rho^2 N^2) factor that is dimensionally
// inconsistent; the standard Boussinesq form is implemented here.
double internal_wave_frequency_squared(double coriolis, double buoyancy_frequency_squared,
                                       double horizontal_wavenumber, double vertical_wavenumber);

// Gradient Richardson number and the bulk Richardson criterion used by KPP.
double gradient_richardson(double gravity, double reference_density, double density_gradient,
                           double shear_u, double shear_v);

// ---------------------------------------------------------------------------
// Turbulent mixing
// ---------------------------------------------------------------------------
// Osborn: kappa = Gamma epsilon / N^2.
double osborn_diffusivity(double mixing_efficiency, double dissipation_rate,
                          double buoyancy_frequency_squared);

// epsilon = 7.5 nu (du'/dz)^2.
double dissipation_from_shear(double kinematic_viscosity, double shear_variance,
                              double constant = 7.5);

// NORi piecewise viscosity closure.
double nori_viscosity(double richardson_number, double shear_viscosity, double convective_viscosity,
                      double neutral_viscosity, double critical_richardson,
                      double transition_width);

// ---------------------------------------------------------------------------
// Carbon cycle and biogeochemistry
// ---------------------------------------------------------------------------
// Martin et al. (1987): F = F100 (z/100)^b.
double martin_export_flux(double flux_at_100m, double depth_m, double exponent);

// F = k K0 (pCO2_sea - pCO2_atm).
double air_sea_co2_flux(double gas_transfer_velocity, double solubility, double pco2_sea,
                        double pco2_atmosphere);

// DIC = [CO2] + [HCO3-] + [CO3^2-].
double dissolved_inorganic_carbon(double co2, double bicarbonate, double carbonate);

// Air-sea heat flux components.
double latent_heat_flux_bulk(double air_density, double latent_heat_vaporization,
                             double transfer_coefficient, double wind_speed_10m,
                             double saturation_specific_humidity,
                             double actual_specific_humidity);
double sensible_heat_flux_bulk(double air_density, double specific_heat,
                               double transfer_coefficient, double wind_speed_10m,
                               double surface_temperature, double air_temperature);
double net_air_sea_heat_flux(double shortwave, double longwave, double latent, double sensible);

// Iron and nutrient co-limitation: min of Michaelis-Menten terms.
double iron_limitation(double nitrogen, double kn, double phosphorus, double kp, double par,
                       double kpar, double iron, double kfe);

// Ligand speciation: Fe' free iron from total ligand, total iron and K_eq.
double free_iron_ligand_speciation(double total_ligand, double total_iron, double equilibrium_constant);

// mCDR potential: (DIC/Alk) (gamma_Alk / gamma_DIC).
double mcdr_potential_ratio(double dic, double alkalinity, double gamma_alkalinity,
                            double gamma_dic);

// ---------------------------------------------------------------------------
// Wind waves
// ---------------------------------------------------------------------------
// Pierson-Moskowitz: S(omega) = alpha g^2 / omega^5 exp[-beta (omega0/omega)^4],
// omega0 = g / U19.5.
double pierson_moskowitz_spectrum(double angular_frequency, double wind_speed_19_5m,
                                  double alpha = 8.1e-3, double beta = 0.74,
                                  double gravity = 9.80665);

// Doppler-shifted intrinsic frequency: omega = sigma + k U.
double doppler_shifted_frequency(double intrinsic_frequency, double wavenumber, double current);

// ---------------------------------------------------------------------------
// Sea ice rheology and momentum
// ---------------------------------------------------------------------------
// Hibler pressure: P = P* h exp(-C (1 - A)).
double hibler_ice_pressure(double p_star, double thickness_m, double concentration, double c);

// ---------------------------------------------------------------------------
// IDEMIX internal-wave energetics
// ---------------------------------------------------------------------------
// mu0 = mu / arccosh(N0/f0); f_e = f arccosh(N/f);
// eps = mu0 f_e (m_*^2 / N^2) E^2.
double idemix_mu0(double mu, double buoyancy_frequency, double coriolis, double n0, double f0);

double idemix_dissipation(double mu0, double coriolis, double buoyancy_frequency,
                          double gm_bandwidth, double wave_energy);

// d(eps)/dE = 2 mu0 f_e (m_*^2/N^2) E.
double idemix_dissipation_derivative(double mu0, double coriolis, double buoyancy_frequency,
                                     double gm_bandwidth, double wave_energy);

// Vertical redistribution term -d/dz(c0 tau_v d(c0 E)/dz) for constant c0.
double idemix_vertical_term(double vertical_group_velocity, double relaxation_time,
                            double second_derivative);

// ---------------------------------------------------------------------------
// KPP bulk Richardson number
// ---------------------------------------------------------------------------
// Rib = [B_r - B(d)] d / (|V_r - V(d)|^2 + V_t^2).
double kpp_bulk_richardson(double reference_buoyancy, double buoyancy_at_depth, double depth,
                           double reference_velocity, double velocity_at_depth,
                           double turbulent_velocity_scale);

// dRib/dVt = -2 Vt [B_r - B(d)] d / (|dV|^2 + Vt^2)^2.
double kpp_bulk_richardson_sensitivity(double reference_buoyancy, double buoyancy_at_depth,
                                       double depth, double reference_velocity,
                                       double velocity_at_depth,
                                       double turbulent_velocity_scale);

// ---------------------------------------------------------------------------
// Stochastic non-hydrostatic pressure
// ---------------------------------------------------------------------------
struct HorizontalVector {
    double x = 0.0;
    double y = 0.0;
};

// Acceleration from the martingale pressure: -(1/rho0) grad_H(dp^sigma).
// Under the strong hydrostatic hypothesis its vertical derivative vanishes.
HorizontalVector stochastic_pressure_acceleration(double dp_dx, double dp_dy,
                                                  double reference_density);

// ---------------------------------------------------------------------------
// Sea-ice Brittle Bingham-Maxwell rheology
// ---------------------------------------------------------------------------
// eta = eta0 (1-d)^alpha exp[-alpha c* (1-A)].
double bbm_effective_viscosity(double undamaged_viscosity, double damage, double damage_exponent,
                               double concentration_parameter, double concentration);

// d(eta)/dd = -alpha eta0 (1-d)^(alpha-1) exp[-alpha c* (1-A)].
double bbm_viscosity_damage_derivative(double undamaged_viscosity, double damage,
                                       double damage_exponent, double concentration_parameter,
                                       double concentration);

// Mohr-Coulomb envelope: tau = c + sigma_N tan(phi).
double mohr_coulomb_yield(double normal_stress, double cohesion, double friction_angle_deg);

// Damaged elasticity: sigma = (1 - D) C epsilon.
double bbm_elastic_stress(double damage, double stiffness, double elastic_strain);

// Friction-element regimes: divergent (sigma_N > 0), elastic, and convergent
// beyond -P_max (viscous flow).
double bbm_friction_stress(double normal_stress, double max_pressure, double elastic_stress,
                           double viscous_stress);

// ---------------------------------------------------------------------------
// TEOS-10 Gibbs function (generic evaluators; coefficient tables are external)
// ---------------------------------------------------------------------------
struct SalineTerm {
    int i = 1;  // salinity power i/2, i >= 1
    int j = 0;  // temperature power
    int k = 0;  // pressure power
    double coefficient = 0.0;
};

// g^W = gu sum_jk g_jk (t/tu)^j (p/pu)^k with g_jk stored row-major, 8 x 7.
double teos10_gibbs_water(double t_celsius, double pressure_pa,
                          const std::vector<double>& g_jk, double gu = 1.0, double tu = 40.0,
                          double pu = 1e8);

// g^S = gu sum { g1jk/2 (S/Su) ln(S/Su) + sum_{i>1} gijk (S/Su)^(i/2) } (t/tu)^j (p/pu)^k.
double teos10_gibbs_saline(double salinity, double t_celsius, double pressure_pa,
                           const std::vector<SalineTerm>& terms, double gu = 1.0,
                           double su = 40.0, double tu = 40.0, double pu = 1e8);

struct GibbsDerivatives {
    double g = 0.0;    // J/kg
    double g_t = 0.0;  // J/kg/K
    double g_p = 0.0;  // m^3/kg
    double g_tt = 0.0;
    double g_tp = 0.0;
    double g_pp = 0.0;
};

// The 41-term IAPWS-2009 pure-water table, stored row-major as 8 x 7 (j, k).
const std::vector<double>& pure_water_gibbs_coefficients();

GibbsDerivatives teos10_pure_water_derivatives(double t_celsius, double pressure_pa,
                                               double gu = 1.0, double tu = 40.0,
                                               double pu = 1e8);

GibbsDerivatives teos10_saline_derivatives(double salinity, double t_celsius,
                                           double pressure_pa,
                                           const std::vector<SalineTerm>& terms, double gu = 1.0,
                                           double su = 40.0, double tu = 40.0,
                                           double pu = 1e8);

// c = g_p sqrt(g_tt / (g_tp^2 - g_tt g_pp)) from a derivative tuple.
double sound_speed_from_derivatives(const GibbsDerivatives& derivatives);

// c = g_P sqrt(g_TT / (g_TP^2 - g_TT g_PP)), with central differences of g(S,t,p).
double gibbs_sound_speed(const std::function<double(double, double, double)>& gibbs,
                         double salinity, double t_celsius, double pressure_pa,
                         double temperature_step = 1e-4, double pressure_step = 1e2);

}  // namespace climate
