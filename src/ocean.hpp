#pragma once

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

}  // namespace climate
