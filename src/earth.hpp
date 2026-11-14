#pragma once

#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Advanced climate sensitivity
// ---------------------------------------------------------------------------
// Planck (blackbody) feedback: lambda_0 = dR/dT = -4 sigma T^3.
double planck_feedback(double surface_temperature_k);

// Pattern-corrected ECS: ECS = -F_2x / (alpha + alpha').
double pattern_corrected_ecs(double historical_feedback, double pattern_correction,
                             double forcing_2xco2 = 3.7);

// Two-layer energy balance model (Gregory 2000; Geoffroy et al. 2013):
//   C   dT/dt  = F - lambda T - eps gamma (T - Td)
//   Cd  dTd/dt = gamma (T - Td)
struct TwoLayerModel {
    double heat_capacity = 8.0;        // W yr m^-2 K^-1 (mixed layer + atmosphere)
    double deep_heat_capacity = 100.0; // W yr m^-2 K^-1
    double feedback = -1.25;           // W m^-2 K^-1
    double uptake_efficiency = 0.7;    // gamma, W m^-2 K^-1
    double efficacy = 1.0;             // epsilon
};

struct TwoLayerState {
    double surface = 0.0;  // K
    double deep = 0.0;     // K
};

// Equilibrium warming of the two-layer model for a step forcing: F / -lambda.
double two_layer_equilibrium(const TwoLayerModel& model, double forcing);

// Integrate the two-layer model with classical RK4 under a step forcing.
TwoLayerState integrate_two_layer(const TwoLayerModel& model, double forcing, double years,
                                  double dt_years = 0.05,
                                  TwoLayerState initial = TwoLayerState{});

// Ocean heat uptake gamma (T - Td), and the uptake efficiency OHU / dT.
double ocean_heat_uptake(const TwoLayerModel& model, const TwoLayerState& state);
double ocean_heat_uptake_efficiency(const TwoLayerModel& model, const TwoLayerState& state);

// Transient climate response estimate for a doubling: F_2x / -(lambda + eps gamma).
double transient_climate_response(const TwoLayerModel& model, double forcing_2xco2 = 3.7);

// Equilibrium-to-transient ratio ECS / TCR.
double equilibration_ratio(const TwoLayerModel& model, double forcing_2xco2 = 3.7);

// ---------------------------------------------------------------------------
// Carbon-cycle feedbacks
// ---------------------------------------------------------------------------
struct CarbonFeedback {
    double beta_land = 0.0;   // dC_land / dC_atmos
    double beta_ocean = 0.0;  // dC_ocean / dC_atmos
    double gamma_land = 0.0;  // dC_land / dT
    double gamma_ocean = 0.0; // dC_ocean / dT

    double beta_total() const { return beta_land + beta_ocean; }
    double gamma_total() const { return gamma_land + gamma_ocean; }
};

// dC_land = beta_land dC_atmos + gamma_land dT, and the ocean analogue.
double land_carbon_change(const CarbonFeedback& feedback, double delta_c_atmos,
                          double delta_temperature);
double ocean_carbon_change(const CarbonFeedback& feedback, double delta_c_atmos,
                           double delta_temperature);

// Airborne fraction = (dC_atmos in ppm * 2.124 PgC/ppm) / emissions in PgC.
double airborne_fraction(double atmospheric_growth_ppm, double emissions_pgc);

// Permafrost pool decay dC/dt = -k(T) C with k(T) = k_ref Q10^((T - T_ref)/10).
double permafrost_decay_rate(double k_reference, double q10, double temperature_c,
                             double reference_temperature_c = 0.0);
double permafrost_remaining_fraction(double k_reference, double q10, double temperature_c,
                                     double years, double reference_temperature_c = 0.0);

// ---------------------------------------------------------------------------
// Longwave radiative transfer (Schwarzschild equation)
// ---------------------------------------------------------------------------
// Spectral Planck radiance, W m^-2 sr^-1 m^-1.
double planck_radiance(double wavelength_m, double temperature_k);

// Invert the Planck function for a monochromatic radiance.
double brightness_temperature(double wavelength_m, double radiance);

// tau = k rho dz, transmittance = exp(-tau).
double optical_depth(double absorption_coefficient, double density, double path_length);
double transmittance(double optical_depth);

// Weighting function W = (1/g) k q exp(-tau) and contribution B * W.
double weighting_function(double absorption_coefficient, double specific_humidity,
                          double optical_depth, double gravity = 9.80665);
double contribution_function(double planck_radiance_value, double weighting);

// Isothermal layer radiance: I = B_sfc exp(-tau) + B_layer (1 - exp(-tau)).
double layer_radiance(double surface_planck, double layer_planck, double optical_depth);

// ---------------------------------------------------------------------------
// Two-stream approximation for a scattering layer over a reflecting surface
//   dF+/dtau =  g1 F+ - g2 F- - g3 w F0 exp(-tau/mu0)
//   dF-/dtau =  g2 F+ - g1 F- - g4 w F0 exp(-tau/mu0)
// ---------------------------------------------------------------------------
struct TwoStreamCoefficients {
    double g1 = 0.0;
    double g2 = 0.0;
    double g3 = 0.0;
    double g4 = 0.0;
};

// Eddington two-stream coefficients; mu0 is the cosine of the solar zenith angle.
TwoStreamCoefficients two_stream_eddington(double single_scattering_albedo, double asymmetry,
                                           double mu0 = 1.0);

struct TwoStreamResult {
    double reflectance = 0.0;           // F+(0) / (mu0 F0)
    double diffuse_transmittance = 0.0; // F-(tau*) / (mu0 F0)
    double direct_transmittance = 0.0;  // exp(-tau*/mu0)
    double total_transmittance = 0.0;   // through the bottom, after surface reflection
    double absorptance = 0.0;
};

// Solves the coupled two-stream equations on a finite-difference grid.
TwoStreamResult two_stream_layer(double optical_depth, double single_scattering_albedo,
                                 double asymmetry, double solar_zenith_deg,
                                 double surface_albedo, int levels = 256);

// ---------------------------------------------------------------------------
// Waves
// ---------------------------------------------------------------------------
double shallow_water_speed(double gravity, double depth);

struct EquatorialWave {
    double frequency = 0.0;    // rad/s
    double phase_speed = 0.0;  // m/s
};

// Kelvin mode: omega = c k.
EquatorialWave kelvin_wave(double zonal_wavenumber, double gravity, double depth);

// Equatorial inertio-gravity / Rossby dispersion relation:
//   omega^3/c^2 - omega (k^2 + (2n+1) beta / c) - beta k = 0.
std::vector<double> equatorial_wave_frequencies(double zonal_wavenumber, double beta,
                                                double gravity_wave_speed, int meridional_mode);
std::vector<EquatorialWave> equatorial_wave_modes(double zonal_wavenumber, double beta,
                                                  double gravity_wave_speed,
                                                  int meridional_mode);

// ---------------------------------------------------------------------------
// Turbulence and boundary layers
// ---------------------------------------------------------------------------
// Refractive-index spectrum phi_n(k) = 0.033 Cn2 k^(-11/3).
double refractive_index_spectrum(double cn2, double wavenumber);

// epsilon = 2 pi f^(5/2) S_u^(3/2) / (U alpha^(3/2)).
double dissipation_rate_from_spectrum(double frequency, double velocity_spectrum,
                                      double wind_speed, double kolmogorov_constant = 0.55);

// Kolmogorov microscale eta = a2 (nu^3 / epsilon)^(1/4).
double kolmogorov_microscale(double kinematic_viscosity, double dissipation_rate,
                             double a2 = 1.0);

// E(k) = alpha epsilon^(2/3) L^(5/3) (kL)^4 / (1 + (kL)^2)^(17/6).
double kolmogorov_spectrum(double wavenumber, double dissipation_rate, double outer_scale,
                           double kolmogorov_constant = 0.55);

// Stable boundary-layer depth ratio L_s / h = 1 / (c_s sqrt(mu)).
double stable_boundary_layer_ratio(double stability_parameter, double c_s = 0.62);

}  // namespace climate
