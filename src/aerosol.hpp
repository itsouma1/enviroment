#ifndef CLIMATE_AEROSOL_HPP
#define CLIMATE_AEROSOL_HPP

#include <complex>
#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Mie theory for homogeneous spheres
// ---------------------------------------------------------------------------
struct MieEfficiencies {
    double q_ext = 0.0;
    double q_sca = 0.0;
    double q_abs = 0.0;
};

// x = pi D / lambda; the imaginary part is the absorption coefficient
// (m = mr - i mi with mi >= 0).
MieEfficiencies mie_efficiencies(double size_parameter, double refractive_real,
                                 double refractive_imaginary);

// sigma = pi r^2 Q for a single sphere.
MieEfficiencies mie_cross_sections(double radius, double wavelength, double refractive_real,
                                   double refractive_imaginary);

// Vertical integral of a polydisperse extinction profile.
double aerosol_optical_depth(const std::vector<double>& extinction,
                             const std::vector<double>& layer_thickness);

// Single-scattering albedo.
double single_scattering_albedo(double scattering, double extinction);

// tau(lambda) = tau_ref (lambda / lambda_ref)^-alpha.
double angstrom_optical_depth(double reference_tau, double wavelength, double reference_wavelength,
                              double angstrom_exponent);

// alpha from two optical depths.
double angstrom_exponent(double tau_1, double tau_2, double wavelength_1, double wavelength_2);

// Charlson-Sokolik-Toon thin-layer direct forcing.
double direct_aerosol_forcing(double solar_constant, double transmission, double cloud_fraction,
                              double surface_albedo, double single_scattering_albedo_value,
                              double asymmetry, double optical_depth);

// Two-stream right-hand sides for the diffuse fluxes.
void two_stream_derivatives(double up, double down, double optical_depth, double single_scattering,
                            double direct_flux, double mu_0, double gamma_1, double gamma_2,
                            double gamma_3, double gamma_4, double& d_up, double& d_down);

// Polarized reflectance from the Q and U Stokes components.
double polarized_reflectance(double stokes_q, double stokes_u, double solar_flux,
                             double solar_zenith_cosine);

// ---------------------------------------------------------------------------
// Aerosol-cloud interactions
// ---------------------------------------------------------------------------
// Kohler single-mode activation N = Na w / (w + d Na).
double kohler_activated_droplets(double aerosol_number, double updraft_velocity,
                                 double activation_coefficient);

// ACI = -d ln r_e / d ln alpha.
double aci_index(double log_effective_radius_change, double log_aerosol_change);

// Cloud optical-depth form at fixed liquid water path.
double aci_optical_depth_index(double log_optical_depth_change, double log_aerosol_change);

// F = F_ari + F_Nd + F_L + F_fc.
double aci_forcing_decomposition(double radiation_term, double droplet_number_term,
                                 double liquid_water_term, double cloud_fraction_term);

// Residual Twomey estimate A_re = ERFaci - A_LWP - A_AMT.
double aci_residual_twomey(double erfaci, double liquid_water_adjustment,
                           double cloud_amount_adjustment);

// Chained machine-learning forcing: -F_down phi CF (d alpha/d ln Nd)(d ln Nd/d ln AI) d ln AI.
double ml_aci_forcing(double mean_downward_flux, double atmospheric_factor, double cloud_fraction,
                      double albedo_sensitivity, double droplet_sensitivity,
                      double log_aerosol_change);

}  // namespace climate

#endif
