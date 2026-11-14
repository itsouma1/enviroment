#pragma once

#include <complex>

namespace climate {

// ---------------------------------------------------------------------------
// Radiative forcing and feedback decomposition
//   N = F + lambda dT,   lambda = lambda_0 + sum_i lambda_i
// ---------------------------------------------------------------------------
struct FeedbackSet {
    double planck = -3.2;       // W m^-2 K^-1, negative (stabilising)
    double water_vapour = 1.8;  // W m^-2 K^-1, positive (amplifying)
    double lapse_rate = -0.6;   // W m^-2 K^-1
    double albedo = 0.35;       // W m^-2 K^-1
    double cloud = 0.4;         // W m^-2 K^-1

    // lambda = lambda_0 + sum_i lambda_i; defaults sum to about -1.25 W m^-2 K^-1,
    // consistent with an ECS near 3 K for F_2x = 3.7 W m^-2.
    double total() const {
        return planck + water_vapour + lapse_rate + albedo + cloud;
    }
};

// Equilibrium warming for a forcing and a net feedback parameter: dT = -F/lambda.
double equilibrium_temperature_change(double effective_forcing, double feedback_parameter);

// Equilibrium climate sensitivity for a doubled CO2 concentration.
double climate_sensitivity_2x(double feedback_parameter, double forcing_2xco2 = 3.7);

// Myhre-style logarithmic CO2 forcing: dF = 5.35 ln(C / C0).
double co2_radiative_forcing(double concentration_ppm, double reference_ppm);

// Aerosol optical depth from a monodisperse population:
//   tau = N * pi r^2 * Qext * L
double aerosol_optical_depth(double number_concentration_m3, double radius_m,
                             double extinction_efficiency, double path_length_m);

// ---------------------------------------------------------------------------
// Mie scattering by a homogeneous sphere (Bohren & Huffman)
// ---------------------------------------------------------------------------
struct MieResult {
    double qext = 0.0;  // extinction efficiency
    double qsca = 0.0;  // scattering efficiency
    double g = 0.0;     // asymmetry parameter
    double single_scattering_albedo = 0.0;
};

// size_parameter = 2 pi r / wavelength; refractive_index = n + i k with k >= 0.
MieResult mie_sphere(double size_parameter, std::complex<double> refractive_index);

// ---------------------------------------------------------------------------
// Wind stress, Ekman transport and geostrophy
// ---------------------------------------------------------------------------
double coriolis_parameter(double latitude_deg);

struct WindStress {
    double tau_x = 0.0;
    double tau_y = 0.0;

    double magnitude() const;
};

// Bulk formula: tau = rho_air * Cd * |U| * U.
WindStress bulk_wind_stress(double air_density, double drag_coefficient, double u10, double v10);

struct EkmanTransport {
    double m_x = 0.0;  // m^2/s
    double m_y = 0.0;

    double magnitude() const;
};

// M = (tau_y, -tau_x) / (rho f).
EkmanTransport ekman_transport(double tau_x, double tau_y, double water_density,
                               double coriolis);

// w_E = (1 / (rho f)) (d tau_y / dx - d tau_x / dy)
double ekman_pumping(double dtau_y_dx, double dtau_x_dy, double water_density, double coriolis);

struct GeostrophicVelocity {
    double u = 0.0;  // m/s, eastward
    double v = 0.0;  // m/s, northward
};

// f u = -(1/rho) dp/dy,  f v = (1/rho) dp/dx.
GeostrophicVelocity geostrophic_velocity(double density, double dp_dx, double dp_dy,
                                         double latitude_deg);

struct ThermalWind {
    double du_dz = 0.0;
    double dv_dz = 0.0;
};

// f du/dz = (g/rho) d rho/dy,  f dv/dz = -(g/rho) d rho/dx.
ThermalWind thermal_wind(double density, double drho_dx, double drho_dy, double latitude_deg);

// ---------------------------------------------------------------------------
// Vegetation index alternatives documented alongside NDVI
// ---------------------------------------------------------------------------
// SAVI = (NIR - RED) / (NIR + RED + L) * (1 + L)
double savi(double red_reflectance, double nir_reflectance, double soil_factor);

// EVI = G (NIR - RED) / (NIR + C1 RED - C2 BLUE + L)
double evi(double red_reflectance, double nir_reflectance, double blue_reflectance,
           double gain = 2.5, double c1 = 6.0, double c2 = 7.5, double canopy_background = 1.0);

}  // namespace climate
