#include "aerosol.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

MieEfficiencies mie_efficiencies(double size_parameter, double refractive_real,
                                 double refractive_imaginary) {
    require(size_parameter > 0.0, "mie_efficiencies: size parameter must be positive");
    require(refractive_imaginary >= 0.0, "mie_efficiencies: imaginary index must be >= 0");
    const double x = size_parameter;
    const std::complex<double> m(refractive_real, -refractive_imaginary);

    const int order_max = static_cast<int>(std::ceil(x + 4.0 * std::cbrt(x) + 2.0));
    // D_n is the logarithmic derivative of psi_n evaluated at z = m x.
    const std::complex<double> z = m * x;
    const int d_order = order_max + 25;
    std::vector<std::complex<double>> logarithmic(d_order + 1);
    logarithmic[d_order] = std::complex<double>(0.0, 0.0);
    std::complex<double> current(0.0, 0.0);
    for (int n = d_order; n >= 1; --n) {
        const std::complex<double> n_over_z = static_cast<double>(n) / z;
        current = n_over_z - std::complex<double>(1.0, 0.0) / (current + n_over_z);
        logarithmic[n - 1] = current;  // D_{n-1} at m x
    }

    // Riccati-Bessel functions: psi by Miller's downward recurrence and chi by
    // the stable upward recurrence, normalised at order zero.
    const int psi_order = order_max + 25;
    std::vector<double> psi(psi_order + 2, 0.0);
    std::vector<double> chi(order_max + 1, 0.0);
    psi[psi_order] = 1e-30;
    for (int n = psi_order; n >= 1; --n)
        psi[n - 1] = (2.0 * n + 1.0) / x * psi[n] - psi[n + 1];
    if (std::fabs(psi[0]) > 0.0) {
        const double scale = std::sin(x) / psi[0];
        for (double& value : psi) value *= scale;
    }
    chi[0] = std::cos(x);
    if (order_max >= 1) chi[1] = std::cos(x) / x + std::sin(x);
    for (int n = 1; n < order_max; ++n)
        chi[n + 1] = (2.0 * n + 1.0) / x * chi[n] - chi[n - 1];

    double q_ext = 0.0;
    double q_sca = 0.0;
    for (int n = 1; n <= order_max; ++n) {
        const std::complex<double> psi_n(psi[n], 0.0);
        const std::complex<double> psi_prev(psi[n - 1], 0.0);
        const std::complex<double> xi_n(psi[n], chi[n]);
        const std::complex<double> xi_prev(psi[n - 1], chi[n - 1]);
        const std::complex<double> ratio_a =
            logarithmic[n] / m + std::complex<double>(static_cast<double>(n) / x, 0.0);
        const std::complex<double> ratio_b =
            m * logarithmic[n] + std::complex<double>(static_cast<double>(n) / x, 0.0);
        const std::complex<double> a = (ratio_a * psi_n - psi_prev) / (ratio_a * xi_n - xi_prev);
        const std::complex<double> b = (ratio_b * psi_n - psi_prev) / (ratio_b * xi_n - xi_prev);
        const double weight = 2.0 * n + 1.0;
        q_ext += weight * (a + b).real();
        q_sca += weight * (std::norm(a) + std::norm(b));
    }
    const double factor = 2.0 / (x * x);
    MieEfficiencies result;
    result.q_ext = factor * q_ext;
    result.q_sca = factor * q_sca;
    result.q_abs = result.q_ext - result.q_sca;
    return result;
}

MieEfficiencies mie_cross_sections(double radius, double wavelength, double refractive_real,
                                   double refractive_imaginary) {
    require(radius > 0.0 && wavelength > 0.0, "mie_cross_sections: radius and wavelength");
    const double size_parameter = 3.14159265358979323846 * 2.0 * radius / wavelength;
    MieEfficiencies efficiencies =
        mie_efficiencies(size_parameter, refractive_real, refractive_imaginary);
    const double area = 3.14159265358979323846 * radius * radius;
    efficiencies.q_ext *= area;
    efficiencies.q_sca *= area;
    efficiencies.q_abs *= area;
    return efficiencies;
}

double aerosol_optical_depth(const std::vector<double>& extinction,
                             const std::vector<double>& layer_thickness) {
    require(extinction.size() == layer_thickness.size(), "aerosol_optical_depth: size mismatch");
    double depth = 0.0;
    for (std::size_t i = 0; i < extinction.size(); ++i) {
        if (layer_thickness[i] < 0.0)
            throw std::runtime_error("aerosol_optical_depth: negative thickness");
        depth += extinction[i] * layer_thickness[i];
    }
    return depth;
}

double single_scattering_albedo(double scattering, double extinction) {
    require(extinction > 0.0, "single_scattering_albedo: extinction must be positive");
    if (scattering < 0.0) throw std::runtime_error("single_scattering_albedo: negative");
    return scattering / extinction;
}

double angstrom_optical_depth(double reference_tau, double wavelength, double reference_wavelength,
                              double angstrom_exponent) {
    require(wavelength > 0.0 && reference_wavelength > 0.0, "angstrom: wavelengths");
    return reference_tau * std::pow(wavelength / reference_wavelength, -angstrom_exponent);
}

double angstrom_exponent(double tau_1, double tau_2, double wavelength_1, double wavelength_2) {
    require(tau_1 > 0.0 && tau_2 > 0.0, "angstrom_exponent: optical depths must be positive");
    require(wavelength_1 > 0.0 && wavelength_2 > 0.0, "angstrom_exponent: wavelengths");
    return -std::log(tau_1 / tau_2) / std::log(wavelength_1 / wavelength_2);
}

double direct_aerosol_forcing(double solar_constant, double transmission, double cloud_fraction,
                              double surface_albedo, double single_scattering_albedo_value,
                              double asymmetry, double optical_depth) {
    if (cloud_fraction < 0.0 || cloud_fraction > 1.0)
        throw std::runtime_error("direct_aerosol_forcing: cloud fraction outside [0, 1]");
    const double scattering_term = single_scattering_albedo_value * (1.0 - asymmetry / 2.0) *
                                   (1.0 - surface_albedo) * (1.0 - surface_albedo);
    const double absorption_term = 4.0 * surface_albedo * (1.0 - single_scattering_albedo_value);
    return -solar_constant / 4.0 * transmission * transmission * (1.0 - cloud_fraction) *
           (scattering_term - absorption_term) * optical_depth;
}

void two_stream_derivatives(double up, double down, double optical_depth, double single_scattering,
                            double direct_flux, double mu_0, double gamma_1, double gamma_2,
                            double gamma_3, double gamma_4, double& d_up, double& d_down) {
    require(mu_0 > 0.0, "two_stream_derivatives: mu_0 must be positive");
    const double direct = direct_flux * std::exp(-optical_depth / mu_0);
    d_up = gamma_1 * up - gamma_2 * down - gamma_3 * single_scattering * direct;
    d_down = gamma_2 * up - gamma_1 * down - gamma_4 * single_scattering * direct;
}

double polarized_reflectance(double stokes_q, double stokes_u, double solar_flux,
                             double solar_zenith_cosine) {
    require(solar_flux > 0.0, "polarized_reflectance: solar flux must be positive");
    require(solar_zenith_cosine > 0.0, "polarized_reflectance: zenith cosine must be positive");
    return 3.14159265358979323846 / (solar_flux * solar_zenith_cosine) *
           std::sqrt(stokes_q * stokes_q + stokes_u * stokes_u);
}

double kohler_activated_droplets(double aerosol_number, double updraft_velocity,
                                 double activation_coefficient) {
    const double denominator = updraft_velocity + activation_coefficient * aerosol_number;
    require(denominator > 0.0, "kohler_activated_droplets: singular denominator");
    return aerosol_number * updraft_velocity / denominator;
}

double aci_index(double log_effective_radius_change, double log_aerosol_change) {
    require(std::fabs(log_aerosol_change) > 0.0, "aci_index: zero aerosol change");
    return -log_effective_radius_change / log_aerosol_change;
}

double aci_optical_depth_index(double log_optical_depth_change, double log_aerosol_change) {
    require(std::fabs(log_aerosol_change) > 0.0, "aci_optical_depth_index: zero aerosol change");
    return log_optical_depth_change / log_aerosol_change;
}

double aci_forcing_decomposition(double radiation_term, double droplet_number_term,
                                 double liquid_water_term, double cloud_fraction_term) {
    return radiation_term + droplet_number_term + liquid_water_term + cloud_fraction_term;
}

double aci_residual_twomey(double erfaci, double liquid_water_adjustment,
                           double cloud_amount_adjustment) {
    return erfaci - liquid_water_adjustment - cloud_amount_adjustment;
}

double ml_aci_forcing(double mean_downward_flux, double atmospheric_factor, double cloud_fraction,
                      double albedo_sensitivity, double droplet_sensitivity,
                      double log_aerosol_change) {
    return -mean_downward_flux * atmospheric_factor * cloud_fraction * albedo_sensitivity *
           droplet_sensitivity * log_aerosol_change;
}

}  // namespace climate
