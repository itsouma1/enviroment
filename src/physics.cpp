#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "common.hpp"

namespace climate {
namespace {

constexpr double CO2_FORCING_COEFFICIENT = 5.35;  // W m^-2 per ln(ppm ratio)
constexpr double GRAVITY = 9.80665;               // m s^-2

}  // namespace

double equilibrium_temperature_change(double effective_forcing, double feedback_parameter) {
    if (feedback_parameter == 0.0)
        throw std::runtime_error("equilibrium_temperature_change: zero feedback parameter");
    return -effective_forcing / feedback_parameter;
}

double climate_sensitivity_2x(double feedback_parameter, double forcing_2xco2) {
    return equilibrium_temperature_change(forcing_2xco2, feedback_parameter);
}

double co2_radiative_forcing(double concentration_ppm, double reference_ppm) {
    if (concentration_ppm <= 0.0 || reference_ppm <= 0.0)
        throw std::runtime_error("co2_radiative_forcing: concentrations must be positive");
    return CO2_FORCING_COEFFICIENT * std::log(concentration_ppm / reference_ppm);
}

double aerosol_optical_depth(double number_concentration_m3, double radius_m,
                             double extinction_efficiency, double path_length_m) {
    if (radius_m < 0.0 || number_concentration_m3 < 0.0 || path_length_m < 0.0)
        throw std::runtime_error("aerosol_optical_depth: inputs must be non-negative");
    return number_concentration_m3 * PI * radius_m * radius_m * extinction_efficiency *
           path_length_m;
}

MieResult mie_sphere(double size_parameter, std::complex<double> refractive_index) {
    const double x = size_parameter;
    if (!(x > 0.0)) throw std::runtime_error("mie_sphere: size parameter must be positive");
    if (refractive_index.imag() < 0.0) refractive_index = std::conj(refractive_index);

    const int nmax = std::max(2, static_cast<int>(std::floor(x + 4.0 * std::cbrt(x) + 2.0)));
    const std::complex<double> mx = refractive_index * x;
    const int nmx = std::max(nmax, static_cast<int>(std::abs(mx))) + 16;

    // Downward logarithmic-derivative recurrence, stable for absorbing spheres.
    std::vector<std::complex<double>> d(static_cast<std::size_t>(nmx) + 2,
                                        std::complex<double>(0.0, 0.0));
    for (int n = nmx; n >= 1; --n) {
        const std::complex<double> ratio = static_cast<double>(n) / mx;
        d[static_cast<std::size_t>(n - 1)] = ratio - 1.0 / (d[static_cast<std::size_t>(n)] + ratio);
    }

    std::vector<double> psi(static_cast<std::size_t>(nmax) + 2, 0.0);
    std::vector<double> chi(static_cast<std::size_t>(nmax) + 2, 0.0);
    psi[0] = std::sin(x);
    psi[1] = std::sin(x) / x - std::cos(x);
    chi[0] = std::cos(x);
    chi[1] = std::cos(x) / x + std::sin(x);
    for (int n = 1; n <= nmax; ++n) {
        const double factor = (2.0 * n + 1.0) / x;
        psi[static_cast<std::size_t>(n + 1)] =
            factor * psi[static_cast<std::size_t>(n)] - psi[static_cast<std::size_t>(n - 1)];
        chi[static_cast<std::size_t>(n + 1)] =
            factor * chi[static_cast<std::size_t>(n)] - chi[static_cast<std::size_t>(n - 1)];
    }

    std::vector<std::complex<double>> a(static_cast<std::size_t>(nmax) + 2,
                                        std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> b(static_cast<std::size_t>(nmax) + 2,
                                        std::complex<double>(0.0, 0.0));
    for (int n = 1; n <= nmax; ++n) {
        const std::size_t index = static_cast<std::size_t>(n);
        // xi_n = psi_n - i chi_n (Bohren & Huffman), the e^{-i omega t} convention.
        const std::complex<double> xi_n(psi[index], -chi[index]);
        const std::complex<double> xi_previous(psi[index - 1], -chi[index - 1]);
        const std::complex<double> da = d[index] / refractive_index + static_cast<double>(n) / x;
        const std::complex<double> db = refractive_index * d[index] + static_cast<double>(n) / x;
        a[index] = (da * psi[index] - psi[index - 1]) / (da * xi_n - xi_previous);
        b[index] = (db * psi[index] - psi[index - 1]) / (db * xi_n - xi_previous);
    }

    double extinction_sum = 0.0;
    double scattering_sum = 0.0;
    double asymmetry_sum = 0.0;
    for (int n = 1; n <= nmax; ++n) {
        const std::size_t index = static_cast<std::size_t>(n);
        const double factor = 2.0 * n + 1.0;
        extinction_sum += factor * (a[index].real() + b[index].real());
        scattering_sum += factor * (std::norm(a[index]) + std::norm(b[index]));

        asymmetry_sum += (static_cast<double>(n) * (n + 2.0) / (n + 1.0)) *
                         (a[index] * std::conj(b[index])).real();
        if (n < nmax) {
            asymmetry_sum += ((2.0 * n + 1.0) / (static_cast<double>(n) * (n + 1.0))) *
                             (a[index] * std::conj(b[index + 1])).real();
        }
    }

    MieResult result;
    result.qext = 2.0 / (x * x) * extinction_sum;
    result.qsca = 2.0 / (x * x) * scattering_sum;
    result.g = result.qsca > 0.0 ? 4.0 / (x * x * result.qsca) * asymmetry_sum : 0.0;
    result.single_scattering_albedo = result.qext > 0.0 ? result.qsca / result.qext : 0.0;
    result.qext = std::max(0.0, result.qext);
    result.qsca = std::max(0.0, std::min(result.qsca, result.qext));
    return result;
}

double coriolis_parameter(double latitude_deg) {
    const double omega = 7.2921159e-5;  // rad/s
    return 2.0 * omega * std::sin(deg_to_rad(latitude_deg));
}

double WindStress::magnitude() const { return std::hypot(tau_x, tau_y); }

WindStress bulk_wind_stress(double air_density, double drag_coefficient, double u10, double v10) {
    const double speed = std::hypot(u10, v10);
    WindStress stress;
    stress.tau_x = air_density * drag_coefficient * speed * u10;
    stress.tau_y = air_density * drag_coefficient * speed * v10;
    return stress;
}

double EkmanTransport::magnitude() const { return std::hypot(m_x, m_y); }

EkmanTransport ekman_transport(double tau_x, double tau_y, double water_density,
                               double coriolis) {
    if (water_density <= 0.0)
        throw std::runtime_error("ekman_transport: density must be positive");
    if (coriolis == 0.0)
        throw std::runtime_error("ekman_transport: Coriolis parameter is zero at the equator");
    EkmanTransport transport;
    transport.m_x = tau_y / (water_density * coriolis);
    transport.m_y = -tau_x / (water_density * coriolis);
    return transport;
}

double ekman_pumping(double dtau_y_dx, double dtau_x_dy, double water_density, double coriolis) {
    if (water_density <= 0.0 || coriolis == 0.0)
        throw std::runtime_error("ekman_pumping: density and Coriolis parameter must be non-zero");
    return (dtau_y_dx - dtau_x_dy) / (water_density * coriolis);
}

GeostrophicVelocity geostrophic_velocity(double density, double dp_dx, double dp_dy,
                                         double latitude_deg) {
    const double f = coriolis_parameter(latitude_deg);
    if (density <= 0.0) throw std::runtime_error("geostrophic_velocity: density must be positive");
    if (f == 0.0)
        throw std::runtime_error("geostrophic_velocity: geostrophy is undefined at the equator");
    GeostrophicVelocity velocity;
    velocity.u = -dp_dy / (density * f);
    velocity.v = dp_dx / (density * f);
    return velocity;
}

ThermalWind thermal_wind(double density, double drho_dx, double drho_dy, double latitude_deg) {
    const double f = coriolis_parameter(latitude_deg);
    if (density <= 0.0) throw std::runtime_error("thermal_wind: density must be positive");
    if (f == 0.0) throw std::runtime_error("thermal_wind: undefined at the equator");
    ThermalWind shear;
    shear.du_dz = GRAVITY * drho_dy / (f * density);
    shear.dv_dz = -GRAVITY * drho_dx / (f * density);
    return shear;
}

double savi(double red_reflectance, double nir_reflectance, double soil_factor) {
    const double denominator = nir_reflectance + red_reflectance + soil_factor;
    if (denominator == 0.0) return std::nan("");
    return (nir_reflectance - red_reflectance) / denominator * (1.0 + soil_factor);
}

double evi(double red_reflectance, double nir_reflectance, double blue_reflectance, double gain,
           double c1, double c2, double canopy_background) {
    const double denominator = nir_reflectance + c1 * red_reflectance - c2 * blue_reflectance +
                               canopy_background;
    if (denominator == 0.0) return std::nan("");
    return gain * (nir_reflectance - red_reflectance) / denominator;
}

}  // namespace climate
