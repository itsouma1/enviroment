#include "ozone.hpp"

#include <cmath>
#include <stdexcept>

#include "linalg.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double chapman_ozone_production(double photolysis_rate, double oxygen_density) {
    require(photolysis_rate >= 0.0 && oxygen_density >= 0.0,
            "chapman_ozone_production: negative input");
    return photolysis_rate * oxygen_density;
}

double chapman_odd_oxygen_loss(double rate_constant, double oxygen_atom, double ozone) {
    require(rate_constant >= 0.0, "chapman_odd_oxygen_loss: negative rate");
    return -rate_constant * oxygen_atom * ozone;
}

double nox_odd_oxygen_loss(double rate_constant, double nitrogen_dioxide, double oxygen_atom) {
    require(rate_constant >= 0.0, "nox_odd_oxygen_loss: negative rate");
    return -2.0 * rate_constant * nitrogen_dioxide * oxygen_atom;
}

double chapman_ozone_ratio(double oxygen_photolysis, double ozone_photolysis) {
    require(ozone_photolysis > 0.0, "chapman_ozone_ratio: J_O3 must be positive");
    if (oxygen_photolysis < 0.0)
        throw std::runtime_error("chapman_ozone_ratio: J_O2 must be >= 0");
    return std::sqrt(oxygen_photolysis / ozone_photolysis);
}

double methane_fractional_release(double entry_methane, double methane) {
    require(entry_methane > 0.0, "methane_fractional_release: entry methane must be positive");
    if (methane < 0.0) throw std::runtime_error("methane_fractional_release: negative methane");
    return (entry_methane - methane) / entry_methane;
}

double age_of_air(double fractional_release) {
    const double a = fractional_release;
    return 0.3 + 15.2 * a - 21.2 * a * a + 10.4 * a * a * a;
}

double stratospheric_water_vapour(double entry_water_vapour, double entry_methane,
                                  double methane) {
    return entry_water_vapour + 2.0 * (entry_methane - methane);
}

double stratospheric_water_vapour_change(double fractional_release, double entry_methane_change) {
    return 2.0 * fractional_release * entry_methane_change;
}

double water_vapour_mass_perturbation(const std::vector<double>& water_vapour,
                                      const std::vector<double>& air_mass,
                                      double water_molar_mass, double air_molar_mass) {
    if (water_vapour.size() != air_mass.size())
        throw std::runtime_error("water_vapour_mass_perturbation: size mismatch");
    require(air_molar_mass > 0.0, "water_vapour_mass_perturbation: molar mass");
    double total = 0.0;
    for (std::size_t i = 0; i < water_vapour.size(); ++i)
        total += water_vapour[i] * air_mass[i] * water_molar_mass / air_molar_mass;
    return total;
}

double oxygen_damping_rate(double transport_damping, double oh, double ho2, double no2,
                           double clo, double bro, double a5, double a7, double b3, double d3,
                           double e3) {
    return transport_damping + a5 * oh + a7 * ho2 + 2.0 * b3 * no2 + 2.0 * d3 * clo +
           2.0 * e3 * bro;
}

double ozone_damping_rate(double transport_damping, double h, double oh, double ho2,
                          double a2, double a6, double a6b) {
    return transport_damping + a2 * h + a6 * oh + a6b * ho2;
}

double odd_oxygen_chemical_tendency(double oxygen_damping, double oxygen_atom,
                                    double ozone_damping, double ozone) {
    return -oxygen_damping * oxygen_atom - ozone_damping * ozone;
}

double effective_ox_damping(double ozone_damping, double oxygen_damping,
                            double transport_damping, double oxygen_atom, double ozone) {
    if (!(ozone > 0.0)) throw std::runtime_error("effective_ox_damping: ozone must be positive");
    return (ozone_damping - transport_damping) +
           oxygen_atom / ozone * (oxygen_damping - transport_damping);
}

double water_vapour_feedback(double radiative_flux_change, double surface_temperature_change) {
    if (!(std::fabs(surface_temperature_change) > 0.0))
        throw std::runtime_error("water_vapour_feedback: zero temperature change");
    return radiative_flux_change / surface_temperature_change;
}

double tropical_adjustment(double atmospheric_flux, double stratospheric_flux) {
    return atmospheric_flux - stratospheric_flux;
}

std::vector<double> temperature_adjustment(const std::vector<std::vector<double>>& jacobian,
                                           const std::vector<double>& heating_change) {
    if (jacobian.size() != heating_change.size())
        throw std::runtime_error("temperature_adjustment: size mismatch");
    std::vector<double> rhs(heating_change.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i) rhs[i] = -heating_change[i];
    std::vector<double> solution;
    if (!solve_linear_system(jacobian, rhs, solution))
        throw std::runtime_error("temperature_adjustment: singular Jacobian");
    return solution;
}

double methane_water_trend(double fractional_release, double entry_methane_trend) {
    return 2.0 * fractional_release * entry_methane_trend;
}

double cold_trap_entry_change(double cold_point_temperature_change) {
    return 0.45 * cold_point_temperature_change;
}

double ozone_kernel_forcing(const std::vector<double>& kernel,
                            const std::vector<double>& ozone_change) {
    if (kernel.size() != ozone_change.size())
        throw std::runtime_error("ozone_kernel_forcing: size mismatch");
    double forcing = 0.0;
    for (std::size_t i = 0; i < kernel.size(); ++i) forcing += kernel[i] * ozone_change[i];
    return forcing;
}

double stratospheric_adjusted_forcing(double shortwave, double longwave_adjusted) {
    return shortwave + longwave_adjusted;
}

}  // namespace climate
