#include "ozone.hpp"

#include <cmath>
#include <stdexcept>

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

}  // namespace climate
