#include "atmos.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <string>

namespace climate {
namespace {

void require_same_size(std::initializer_list<std::size_t> sizes, const char* name) {
    const std::size_t first = *sizes.begin();
    for (const std::size_t size : sizes) {
        if (size != first)
            throw std::runtime_error(std::string(name) + ": inputs must have equal length");
    }
}

}  // namespace

CapeResult cape_cin(const std::vector<ParcelLevel>& levels, double gravity) {
    if (levels.size() < 2)
        throw std::runtime_error("cape_cin: need at least two levels");
    for (std::size_t i = 1; i < levels.size(); ++i) {
        if (levels[i].height_m <= levels[i - 1].height_m)
            throw std::runtime_error("cape_cin: levels must be sorted by increasing height");
    }

    CapeResult result;
    bool in_positive = false;
    bool seen_lfc = false;
    for (std::size_t i = 0; i + 1 < levels.size(); ++i) {
        const double dz = levels[i + 1].height_m - levels[i].height_m;
        const double parcel_mid =
            0.5 * (levels[i].parcel_virtual_temperature_k +
                   levels[i + 1].parcel_virtual_temperature_k);
        const double environment_mid =
            0.5 * (levels[i].environment_virtual_temperature_k +
                   levels[i + 1].environment_virtual_temperature_k);
        if (environment_mid <= 0.0)
            throw std::runtime_error("cape_cin: environment temperature must be positive");
        const double buoyancy = gravity * (parcel_mid - environment_mid) / environment_mid;
        const double contribution = buoyancy * dz;

        if (buoyancy > 0.0) {
            if (!in_positive) {
                in_positive = true;
                if (!seen_lfc) {
                    seen_lfc = true;
                    result.found_lfc = true;
                    result.lfc_m = levels[i].height_m;
                }
            }
            result.cape_j_per_kg += contribution;
        } else {
            if (in_positive) {
                // The parcel has reached its level of neutral buoyancy.
                in_positive = false;
                result.found_lnb = true;
                result.lnb_m = levels[i].height_m;
            } else if (!seen_lfc) {
                result.cin_j_per_kg += contribution;
            }
        }
    }
    return result;
}

double photolysis_rate(const std::vector<double>& wavelength_m,
                       const std::vector<double>& quantum_yield,
                       const std::vector<double>& cross_section_m2,
                       const std::vector<double>& actinic_flux) {
    require_same_size({wavelength_m.size(), quantum_yield.size(), cross_section_m2.size(),
                       actinic_flux.size()},
                      "photolysis_rate");
    if (wavelength_m.size() < 2) throw std::runtime_error("photolysis_rate: need at least two bins");
    double total = 0.0;
    for (std::size_t i = 0; i + 1 < wavelength_m.size(); ++i) {
        const double dlambda = wavelength_m[i + 1] - wavelength_m[i];
        if (dlambda <= 0.0)
            throw std::runtime_error("photolysis_rate: wavelengths must increase");
        const auto product = [&](std::size_t j) {
            return quantum_yield[j] * cross_section_m2[j] * actinic_flux[j];
        };
        total += 0.5 * (product(i) + product(i + 1)) * dlambda;
    }
    return total;
}

double falloff_rate_constant(double temperature_k, double number_density_cm3,
                             double k_low_300, double k_high_300,
                             double low_temperature_exponent,
                             double high_temperature_exponent) {
    if (temperature_k <= 0.0)
        throw std::runtime_error("falloff_rate_constant: temperature must be positive");
    if (number_density_cm3 < 0.0)
        throw std::runtime_error("falloff_rate_constant: [M] must be non-negative");

    const double ratio = temperature_k / 300.0;
    const double k_low = k_low_300 * std::pow(ratio, low_temperature_exponent) *
                         number_density_cm3;
    const double k_high = k_high_300 * std::pow(ratio, high_temperature_exponent);
    if (k_high == 0.0)
        throw std::runtime_error("falloff_rate_constant: k_high must be non-zero");
    if (k_low == 0.0) return 0.0;

    const double logarithm = std::log10(k_low / k_high);
    const double z = 1.0 / (1.0 + logarithm * logarithm);
    return k_low / (1.0 + k_low / k_high) * std::pow(0.6, z);
}

double aqueous_sulfate_oxidation_rate(double k0, double k1, double k2,
                                      double dissolved_so2, double bisulfite,
                                      double sulfite, double oxidant) {
    if (dissolved_so2 < 0.0 || bisulfite < 0.0 || sulfite < 0.0 || oxidant < 0.0)
        throw std::runtime_error("aqueous_sulfate_oxidation_rate: concentrations must be >= 0");
    return k0 * dissolved_so2 + k1 * bisulfite + k2 * sulfite * oxidant;
}

double cloud_radiative_susceptibility(double dcre_dln_ai, double delta_ln_ai) {
    return dcre_dln_ai * delta_ln_ai;
}

double aci_radiative_forcing(double mean_cloud_fraction, double mean_cloudy_shortwave,
                             double dcf_dln_ai, double dsw_dln_ai, double delta_ln_ai) {
    const double susceptibility = dcf_dln_ai * mean_cloudy_shortwave +
                                  mean_cloud_fraction * dsw_dln_ai;
    return susceptibility * delta_ln_ai;
}

double free_molecular_coagulation(double radius1_m, double radius2_m, double temperature_k,
                                  double particle_density_kg_m3) {
    if (radius1_m <= 0.0 || radius2_m <= 0.0)
        throw std::runtime_error("free_molecular_coagulation: radii must be positive");
    if (temperature_k <= 0.0)
        throw std::runtime_error("free_molecular_coagulation: temperature must be positive");
    if (particle_density_kg_m3 <= 0.0)
        throw std::runtime_error("free_molecular_coagulation: density must be positive");
    constexpr double BOLTZMANN = 1.380649e-23;
    const double thermal = std::sqrt(6.0 * BOLTZMANN * temperature_k / particle_density_kg_m3);
    const double inverse_cubes =
        std::sqrt(1.0 / std::pow(radius1_m, 3.0) + 1.0 / std::pow(radius2_m, 3.0));
    const double cross = (radius1_m + radius2_m) * (radius1_m + radius2_m);
    return thermal * inverse_cubes * cross;
}

double coagulation_sink(const std::vector<double>& coefficients,
                        const std::vector<double>& number_concentrations) {
    require_same_size({coefficients.size(), number_concentrations.size()}, "coagulation_sink");
    double sink = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        if (number_concentrations[i] < 0.0)
            throw std::runtime_error("coagulation_sink: concentrations must be >= 0");
        sink += coefficients[i] * number_concentrations[i];
    }
    return sink;
}

double ion_equilibrium_number(double ionization_rate, double recombination_coefficient,
                              double attachment_coefficient, double aerosol_concentration) {
    if (ionization_rate < 0.0) throw std::runtime_error("ion_equilibrium_number: q must be >= 0");
    if (recombination_coefficient <= 0.0)
        throw std::runtime_error("ion_equilibrium_number: alpha must be positive");
    if (attachment_coefficient < 0.0 || aerosol_concentration < 0.0)
        throw std::runtime_error("ion_equilibrium_number: attachment inputs must be >= 0");

    // q - alpha n^2 - beta Z n = 0
    const double linear = attachment_coefficient * aerosol_concentration;
    const double discriminant = linear * linear + 4.0 * recombination_coefficient * ionization_rate;
    return (-linear + std::sqrt(discriminant)) / (2.0 * recombination_coefficient);
}

double recombination_limited_ions(double ionization_rate, double recombination_coefficient) {
    if (ionization_rate < 0.0)
        throw std::runtime_error("recombination_limited_ions: q must be >= 0");
    if (recombination_coefficient <= 0.0)
        throw std::runtime_error("recombination_limited_ions: alpha must be positive");
    return std::sqrt(ionization_rate / recombination_coefficient);
}

}  // namespace climate
