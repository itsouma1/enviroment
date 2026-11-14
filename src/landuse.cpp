#include "landuse.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double mixed_albedo(const std::vector<AlbedoFraction>& components) {
    double albedo = 0.0;
    for (const auto& component : components) {
        albedo += component.land_cover_fraction * component.snow_fraction *
                  component.radiation_fraction * component.albedo;
    }
    return albedo;
}

double annual_albedo_change(const std::vector<double>& new_monthly,
                            const std::vector<double>& old_monthly) {
    require(!new_monthly.empty() && new_monthly.size() == old_monthly.size(),
            "annual_albedo_change: size mismatch");
    double change = 0.0;
    for (std::size_t i = 0; i < new_monthly.size(); ++i)
        change += new_monthly[i] - old_monthly[i];
    return change / static_cast<double>(new_monthly.size());
}

double albedo_radiative_forcing(double kernel, double albedo_change) {
    return kernel * albedo_change;
}

double clear_sky_albedo_forcing(double albedo_change, double downwelling_shortwave) {
    return -0.7 * albedo_change * downwelling_shortwave;
}

double global_conversion_forcing(const std::vector<std::vector<double>>& albedo_change_by_month,
                                 const std::vector<std::vector<double>>& converted_fraction,
                                 const std::vector<std::vector<double>>& cell_area,
                                 const std::vector<std::vector<double>>& kernel,
                                 double earth_area) {
    require(earth_area > 0.0, "global_conversion_forcing: earth area must be positive");
    require(albedo_change_by_month.size() == converted_fraction.size() &&
                albedo_change_by_month.size() == cell_area.size() &&
                albedo_change_by_month.size() == kernel.size(),
            "global_conversion_forcing: month count mismatch");
    double forcing = 0.0;
    for (std::size_t month = 0; month < albedo_change_by_month.size(); ++month) {
        require(albedo_change_by_month[month].size() == converted_fraction[month].size() &&
                    albedo_change_by_month[month].size() == cell_area[month].size() &&
                    albedo_change_by_month[month].size() == kernel[month].size(),
                "global_conversion_forcing: cell count mismatch");
        for (std::size_t cell = 0; cell < albedo_change_by_month[month].size(); ++cell) {
            forcing += albedo_change_by_month[month][cell] * converted_fraction[month][cell] *
                       cell_area[month][cell] * kernel[month][cell];
        }
    }
    return forcing / earth_area;
}

double canopy_resistance(double leaf_stomatal_resistance, double leaf_area_index) {
    if (!(leaf_area_index > 0.0))
        throw std::runtime_error("canopy_resistance: LAI must be positive");
    return leaf_stomatal_resistance / (0.5 * leaf_area_index);
}

double canopy_penman_monteith(double slope, double available_energy, double air_density,
                              double specific_heat, double vapour_pressure_deficit,
                              double aerodynamic_resistance, double psychrometric_constant,
                              double canopy_resistance_value) {
    require(aerodynamic_resistance > 0.0,
            "canopy_penman_monteith: aerodynamic resistance");
    const double numerator = slope * available_energy +
                             air_density * specific_heat * vapour_pressure_deficit /
                                 aerodynamic_resistance;
    const double denominator =
        slope + psychrometric_constant * (1.0 + canopy_resistance_value / aerodynamic_resistance);
    require(std::fabs(denominator) > 0.0, "canopy_penman_monteith: zero denominator");
    return numerator / denominator;
}

double surface_temperature_decomposition(
    double incoming_shortwave, double albedo_change, double albedo,
    double incoming_shortwave_change, double surface_emissivity, double air_temperature,
    double air_emissivity_change, double latent_heat_of_vaporisation,
    double evapotranspiration_change, double air_density, double drag_coefficient,
    double surface_temperature, double air_temperature_change, double aerodynamic_resistance,
    double aerodynamic_resistance_change, double stefan_boltzmann) {
    require(aerodynamic_resistance > 0.0,
            "surface_temperature_decomposition: aerodynamic resistance");
    const double turbulent = air_density * drag_coefficient / aerodynamic_resistance;
    const double surface_radiative = 4.0 * surface_emissivity * stefan_boltzmann *
                                     std::pow(surface_temperature, 3.0);
    const double factor = turbulent + surface_radiative;
    require(std::fabs(factor) > 0.0,
            "surface_temperature_decomposition: zero energy redistribution factor");
    const double bracket =
        -incoming_shortwave * albedo_change + (1.0 - albedo) * incoming_shortwave_change +
        surface_emissivity * stefan_boltzmann * std::pow(air_temperature, 4.0) *
            air_emissivity_change -
        latent_heat_of_vaporisation * evapotranspiration_change +
        air_density * drag_coefficient * (surface_temperature - air_temperature) /
            (aerodynamic_resistance * aerodynamic_resistance) * aerodynamic_resistance_change;
    const double air_term =
        (turbulent + 4.0 * surface_emissivity * stefan_boltzmann *
                          std::pow(air_temperature, 3.0)) /
        factor * air_temperature_change;
    return bracket / factor + air_term;
}

double bookkeeping_carbon_flux(double transition, double vegetation_carbon,
                               double soil_carbon, double vegetation_response,
                               double soil_response, double slash_response) {
    return transition * (vegetation_carbon * vegetation_response +
                         soil_carbon * soil_response + vegetation_carbon * slash_response);
}

double stock_difference(double final_stock, double initial_stock, double interval) {
    if (!(interval > 0.0)) throw std::runtime_error("stock_difference: interval must be positive");
    return (final_stock - initial_stock) / interval;
}

double immediate_emission_factor(double annual_emission_factor, double conversion_time) {
    if (!(conversion_time > 0.0))
        throw std::runtime_error("immediate_emission_factor: conversion time must be positive");
    return annual_emission_factor / conversion_time;
}

double gain_loss_stock_change(double increase, double decrease, double weight, double stock_after,
                              double stock_before, double area) {
    return (increase - decrease + weight * (stock_after - stock_before)) * area;
}

void land_use_pool_tendencies(double& vegetation, double& litter, double& soil,
                              double& product, double gpp, double vegetation_respiration,
                              double litter_flux, double wood_flux, double litter_respiration,
                              double soil_flux, double soil_respiration,
                              double product_release) {
    vegetation = gpp - vegetation_respiration - litter_flux - wood_flux;
    litter = litter_flux - litter_respiration - soil_flux;
    soil = soil_flux - soil_respiration;
    product = wood_flux - product_release;
}

double net_land_use_flux(double vegetation_star, double litter_star, double soil_star,
                         double product_star, double vegetation, double litter, double soil) {
    return (vegetation_star + litter_star + soil_star + product_star) -
           (vegetation + litter + soil);
}

}  // namespace climate
