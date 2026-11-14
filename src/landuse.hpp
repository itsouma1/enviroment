#ifndef CLIMATE_LANDUSE_HPP
#define CLIMATE_LANDUSE_HPP

#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Albedo
// ---------------------------------------------------------------------------
struct AlbedoFraction {
    double land_cover_fraction = 0.0;
    double snow_fraction = 0.0;
    double radiation_fraction = 0.0;
    double albedo = 0.0;
};

// alpha = sum f_l f_s f_r alpha_{l,s,r}.
double mixed_albedo(const std::vector<AlbedoFraction>& components);

// Annual mean of the monthly albedo differences.
double annual_albedo_change(const std::vector<double>& new_monthly,
                            const std::vector<double>& old_monthly);

// Delta F = K Delta alpha.
double albedo_radiative_forcing(double kernel, double albedo_change);

// RFAlb_clear = -0.7 Delta alpha SWDW_TOA.
double clear_sky_albedo_forcing(double albedo_change, double downwelling_shortwave);

// RF = (1/A_earth) sum_m sum_i dalpha * converted_fraction * area * kernel.
double global_conversion_forcing(const std::vector<std::vector<double>>& albedo_change_by_month,
                                 const std::vector<std::vector<double>>& converted_fraction,
                                 const std::vector<std::vector<double>>& cell_area,
                                 const std::vector<std::vector<double>>& kernel,
                                 double earth_area);

// ---------------------------------------------------------------------------
// Evapotranspiration and surface temperature
// ---------------------------------------------------------------------------
// r_c = r_xi / (0.5 LAI).
double canopy_resistance(double leaf_stomatal_resistance, double leaf_area_index);

// Canopy-scale Penman-Monteith latent heat flux.
double canopy_penman_monteith(double slope, double available_energy, double air_density,
                              double specific_heat, double vapour_pressure_deficit,
                              double aerodynamic_resistance, double psychrometric_constant,
                              double canopy_resistance_value);

// Eight-term surface temperature decomposition.
double surface_temperature_decomposition(
    double incoming_shortwave, double albedo_change, double albedo,
    double incoming_shortwave_change, double surface_emissivity,
    double air_temperature, double air_emissivity_change, double latent_heat_of_vaporisation,
    double evapotranspiration_change, double air_density, double drag_coefficient,
    double surface_temperature, double air_temperature_change, double aerodynamic_resistance,
    double aerodynamic_resistance_change, double stefan_boltzmann = 5.670374419e-8);

// ---------------------------------------------------------------------------
// Carbon accounting
// ---------------------------------------------------------------------------
// Bookkeeping flux for one transition.
double bookkeeping_carbon_flux(double transition, double vegetation_carbon,
                               double soil_carbon, double vegetation_response,
                               double soil_response, double slash_response);

// Delta C = (C2 - C1)/(t2 - t1).
double stock_difference(double final_stock, double initial_stock, double interval);

// IEF_20y = EF_annual / 20.
double immediate_emission_factor(double annual_emission_factor, double conversion_time = 20.0);

// Gain-loss stock change with effectiveness weights.
double gain_loss_stock_change(double increase, double decrease, double weight, double stock_after,
                              double stock_before, double area);

// Pool tendencies under land-use change.
void land_use_pool_tendencies(double& vegetation, double& litter, double& soil,
                              double& product, double gpp, double vegetation_respiration,
                              double litter_flux, double wood_flux, double litter_respiration,
                              double soil_flux, double soil_respiration,
                              double product_release);

// Net land-use flux from the perturbed and control pools.
double net_land_use_flux(double vegetation_star, double litter_star, double soil_star,
                         double product_star, double vegetation, double litter, double soil);

}  // namespace climate

#endif
