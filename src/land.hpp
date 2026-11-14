#pragma once

namespace climate {

// ---------------------------------------------------------------------------
// Land surface energy balance
// ---------------------------------------------------------------------------
// Rn = (1 - albedo) S_in + eps_s sigma (eps_a T_a^4 - T_s^4)
double land_net_radiation(double albedo, double shortwave_down, double emissivity_surface,
                          double emissivity_air, double air_temperature_k,
                          double surface_temperature_k);

// H = rho cp (T_s - T_a) / r_a
double sensible_heat_flux(double air_density, double specific_heat, double surface_temperature_k,
                          double air_temperature_k, double aerodynamic_resistance);

// r_a = 1/(k u*) [ ln((z_w - d)/z_oh) - psi_h((z_w - d)/L) ]
double aerodynamic_resistance(double von_karman, double friction_velocity, double measurement_height,
                              double displacement_height, double scalar_roughness_length,
                              double stability_correction);

// ---------------------------------------------------------------------------
// Penman-Monteith evapotranspiration
// ---------------------------------------------------------------------------
double penman_monteith_latent_heat(double delta, double net_radiation_minus_ground,
                                   double air_density, double specific_heat,
                                   double humidity_deficit, double aerodynamic_resistance,
                                   double psychrometric_constant, double surface_resistance);

double penman_monteith_sensible_heat(double delta, double net_radiation_minus_ground,
                                     double air_density, double specific_heat,
                                     double humidity_deficit, double aerodynamic_resistance,
                                     double psychrometric_constant, double surface_resistance);

// r_s = r_min ((theta_sat - theta_res) / (theta_top - theta_res))^3
double surface_resistance_soil(double minimum_resistance, double saturation_moisture,
                               double residual_moisture, double top_layer_moisture);

// ---------------------------------------------------------------------------
// Photosynthesis and stomatal conductance
// ---------------------------------------------------------------------------
// Rubisco-limited: Ac = Vcmax (Ci - Gamma*) / (Ci + Kc (1 + Oi/Ko))
double farquhar_rubisco(double vcmax, double intercellular_co2, double compensation_point,
                        double kc, double intercellular_o2, double ko);

// Electron-transport-limited: Aj = J (Ci - Gamma*) / (4 Ci + 8 Gamma*)
double farquhar_electron_transport(double electron_transport, double intercellular_co2,
                                   double compensation_point);

// An = min(Ac, Aj) - Rd
double farquhar_assimilation(double rubisco_rate, double electron_transport_rate,
                             double dark_respiration);

// Ball-Berry-Leuning: gs = g1 An Patm / ((Cs - Gamma*) (1 + VPD/D0))
double ball_berry_conductance(double g1, double assimilation, double atmospheric_pressure,
                              double surface_co2, double compensation_point, double vpd,
                              double vpd_scale);

// Medlyn: gs = (1 + g1/sqrt(VPD)) An / Cs
double medlyn_conductance(double g1, double assimilation, double surface_co2, double vpd);

// Leuning (1995): gs = g0 + a1 An/(Ct - Gamma*) f(VPD) Patmr
double leuning_conductance(double cuticular_conductance, double a1, double assimilation,
                           double leaf_internal_co2, double compensation_point,
                           double vpd_response, double pressure_scaling);

// ---------------------------------------------------------------------------
// Carbon balance
// ---------------------------------------------------------------------------
// LPJ: NPP = 0.75 (GPP - Rm)
double npp_lpj(double gpp, double maintenance_respiration);

// Rm = c_RES,N Q10^((T - 20)/10) N_L [(C_L + C_R) + c_sapwood (C_WL + C_WR)]
double maintenance_respiration(double respiration_coefficient, double q10, double temperature_c,
                               double leaf_nitrogen, double leaf_carbon, double root_carbon,
                               double sapwood_coefficient, double wood_leaf_carbon,
                               double wood_root_carbon);

// LAI = C_leaf * SLA
double lai_from_leaf_carbon(double leaf_carbon, double specific_leaf_area);

// ---------------------------------------------------------------------------
// Plant hydraulics and soil water
// ---------------------------------------------------------------------------
// G = 1 - 1/(1 + exp(a (psi_soil - P50)))
double xylem_conductance(double soil_water_potential, double p50, double shape = 5.0);

// dW = (P - E - R - D) dt
double soil_water_step(double storage, double precipitation, double evaporation, double runoff,
                       double drainage, double dt);

// I = I_ds + (theta - theta_res) * thermal_scaling
double soil_thermal_inertia(double dry_soil_inertia, double volumetric_moisture,
                            double residual_moisture, double thermal_scaling);

// Arrhenius temperature response of Vcmax: V25 exp(dH (T - T0) / (T0 R T))
double vcmax_temperature(double vcmax_25, double temperature_c, double activation_energy = 65330.0);

// Dynamic vegetation and advanced land surface (notes 128-144 KB) ----------

// B_turnover = B_l / (29664 (A_SL * 10)^-0.909).
double leaf_turnover(double specific_leaf_area, double leaf_biomass);

// E_i = f_top E_0 + f_root(i) E_1.
double root_uptake(double top_fraction, double surface_potential, double root_fraction,
                   double root_potential);

// T = (theta/theta_s)^n f_c p.
double transpiration_moisture(double moisture, double saturation, double exponent,
                              double canopy_fraction, double potential);

// A_net (1 - E/E_crit), the hydraulic-penalty stomatal objective.
double optimal_stomatal_objective(double assimilation, double transpiration,
                                  double critical_transpiration);

// F = A - alpha J_max - gamma dpsi^2.
double hydraulic_profit(double assimilation, double jmax, double carbon_cost,
                        double hydraulic_cost, double water_potential_drop);

// p = a + b t_growth.
double growth_temperature_acclimation(double intercept, double slope,
                                      double growth_temperature_c);

// L_m = 1 - (A_a - g_m)/(A_i - g_m), expressed in percent.
double mesophyll_limitation(double assimilation_ambient, double assimilation_intercellular,
                            double mesophyll_conductance);

// L_in = min(L_N, L_P).
double nutrient_colimitation(double nitrogen_limitation, double phosphorus_limitation);

// max(C_net / CN - N_retrans, 0).
double nutrient_demand(double net_carbon, double stoichiometric_ratio, double retranslocation);

// dot M_f = r_m rho_f M_f; the daily accumulation is the sum over time steps.
double nitrogen_fixation(double rate, double density, double mass, double dt);

// Canopy air-space exchange E' = -rho_s (q_s - q')/r_q'.
double canopy_air_exchange(double air_density, double canopy_specific_humidity,
                           double air_specific_humidity, double resistance);

// Multilayer scalar conservation: dz dC/dt = -(F_up - F_down) + S dz.
double multilayer_scalar_tendency(double flux_up, double flux_down, double source,
                                  double layer_thickness);

// Soil production P = P_0 exp(-h/h_p*).
double soil_production(double reference_rate, double regolith_thickness, double decay_depth);

// P_0 = k_rp q exp[-E/R (1/T - 1/T0)].
double soil_production_reference(double coefficient, double runoff, double activation_energy,
                                 double temperature_k, double reference_temperature = 288.15);

}  // namespace climate
