#pragma once

#include <vector>

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

// Green-Ampt infiltration and perirhizal hydraulics -------------------------

// f = Ks (1 + S (theta_f - theta_i) / F).
double green_ampt_infiltration(double saturated_conductivity, double suction_head,
                               double final_moisture, double initial_moisture,
                               double cumulative_infiltration);

// df/dF = -Ks S (theta_f - theta_i) / F^2.
double green_ampt_infiltration_sensitivity(double saturated_conductivity, double suction_head,
                                           double final_moisture, double initial_moisture,
                                           double cumulative_infiltration);

// t_p = S Ks (theta_f - theta_i) / (P (P - Ks)).
double green_ampt_ponding_time(double saturated_conductivity, double suction_head,
                               double final_moisture, double initial_moisture, double rainfall);

// B = 2(p^2 - 1) / (1 - (0.53 p)^2 + 2 p^2 [ln p + ln 0.53]).
double perirhizal_geometry(double radius_ratio);

// Mean conductivity between two matric potentials.
double perirhizal_effective_conductivity(double cumulative_at_bulk, double cumulative_at_root,
                                         double bulk_potential, double root_potential);

// K_prhiz = 2 pi l_root B k_eff.
double perirhizal_conductance(double root_length_density, double geometry_factor,
                              double effective_conductivity);

// K_srs = K_rs K_prhiz / (SUF K_rs + K_prhiz).
double root_system_conductance(double root_system_conductance, double perirhizal_conductance,
                               double standard_uptake_fraction);

// dK_srs/dK_prhiz.
double root_system_conductance_sensitivity(double root_system_conductance,
                                           double perirhizal_conductance,
                                           double standard_uptake_fraction);

// Multilayer canopy ---------------------------------------------------------

// F = (C_i - C_{i+1}) / r.
double canopy_layer_flux(double concentration_i, double concentration_next, double resistance);

// A_n,i = -(LAI/SF) dC / (r_b + 1/(SF/tau_cut + g_s)).
double canopy_layer_assimilation(double leaf_area_index, double sidedness_factor,
                                 double ambient_co2, double intercellular_co2,
                                 double boundary_resistance, double cuticular_resistance,
                                 double stomatal_conductance);

// dA_n/dg_s.
double canopy_layer_assimilation_sensitivity(double leaf_area_index, double sidedness_factor,
                                             double ambient_co2, double intercellular_co2,
                                             double boundary_resistance,
                                             double cuticular_resistance,
                                             double stomatal_conductance);

// Critical zone (notes 144-160 KB) -----------------------------------------

// K = kd (1 - exp(-kw q)) exp[-Ea/R (1/T - 1/T0)].
double mineral_dissolution_constant(double dissolution_coefficient, double runoff,
                                    double saturation_parameter, double activation_energy,
                                    double temperature_k, double reference_temperature = 288.15);

// F = chi int K tau^sigma x_p dz, trapezoidal over the column.
double silicate_weathering_flux(double calcium_magnesium_fraction,
                                const std::vector<double>& dissolution,
                                const std::vector<double>& particle_age,
                                const std::vector<double>& mineral_fraction, double dz,
                                double age_exponent, double reference_age = 1.0);

// pCO2_max = 1 + 0.302484 q^0.8.
double soil_pco2_maximum(double runoff);

// pCO2 = atm + max / (1 + exp(1.315 - 0.116 (T - 273.15))).
double soil_pco2(double atmospheric_pco2, double maximum_soil_pco2, double temperature_k);

// F_carw = k_carb x_L q [Ca2+]_eq.
double carbonate_weathering_flux(double rate_constant, double lithology_fraction,
                                 double runoff, double equilibrium_calcium);

// E = ke q^0.5 S.
double stream_power_erosion(double erodibility, double runoff, double slope);

// F_focw = 0.5 sum_l x_L chi_foc E.
double petrogenic_carbon_flux(double lithology_fraction, double bedrock_fraction,
                              double erosion);

// F_sulw = sum_l x_L chi_S E.
double sulfide_weathering_flux(double lithology_fraction, double sulfur_fraction,
                               double erosion);

// U = Vmax DOC/(K_DOC + DOC) O2/(K_O2 + O2) B.
double microbial_uptake(double maximum_rate, double dissolved_organic_carbon,
                        double doc_half_saturation, double oxygen, double oxygen_half_saturation,
                        double active_biomass);

// Vmax(T) = Vref exp[-Ea/R (1/T - 1/Tref)].
double microbial_rate_temperature(double reference_rate, double activation_energy,
                                  double temperature_k, double reference_temperature = 288.15);

// MAOC_max = Ceq + (C0 - Ceq) exp(-k t).
double maoc_saturation(double equilibrium_carbon, double initial_carbon, double turnover,
                       double time);

// 1 - (0.17 + 0.68 sand) - (0.003 + 0.032 clay).
double soil_microbial_stabilization(double sand_fraction, double clay_fraction);

// F_arc = (1 - chi_sub,c) F_sub.
double arc_volcanic_flux(double subducted_flux, double subducted_fraction);

// F_degas = f_UM chi_degas C_mantle.
double mantle_degassing_flux(double melt_fraction, double degassing_efficiency,
                             double mantle_carbon);

// P_org = C_org / CPR.
double organic_phosphorus(double organic_carbon, double carbon_phosphorus_ratio);

// P_in = P_in,org C_roots/(C_roots + K_n).
double plant_phosphorus_uptake(double organic_phosphorus_input, double root_carbon,
                               double half_saturation);

// EPIC partition fractions.
double epic_stable_phosphorus_fraction(double saturation);
double epic_labile_phosphorus_loss(double saturation);

// s_W = sqrt(D_eff t).
double weathering_front_thickness(double effective_diffusivity, double time);

// Sensitivities of the empirical laws --------------------------------------

// d(pCO2_max)/dq = 0.8 * 0.302484 q^-0.2 (atm per runoff unit).
double soil_pco2_maximum_sensitivity(double runoff);

// dE/dq = 0.5 ke q^-0.5 S and dE/dS = ke q^0.5.
double stream_power_erosion_runoff_sensitivity(double erodibility, double runoff,
                                               double slope);
double stream_power_erosion_slope_sensitivity(double erodibility, double runoff);

// dP_0/dq = k_rp exp[-E/R (1/T - 1/T0)].
double soil_production_reference_sensitivity(double coefficient, double activation_energy,
                                             double temperature_k,
                                             double reference_temperature = 288.15);

// van Genuchten-Mualem closure ---------------------------------------------

// Se = [1 + (alpha |h|)^n]^-m with m = 1 - 1/n.
double van_genuchten_effective_saturation(double matric_potential, double alpha,
                                          double pore_size_index);

// K = Ks Se^l [1 - (1 - Se^(1/m))^m]^2.
double van_genuchten_conductivity(double matric_potential, double saturated_conductivity,
                                  double alpha, double pore_size_index, double tortuosity = 0.5);

// Closed-form matric flux potential (cumulative from h to saturation).
double van_genuchten_matric_flux_potential(double matric_potential,
                                           double saturated_conductivity, double alpha,
                                           double pore_size_index);

// K_prhiz = 2 pi l_root B (Phi(bs) - Phi(sr)) / (bs - sr).
double perirhizal_conductance_van_genuchten(double root_length_density, double geometry_factor,
                                            double bulk_potential, double root_potential,
                                            double saturated_conductivity, double alpha,
                                            double pore_size_index);

}  // namespace climate
