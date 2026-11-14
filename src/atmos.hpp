#pragma once

#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Convection: CAPE and CIN
//   CAPE = int g (Tv_parcel - Tv_env)/Tv_env dz from LFC to LNB
// ---------------------------------------------------------------------------
struct ParcelLevel {
    double height_m = 0.0;
    double parcel_virtual_temperature_k = 0.0;
    double environment_virtual_temperature_k = 0.0;
};

struct CapeResult {
    double cape_j_per_kg = 0.0;
    double cin_j_per_kg = 0.0;
    double lfc_m = 0.0;  // level of free convection
    double lnb_m = 0.0;  // level of neutral buoyancy
    bool found_lfc = false;
    bool found_lnb = false;
};

// Trapezoidal integration of the buoyancy integral; levels must be sorted by height.
CapeResult cape_cin(const std::vector<ParcelLevel>& levels, double gravity = 9.80665);

// ---------------------------------------------------------------------------
// Photolysis and gas-phase chemistry
// ---------------------------------------------------------------------------
// J = integral phi(lambda) sigma(lambda) F(lambda) dlambda (trapezoidal).
double photolysis_rate(const std::vector<double>& wavelength_m,
                       const std::vector<double>& quantum_yield,
                       const std::vector<double>& cross_section_m2,
                       const std::vector<double>& actinic_flux);

// Falloff rate constant:
//   k = k_low [M] / (1 + k_low [M] / k_high) * 0.6^Z,
//   Z = (1 + (log10(k_low [M] / k_high))^2)^-1.
double falloff_rate_constant(double temperature_k, double number_density_cm3,
                             double k_low_300 = 3.3e-31, double k_high_300 = 1.6e-12,
                             double low_temperature_exponent = -4.3,
                             double high_temperature_exponent = 0.0);

// Aqueous-phase sulfate oxidation:
//   R = k0 [SO2.H2O] + k1 [HSO3-] + k2 [SO3^2-] [Oxidant]
double aqueous_sulfate_oxidation_rate(double k0, double k1, double k2,
                                      double dissolved_so2, double bisulfite,
                                      double sulfite, double oxidant);

// ---------------------------------------------------------------------------
// Aerosol-cloud interactions
// ---------------------------------------------------------------------------
// ERFaci = dCRE/dln(AI) * dln(AI), with CRE = CF * SW_cloudy so that
// dCRE/dlnAI = dCF/dlnAI * SW + CF * dSW/dlnAI.
double cloud_radiative_susceptibility(double dcre_dln_ai, double delta_ln_ai);

double aci_radiative_forcing(double mean_cloud_fraction, double mean_cloudy_shortwave,
                             double dcf_dln_ai, double dsw_dln_ai, double delta_ln_ai);

// ---------------------------------------------------------------------------
// Aerosol microphysics
// ---------------------------------------------------------------------------
// Free-molecular coagulation coefficient for radii r1, r2 and particle density rho_p.
double free_molecular_coagulation(double radius1_m, double radius2_m, double temperature_k,
                                  double particle_density_kg_m3);

// Coagulation sink CS = sum_j K(d, d_j) N_j.
double coagulation_sink(const std::vector<double>& coefficients,
                        const std::vector<double>& number_concentrations);

// ---------------------------------------------------------------------------
// Atmospheric electricity
// ---------------------------------------------------------------------------
// Steady state of dn/dt = q - alpha n^2 - n beta Z.
double ion_equilibrium_number(double ionization_rate, double recombination_coefficient,
                              double attachment_coefficient, double aerosol_concentration);

// Recombination-limited limit n = sqrt(q / alpha).
double recombination_limited_ions(double ionization_rate, double recombination_coefficient);

}  // namespace climate
