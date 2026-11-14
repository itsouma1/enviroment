#ifndef CLIMATE_OZONE_HPP
#define CLIMATE_OZONE_HPP

#include <vector>

namespace climate {

// Chapman cycle -------------------------------------------------------------

// P(O3) = J_O2 [O2].
double chapman_ozone_production(double photolysis_rate, double oxygen_density);

// d([O] + [O3])/dt = -k4 [O][O3].
double chapman_odd_oxygen_loss(double rate_constant, double oxygen_atom, double ozone);

// Catalytic cycles ----------------------------------------------------------

// NOx: d([O] + [O3])/dt = -2 k6 [NO2][O].
double nox_odd_oxygen_loss(double rate_constant, double nitrogen_dioxide, double oxygen_atom);

// Pure-Chapman equilibrium scaling [O3] proportional to sqrt(J_O2 / J_O3).
double chapman_ozone_ratio(double oxygen_photolysis, double ozone_photolysis);

// Fractional release of methane --------------------------------------------

// alpha = (CH4_entry - CH4)/(CH4_entry).
double methane_fractional_release(double entry_methane, double methane);

// AoA = 0.3 + 15.2 alpha - 21.2 alpha^2 + 10.4 alpha^3.
double age_of_air(double fractional_release);

// H2O = entry + 2 (CH4_entry - CH4).
double stratospheric_water_vapour(double entry_water_vapour, double entry_methane,
                                  double methane);

// Delta SWV = 2 alpha Delta CH4_entry.
double stratospheric_water_vapour_change(double fractional_release, double entry_methane_change);

// Column mass perturbation: sum SWV m_air M_H2O / M_air.
double water_vapour_mass_perturbation(const std::vector<double>& water_vapour,
                                      const std::vector<double>& air_mass,
                                      double water_molar_mass, double air_molar_mass);

// Linear damping closure (Chapman + 2) -------------------------------------

double oxygen_damping_rate(double transport_damping, double oh, double ho2, double no2,
                           double clo, double bro, double a5, double a7, double b3, double d3,
                           double e3);

double ozone_damping_rate(double transport_damping, double h, double oh, double ho2,
                          double a2, double a6, double a6b);

double odd_oxygen_chemical_tendency(double oxygen_damping, double oxygen_atom,
                                    double ozone_damping, double ozone);

double effective_ox_damping(double ozone_damping, double oxygen_damping,
                            double transport_damping, double oxygen_atom, double ozone);

// Feedback and adjustment ---------------------------------------------------

double water_vapour_feedback(double radiative_flux_change, double surface_temperature_change);

double tropical_adjustment(double atmospheric_flux, double stratospheric_flux);

// Delta T = -J^-1 Delta H.
std::vector<double> temperature_adjustment(const std::vector<std::vector<double>>& jacobian,
                                           const std::vector<double>& heating_change);

// d[H2O]/dt = 2 alpha d[CH4]_entry/dt (lagged by the age of air).
double methane_water_trend(double fractional_release, double entry_methane_trend);

// Delta H2O_entry = 0.45 ppmv/K Delta T_CPT.
double cold_trap_entry_change(double cold_point_temperature_change);

// Kernel sum and stratospheric-adjusted forcing.
double ozone_kernel_forcing(const std::vector<double>& kernel,
                            const std::vector<double>& ozone_change);

double stratospheric_adjusted_forcing(double shortwave, double longwave_adjusted);

}  // namespace climate

#endif
