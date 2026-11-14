#ifndef CLIMATE_PDE_HPP
#define CLIMATE_PDE_HPP

#include <cstddef>
#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Blatter-Pattyn first-order ice flow on an x-z slice
// ---------------------------------------------------------------------------
// 4 d/dx(mu du/dx) + d/dz(mu du/dz) = rho g ds/dx,
// mu = 0.5 B eps_e^(1/n - 1), eps_e^2 = eps_xx^2 + eps_xz^2.
struct BlatterPattynParameters {
    std::size_t nx = 61;
    std::size_t nz = 31;
    double length = 50000.0;       // m
    double thickness = 1000.0;     // m
    double surface_slope = 0.01;   // ds/dx
    double rate_factor = 1e-16;    // Pa^-n s^-1
    double exponent = 3.0;
    double density = 917.0;        // kg/m^3
    double gravity = 9.81;
    double bed_slip = 0.0;         // m/s, uniform sliding velocity
    std::size_t max_iterations = 600;
    double tolerance = 1e-7;
    double relaxation = 1.0;
    std::size_t sor_sweeps = 400;
};

struct BlatterPattynResult {
    std::vector<double> velocity;  // row-major [i * nz + j], j = 0 at the bed
    std::size_t iterations = 0;
    double residual = 0.0;
    double dx = 0.0;
    double dz = 0.0;
};

BlatterPattynResult solve_blatter_pattyn(const BlatterPattynParameters& parameters);

// ---------------------------------------------------------------------------
// One-dimensional Stefan problem (enthalpy method)
// ---------------------------------------------------------------------------
struct StefanParameters {
    std::size_t cells = 101;
    double length = 10.0;             // m
    double initial_temperature = 265.0;
    double surface_temperature = 275.0;
    double basal_temperature = 265.0;
    double thermal_conductivity = 2.0;
    double density = 1000.0;
    double specific_heat = 2000.0;
    double latent_heat = 3.34e5;
    double freezing_temperature = 273.15;
    double time_step = 3600.0;        // s per reported step
    std::size_t steps = 24;
};

struct StefanResult {
    std::vector<double> temperature;
    std::vector<double> liquid_fraction;
    std::vector<double> enthalpy;
    double front_depth = 0.0;
    double dx = 0.0;
    std::size_t substeps = 0;
};

StefanResult solve_stefan_1d(const StefanParameters& parameters);

// ---------------------------------------------------------------------------
// Phase-field fracture, 1-D
// ---------------------------------------------------------------------------
// dphi/dt = -M[-2(1-phi) W+ + Gc/ell (phi - ell^2 lap(phi))].
struct PhaseFieldParameters {
    std::size_t cells = 101;
    double length = 1.0;
    double time_step = 1e-6;
    std::size_t steps = 100;
    double mobility = 1e-6;
    double fracture_toughness = 1000.0;
    double length_scale = 0.02;
    double tensile_energy = 1e5;
    std::vector<double> initial_phase;  // empty means all ones
};

std::vector<double> solve_phase_field_1d(const PhaseFieldParameters& parameters);

// ---------------------------------------------------------------------------
// Firn densification column
// ---------------------------------------------------------------------------
std::vector<double> evolve_firn_column(std::vector<double> density, double layer_thickness,
                                       double time_step, std::size_t steps, double grain_radius,
                                       double temperature_k,
                                       double activation_energy = 60000.0,
                                       double ice_density = 917.0);

// ---------------------------------------------------------------------------
// Weathering column: advection-reaction for primary minerals and particle age
// ---------------------------------------------------------------------------
// d x_p/dt + v dz x_p = -K tau^sigma x_p,  d tau/dt + v dz tau = 1,
// with v = P0 exp(-h/hp*) directed upward and z = 0 at the bedrock interface.
struct WeatheringParameters {
    std::size_t cells = 101;
    double thickness = 10.0;
    double production_rate = 1e-4;  // m/s (or m/yr, consistent with dt)
    double decay_depth = 2.0;
    double time_step = 100.0;
    std::size_t steps = 100;
    std::vector<double> dissolution;  // K at each node
    double age_exponent = -0.5;
};

struct WeatheringResult {
    std::vector<double> mineral_fraction;
    std::vector<double> particle_age;
    double elapsed = 0.0;
    std::size_t substeps = 0;
};

WeatheringResult solve_weathering_column(const WeatheringParameters& parameters);

}  // namespace climate

#endif
