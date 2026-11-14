#include "pde.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double clamp(double value, double low, double high) {
    return std::min(std::max(value, low), high);
}

}  // namespace

BlatterPattynResult solve_blatter_pattyn(const BlatterPattynParameters& parameters) {
    require(parameters.nx >= 4 && parameters.nz >= 4, "blatter_pattyn: grid too small");
    require(parameters.length > 0.0, "blatter_pattyn: length must be positive");
    require(parameters.thickness > 0.0, "blatter_pattyn: thickness must be positive");
    require(parameters.rate_factor > 0.0, "blatter_pattyn: rate factor must be positive");

    const std::size_t nx = parameters.nx;
    const std::size_t nz = parameters.nz;
    const double dx = parameters.length / static_cast<double>(nx - 1);
    const double dz = parameters.thickness / static_cast<double>(nz - 1);
    const double n = parameters.exponent;
    const double hardness = std::pow(parameters.rate_factor, -1.0 / n);
    const double driving = parameters.density * parameters.gravity * parameters.surface_slope;

    // Start from the shallow-ice velocity scale so Picard iteration does not
    // stall on the low-strain / high-viscosity branch.
    const double scale = parameters.surface_slope > 0.0
                             ? 2.0 * parameters.rate_factor / (n + 1.0) *
                                   std::pow(driving, n) *
                                   std::pow(parameters.thickness, n + 1.0)
                             : 0.0;
    std::vector<double> u(nx * nz, parameters.bed_slip);
    std::vector<double> mu(nx * nz, 1.0);
    std::vector<double> next(nx * nz, parameters.bed_slip);

    auto at = [nz](std::size_t i, std::size_t j) { return i * nz + j; };
    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < nz; ++j) {
            const double fraction = static_cast<double>(j) / static_cast<double>(nz - 1);
            u[at(i, j)] = parameters.bed_slip + scale * fraction * fraction;
        }
    }

    double residual = 0.0;
    std::size_t iteration = 0;
    for (; iteration < parameters.max_iterations; ++iteration) {
        // Viscosity on the current velocity field.
        for (std::size_t i = 0; i < nx; ++i) {
            for (std::size_t j = 0; j < nz; ++j) {
                const std::size_t iw = i == 0 ? 0 : i - 1;
                const std::size_t ie = i + 1 < nx ? i + 1 : i;
                const std::size_t js = j == 0 ? 0 : j - 1;
                const std::size_t jn = j + 1 < nz ? j + 1 : j;
                const double du_dx = (u[at(ie, j)] - u[at(iw, j)]) /
                                     ((ie - iw) * dx);
                const double du_dz = (u[at(i, jn)] - u[at(i, js)]) /
                                     ((jn - js) * dz);
                const double eps_e = std::sqrt(du_dx * du_dx + 0.25 * du_dz * du_dz);
                const double safe = std::max(eps_e, 1e-12);
                mu[at(i, j)] = 0.5 * hardness * std::pow(safe, 1.0 / n - 1.0);
            }
        }

        // The surface viscosity is taken from the level below: at the
        // stress-free surface both stress and strain rate vanish together,
        // so the ratio would otherwise be dominated by discretisation noise.
        for (std::size_t i = 0; i < nx; ++i) mu[at(i, nz - 1)] = mu[at(i, nz - 2)];

        // Gauss-Seidel/SOR sweeps on the linearised operator.
        next = u;
        for (std::size_t sweep = 0; sweep < parameters.sor_sweeps; ++sweep) {
            for (std::size_t i = 1; i + 1 < nx; ++i) {
                for (std::size_t j = 1; j < nz; ++j) {
                    const double centre_mu = mu[at(i, j)];
                    double east = 0.0, west = 0.0, north = 0.0, south = 0.0;
                    double diagonal = 0.0;
                    double rhs = driving;

                    if (i + 1 < nx) {
                        const double face = 0.5 * (centre_mu + mu[at(i + 1, j)]);
                        east = 4.0 * face / (dx * dx);
                        diagonal += east;
                    } else {
                        rhs += 0.0;  // zero-gradient outflow: u[nx-1] = u[nx-2]
                    }
                    {
                        const double face = 0.5 * (centre_mu + mu[at(i - 1, j)]);
                        west = 4.0 * face / (dx * dx);
                        diagonal += west;
                    }
                    if (j + 1 < nz) {
                        const double face = 0.5 * (centre_mu + mu[at(i, j + 1)]);
                        north = face / (dz * dz);
                        diagonal += north;
                    }
                    // At the stress-free surface the ghost value equals the
                    // surface value, so there is no extra north contribution.
                    {
                        const double face = 0.5 * (centre_mu + mu[at(i, j - 1)]);
                        south = face / (dz * dz);
                        diagonal += south;
                    }

                    double value = rhs;
                    if (i + 1 < nx) value += east * next[at(i + 1, j)];
                    else value += east * next[at(i - 1, j)];
                    value += west * next[at(i - 1, j)];
                    if (j + 1 < nz) value += north * next[at(i, j + 1)];
                    value += south * next[at(i, j - 1)];

                    const double target = value / diagonal;
                    next[at(i, j)] = (1.0 - parameters.relaxation) * next[at(i, j)] +
                                     parameters.relaxation * target;
                }
            }
            // Boundary conditions: inflow and bed fixed, outflow zero gradient.
            for (std::size_t j = 0; j < nz; ++j) {
                next[at(0, j)] = 0.0;
                next[at(nx - 1, j)] = next[at(nx - 2, j)];
            }
            for (std::size_t i = 0; i < nx; ++i) next[at(i, 0)] = parameters.bed_slip;
        }

        residual = 0.0;
        double magnitude = 1.0;
        for (std::size_t k = 0; k < u.size(); ++k) {
            residual = std::max(residual, std::fabs(next[k] - u[k]));
            magnitude = std::max(magnitude, std::fabs(next[k]));
        }
        u.swap(next);
        if (residual < parameters.tolerance * magnitude) {
            ++iteration;
            break;
        }
    }

    BlatterPattynResult result;
    result.velocity = std::move(u);
    result.iterations = iteration;
    result.residual = residual;
    result.dx = dx;
    result.dz = dz;
    return result;
}

StefanResult solve_stefan_1d(const StefanParameters& parameters) {
    require(parameters.cells >= 3, "stefan: at least three cells");
    require(parameters.length > 0.0, "stefan: length must be positive");
    require(parameters.density > 0.0 && parameters.specific_heat > 0.0,
            "stefan: density and heat capacity must be positive");
    require(parameters.thermal_conductivity > 0.0, "stefan: conductivity must be positive");

    const std::size_t n = parameters.cells;
    const double dx = parameters.length / static_cast<double>(n - 1);
    const double sensible = parameters.density * parameters.specific_heat;
    const double latent = parameters.density * parameters.latent_heat;
    const double freezing = parameters.freezing_temperature;

    auto enthalpy_of = [&](double temperature) {
        const double offset = sensible * (temperature - freezing);
        return temperature <= freezing ? offset : latent + offset;
    };
    auto temperature_of = [&](double enthalpy) {
        if (enthalpy <= 0.0) return freezing + enthalpy / sensible;
        if (enthalpy <= latent) return freezing;
        return freezing + (enthalpy - latent) / sensible;
    };
    auto liquid_of = [&](double enthalpy) {
        if (enthalpy <= 0.0) return 0.0;
        if (enthalpy >= latent) return 1.0;
        return enthalpy / latent;
    };

    std::vector<double> enthalpy(n), temperature(n), liquid(n);
    for (std::size_t i = 0; i < n; ++i) {
        temperature[i] = parameters.initial_temperature;
        enthalpy[i] = enthalpy_of(parameters.initial_temperature);
        liquid[i] = liquid_of(enthalpy[i]);
    }

    const double alpha = parameters.thermal_conductivity / sensible;
    const double stable_dt = 0.4 * dx * dx / alpha;
    const std::size_t substeps =
        static_cast<std::size_t>(std::ceil(parameters.time_step / stable_dt));
    const double dt = parameters.time_step / static_cast<double>(substeps);

    for (std::size_t step = 0; step < parameters.steps; ++step) {
        for (std::size_t sub = 0; sub < substeps; ++sub) {
            for (std::size_t i = 1; i + 1 < n; ++i) temperature[i] = temperature_of(enthalpy[i]);
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const double flux = parameters.thermal_conductivity *
                                    (temperature[i + 1] - 2.0 * temperature[i] +
                                     temperature[i - 1]) /
                                    (dx * dx);
                enthalpy[i] += dt * flux;
            }
            enthalpy[0] = enthalpy_of(parameters.surface_temperature);
            enthalpy[n - 1] = enthalpy_of(parameters.basal_temperature);
            temperature[0] = parameters.surface_temperature;
            temperature[n - 1] = parameters.basal_temperature;
            for (std::size_t i = 1; i + 1 < n; ++i) temperature[i] = temperature_of(enthalpy[i]);
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        liquid[i] = liquid_of(enthalpy[i]);
        temperature[i] = temperature_of(enthalpy[i]);
    }

    StefanResult result;
    result.temperature = std::move(temperature);
    result.liquid_fraction = std::move(liquid);
    result.enthalpy = std::move(enthalpy);
    result.dx = dx;
    result.substeps = substeps;
    double depth = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        if (result.liquid_fraction[i] > 0.5) depth = static_cast<double>(i) * dx;
    result.front_depth = depth;
    return result;
}

std::vector<double> solve_phase_field_1d(const PhaseFieldParameters& parameters) {
    require(parameters.cells >= 3, "phase_field: at least three cells");
    require(parameters.length > 0.0, "phase_field: length must be positive");
    require(parameters.length_scale > 0.0, "phase_field: length scale must be positive");

    const std::size_t n = parameters.cells;
    const double dx = parameters.length / static_cast<double>(n - 1);
    std::vector<double> phi(n, 1.0);
    if (!parameters.initial_phase.empty()) {
        require(parameters.initial_phase.size() == n, "phase_field: initial size mismatch");
        phi = parameters.initial_phase;
    }

    std::vector<double> next(n, 0.0);
    for (std::size_t step = 0; step < parameters.steps; ++step) {
        next = phi;
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const double laplacian = (phi[i + 1] - 2.0 * phi[i] + phi[i - 1]) / (dx * dx);
            const double variational =
                -2.0 * (1.0 - phi[i]) * parameters.tensile_energy +
                parameters.fracture_toughness / parameters.length_scale *
                    (phi[i] - parameters.length_scale * parameters.length_scale * laplacian);
            const double updated = phi[i] - parameters.time_step * parameters.mobility * variational;
            next[i] = clamp(updated, 0.0, 1.0);
        }
        next[0] = next[1];
        next[n - 1] = next[n - 2];
        phi.swap(next);
    }
    return phi;
}

std::vector<double> evolve_firn_column(std::vector<double> density, double layer_thickness,
                                       double time_step, std::size_t steps, double grain_radius,
                                       double temperature_k, double activation_energy,
                                       double ice_density) {
    require(density.size() >= 2, "firn: at least two layers");
    require(layer_thickness > 0.0, "firn: layer thickness must be positive");
    require(grain_radius > 0.0, "firn: grain radius must be positive");
    require(temperature_k > 0.0, "firn: temperature must be positive");
    const double gravity = 9.80665;
    const double gas_constant = 8.314462618;
    const double rate = std::exp(-activation_energy / (gas_constant * temperature_k));

    const std::size_t n = density.size();
    std::vector<double> overburden(n, 0.0);
    for (std::size_t step = 0; step < steps; ++step) {
        // density[0] is the surface layer, density[n-1] the deepest.
        overburden[n - 1] = 0.0;
        for (std::size_t i = n - 1; i > 0; --i)
            overburden[i - 1] = overburden[i] + density[i] * gravity * layer_thickness;
        for (std::size_t i = 0; i < n; ++i) {
            const double kc = density[i] < 550.0 ? 9.2e-9 : 3.7e-9;
            const double rate_density =
                kc * overburden[i] / (grain_radius * grain_radius) * (ice_density - density[i]) *
                rate;
            density[i] = std::min(ice_density, density[i] + time_step * rate_density);
        }
    }
    return density;
}

WeatheringResult solve_weathering_column(const WeatheringParameters& parameters) {
    require(parameters.cells >= 3, "weathering: at least three cells");
    require(parameters.thickness > 0.0, "weathering: thickness must be positive");
    const std::size_t n = parameters.cells;
    const double dz = parameters.thickness / static_cast<double>(n - 1);
    const double velocity = parameters.production_rate;
    require(velocity >= 0.0, "weathering: production rate must be >= 0");
    if (!parameters.dissolution.empty()) {
        require(parameters.dissolution.size() == n, "weathering: dissolution size mismatch");
    }

    std::vector<double> mineral(n, 1.0);
    std::vector<double> age(n, 0.0);
    // The advective CFL limit; fall back to one substep if v = 0.
    const double cfl_dt = velocity > 0.0 ? 0.5 * dz / velocity : parameters.time_step;
    const std::size_t substeps = static_cast<std::size_t>(
        std::ceil(parameters.time_step / std::max(cfl_dt, 1e-30)));
    const double dt = parameters.time_step / static_cast<double>(substeps);

    std::vector<double> next_mineral(n, 0.0);
    std::vector<double> next_age(n, 0.0);
    for (std::size_t step = 0; step < parameters.steps; ++step) {
        for (std::size_t sub = 0; sub < substeps; ++sub) {
            next_mineral = mineral;
            next_age = age;
            for (std::size_t i = 1; i < n; ++i) {
                const double advection = velocity * (mineral[i] - mineral[i - 1]) / dz;
                const double k = parameters.dissolution.empty() ? 0.0 : parameters.dissolution[i];
                const double age_value = std::max(age[i], 0.0);
                const double normalised_age =
                    parameters.reference_age > 0.0 ? age_value / parameters.reference_age : age_value;
                const double reactivity =
                    normalised_age > 0.0 ? std::pow(normalised_age, parameters.age_exponent) : 0.0;
                const double reaction = k * reactivity * mineral[i];
                next_mineral[i] = mineral[i] - dt * advection - dt * reaction;
                const double age_advection = velocity * (age[i] - age[i - 1]) / dz;
                next_age[i] = std::max(0.0, age[i] + dt * (1.0 - age_advection));
                if (next_mineral[i] < 0.0) next_mineral[i] = 0.0;
            }
            next_mineral[0] = 1.0;
            next_age[0] = 0.0;
            mineral.swap(next_mineral);
            age.swap(next_age);
        }
    }

    WeatheringResult result;
    result.mineral_fraction = std::move(mineral);
    result.particle_age = std::move(age);
    result.elapsed = parameters.time_step * static_cast<double>(parameters.steps);
    result.substeps = substeps;
    return result;
}

}  // namespace climate
