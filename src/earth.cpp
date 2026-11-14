#include "earth.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "common.hpp"

namespace climate {
namespace {

constexpr double GRAVITY = 9.80665;       // m s^-2
constexpr double PLANCK = 6.62607015e-34; // J s
constexpr double LIGHT_SPEED = 2.99792458e8;
constexpr double BOLTZMANN = 1.380649e-23;

}  // namespace

double planck_feedback(double surface_temperature_k) {
    if (surface_temperature_k <= 0.0)
        throw std::runtime_error("planck_feedback: temperature must be positive");
    return -4.0 * STEFAN_BOLTZMANN * surface_temperature_k * surface_temperature_k *
           surface_temperature_k;
}

double pattern_corrected_ecs(double historical_feedback, double pattern_correction,
                             double forcing_2xco2) {
    const double total = historical_feedback + pattern_correction;
    if (total == 0.0)
        throw std::runtime_error("pattern_corrected_ecs: zero net feedback");
    return -forcing_2xco2 / total;
}

namespace {

TwoLayerState two_layer_derivative(const TwoLayerModel& model, double forcing,
                                   const TwoLayerState& state) {
    const double exchange = model.uptake_efficiency * model.efficacy *
                            (state.surface - state.deep);
    TwoLayerState derivative;
    // The notes write C dT/dt = F - lambda T with lambda positive; this model
    // stores the IPCC-style negative lambda, so the term enters as +lambda T.
    derivative.surface = (forcing + model.feedback * state.surface - exchange) /
                         model.heat_capacity;
    derivative.deep = model.uptake_efficiency * (state.surface - state.deep) /
                      model.deep_heat_capacity;
    return derivative;
}

TwoLayerState add_scaled(const TwoLayerState& state, const TwoLayerState& derivative,
                         double scale) {
    TwoLayerState result;
    result.surface = state.surface + scale * derivative.surface;
    result.deep = state.deep + scale * derivative.deep;
    return result;
}

}  // namespace

double two_layer_equilibrium(const TwoLayerModel& model, double forcing) {
    if (model.feedback == 0.0)
        throw std::runtime_error("two_layer_equilibrium: zero feedback parameter");
    return -forcing / model.feedback;
}

TwoLayerState integrate_two_layer(const TwoLayerModel& model, double forcing, double years,
                                  double dt_years, TwoLayerState initial) {
    if (model.heat_capacity <= 0.0 || model.deep_heat_capacity <= 0.0)
        throw std::runtime_error("integrate_two_layer: heat capacities must be positive");
    if (dt_years <= 0.0) throw std::runtime_error("integrate_two_layer: dt must be positive");
    const long steps = static_cast<long>(std::ceil(years / dt_years));
    const double dt = years / static_cast<double>(std::max<long>(steps, 1));

    TwoLayerState state = initial;
    for (long step = 0; step < steps; ++step) {
        const TwoLayerState k1 = two_layer_derivative(model, forcing, state);
        const TwoLayerState k2 =
            two_layer_derivative(model, forcing, add_scaled(state, k1, dt / 2.0));
        const TwoLayerState k3 =
            two_layer_derivative(model, forcing, add_scaled(state, k2, dt / 2.0));
        const TwoLayerState k4 = two_layer_derivative(model, forcing, add_scaled(state, k3, dt));
        state.surface += dt / 6.0 * (k1.surface + 2.0 * k2.surface + 2.0 * k3.surface + k4.surface);
        state.deep += dt / 6.0 * (k1.deep + 2.0 * k2.deep + 2.0 * k3.deep + k4.deep);
    }
    return state;
}

double ocean_heat_uptake(const TwoLayerModel& model, const TwoLayerState& state) {
    return model.uptake_efficiency * (state.surface - state.deep);
}

double ocean_heat_uptake_efficiency(const TwoLayerModel& model, const TwoLayerState& state) {
    if (state.surface == 0.0)
        throw std::runtime_error("ocean_heat_uptake_efficiency: zero surface warming");
    return ocean_heat_uptake(model, state) / state.surface;
}

double transient_climate_response(const TwoLayerModel& model, double forcing_2xco2) {
    // With the negative feedback convention the damping is -lambda + eps gamma.
    const double denominator = -model.feedback + model.efficacy * model.uptake_efficiency;
    if (denominator == 0.0)
        throw std::runtime_error("transient_climate_response: zero denominator");
    return forcing_2xco2 / denominator;
}

double equilibration_ratio(const TwoLayerModel& model, double forcing_2xco2) {
    const double tcr = transient_climate_response(model, forcing_2xco2);
    if (tcr == 0.0) throw std::runtime_error("equilibration_ratio: zero TCR");
    return two_layer_equilibrium(model, forcing_2xco2) / tcr;
}

double land_carbon_change(const CarbonFeedback& feedback, double delta_c_atmos,
                          double delta_temperature) {
    return feedback.beta_land * delta_c_atmos + feedback.gamma_land * delta_temperature;
}

double ocean_carbon_change(const CarbonFeedback& feedback, double delta_c_atmos,
                           double delta_temperature) {
    return feedback.beta_ocean * delta_c_atmos + feedback.gamma_ocean * delta_temperature;
}

double airborne_fraction(double atmospheric_growth_ppm, double emissions_pgc) {
    if (emissions_pgc == 0.0)
        throw std::runtime_error("airborne_fraction: emissions must be non-zero");
    return atmospheric_growth_ppm * CARBON_PPM_TO_PGC / emissions_pgc;
}

double permafrost_decay_rate(double k_reference, double q10, double temperature_c,
                             double reference_temperature_c) {
    if (k_reference < 0.0) throw std::runtime_error("permafrost_decay_rate: k must be >= 0");
    if (q10 <= 0.0) throw std::runtime_error("permafrost_decay_rate: Q10 must be positive");
    return k_reference *
           std::pow(q10, (temperature_c - reference_temperature_c) / 10.0);
}

double permafrost_remaining_fraction(double k_reference, double q10, double temperature_c,
                                     double years, double reference_temperature_c) {
    if (years < 0.0) throw std::runtime_error("permafrost_remaining_fraction: time must be >= 0");
    const double rate = permafrost_decay_rate(k_reference, q10, temperature_c,
                                              reference_temperature_c);
    return std::exp(-rate * years);
}

double planck_radiance(double wavelength_m, double temperature_k) {
    if (!(wavelength_m > 0.0)) throw std::runtime_error("planck_radiance: wavelength must be positive");
    if (!(temperature_k > 0.0)) throw std::runtime_error("planck_radiance: temperature must be positive");
    const double lambda5 = std::pow(wavelength_m, 5.0);
    const double numerator = 2.0 * PLANCK * LIGHT_SPEED * LIGHT_SPEED;
    const double exponent = PLANCK * LIGHT_SPEED / (wavelength_m * BOLTZMANN * temperature_k);
    // expm1 keeps precision for long wavelengths / low temperatures.
    return numerator / lambda5 / std::expm1(exponent);
}

double brightness_temperature(double wavelength_m, double radiance) {
    if (!(wavelength_m > 0.0)) throw std::runtime_error("brightness_temperature: wavelength must be positive");
    if (!(radiance > 0.0)) throw std::runtime_error("brightness_temperature: radiance must be positive");
    const double numerator = 2.0 * PLANCK * LIGHT_SPEED * LIGHT_SPEED;
    const double inside = 1.0 + numerator / (std::pow(wavelength_m, 5.0) * radiance);
    return PLANCK * LIGHT_SPEED / (wavelength_m * BOLTZMANN * std::log(inside));
}

double optical_depth(double absorption_coefficient, double density, double path_length) {
    if (absorption_coefficient < 0.0 || density < 0.0 || path_length < 0.0)
        throw std::runtime_error("optical_depth: inputs must be non-negative");
    return absorption_coefficient * density * path_length;
}

double transmittance(double tau) {
    if (tau < 0.0) throw std::runtime_error("transmittance: optical depth must be >= 0");
    return std::exp(-tau);
}

double weighting_function(double absorption_coefficient, double specific_humidity,
                          double tau, double gravity) {
    if (tau < 0.0) throw std::runtime_error("weighting_function: optical depth must be >= 0");
    if (gravity == 0.0) throw std::runtime_error("weighting_function: gravity must be non-zero");
    return absorption_coefficient * specific_humidity * std::exp(-tau) / gravity;
}

double contribution_function(double planck_radiance_value, double weighting) {
    return planck_radiance_value * weighting;
}

double layer_radiance(double surface_planck, double layer_planck, double tau) {
    if (tau < 0.0) throw std::runtime_error("layer_radiance: optical depth must be >= 0");
    const double transmission = std::exp(-tau);
    return surface_planck * transmission + layer_planck * (1.0 - transmission);
}

TwoStreamCoefficients two_stream_eddington(double single_scattering_albedo, double asymmetry,
                                           double mu0) {
    const double w = single_scattering_albedo;
    if (w < 0.0 || w > 1.0)
        throw std::runtime_error("two_stream_eddington: single-scattering albedo must be in [0,1]");
    if (asymmetry < -1.0 || asymmetry > 1.0)
        throw std::runtime_error("two_stream_eddington: asymmetry must be in [-1,1]");
    if (mu0 <= 0.0 || mu0 > 1.0)
        throw std::runtime_error("two_stream_eddington: mu0 must be in (0,1]");
    TwoStreamCoefficients coefficients;
    coefficients.g1 = (7.0 - w * (4.0 + 3.0 * asymmetry)) / 4.0;
    coefficients.g2 = -(1.0 - w * (4.0 - 3.0 * asymmetry)) / 4.0;
    coefficients.g3 = (2.0 - 3.0 * mu0 * asymmetry) / 4.0;
    coefficients.g4 = 1.0 - coefficients.g3;
    return coefficients;
}

namespace {

// Dense Gaussian elimination with partial pivoting for the two-stream BVP.
std::vector<double> solve_dense(std::vector<std::vector<double>> matrix,
                                std::vector<double> rhs) {
    const std::size_t n = rhs.size();
    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        double best = std::fabs(matrix[column][column]);
        for (std::size_t row = column + 1; row < n; ++row) {
            const double candidate = std::fabs(matrix[row][column]);
            if (candidate > best) {
                best = candidate;
                pivot = row;
            }
        }
        if (best < 1e-14)
            throw std::runtime_error("two_stream_layer: singular system; check optical depth");
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(rhs[pivot], rhs[column]);
        }
        const double diagonal = matrix[column][column];
        for (std::size_t row = column + 1; row < n; ++row) {
            const double factor = matrix[row][column] / diagonal;
            if (factor == 0.0) continue;
            for (std::size_t k = column; k < n; ++k) matrix[row][k] -= factor * matrix[column][k];
            rhs[row] -= factor * rhs[column];
        }
    }
    std::vector<double> solution(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double sum = rhs[i];
        for (std::size_t j = i + 1; j < n; ++j) sum -= matrix[i][j] * solution[j];
        solution[i] = sum / matrix[i][i];
    }
    return solution;
}

}  // namespace

TwoStreamResult two_stream_layer(double tau_star, double single_scattering_albedo,
                                 double asymmetry, double solar_zenith_deg,
                                 double surface_albedo, int levels) {
    if (tau_star < 0.0) throw std::runtime_error("two_stream_layer: optical depth must be >= 0");
    if (surface_albedo < 0.0 || surface_albedo > 1.0)
        throw std::runtime_error("two_stream_layer: surface albedo must be in [0,1]");
    if (levels < 4) throw std::runtime_error("two_stream_layer: need at least 4 levels");

    const double mu0 = std::cos(deg_to_rad(solar_zenith_deg));
    if (mu0 <= 0.0)
        throw std::runtime_error("two_stream_layer: the sun is below the horizon");

    TwoStreamResult result;
    result.direct_transmittance = std::exp(-tau_star / mu0);

    // No scattering: pure Beer-Lambert extinction of the direct beam.
    if (single_scattering_albedo <= 0.0 || tau_star == 0.0) {
        result.total_transmittance = (1.0 - surface_albedo) * result.direct_transmittance;
        result.diffuse_transmittance = 0.0;
        result.reflectance = 0.0;
        result.absorptance = 1.0 - result.total_transmittance;
        return result;
    }

    const TwoStreamCoefficients coefficients =
        two_stream_eddington(single_scattering_albedo, asymmetry, mu0);

    const std::size_t n = static_cast<std::size_t>(levels);
    const double h = tau_star / static_cast<double>(levels);
    const double source_scale = single_scattering_albedo;  // F0 = 1, so mu0 F0 = mu0
    const auto source = [&](double tau) {
        return std::exp(-tau / mu0);
    };

    // Unknowns: F+[0..N] then F-[0..N], 2(N+1) in total.
    const std::size_t unknowns = 2 * (n + 1);
    std::vector<std::vector<double>> matrix(unknowns, std::vector<double>(unknowns, 0.0));
    std::vector<double> rhs(unknowns, 0.0);
    const auto plus = [&](std::size_t i) { return i; };
    const auto minus = [&](std::size_t i) { return n + 1 + i; };

    // Top boundary: no diffuse downward radiation; the outgoing upward flux is
    // the reflected beam, solved from a one-sided difference at tau = 0.
    matrix[minus(0)][minus(0)] = 1.0;
    rhs[minus(0)] = 0.0;
    {
        const double scale = source_scale * source(0.0);
        auto& row = matrix[plus(0)];
        row[plus(1)] += 1.0 / h;
        row[plus(0)] -= 1.0 / h + coefficients.g1;
        row[minus(0)] += coefficients.g2;
        rhs[plus(0)] = -coefficients.g3 * scale;
    }


    // Interior: central differences for F+, upwind differences for F-.
    for (std::size_t i = 1; i < n; ++i) {
        const double tau = static_cast<double>(i) * h;
        const double scale = source_scale * source(tau);

        auto& row_plus = matrix[plus(i)];
        row_plus[plus(i + 1)] += 1.0 / (2.0 * h);
        row_plus[plus(i - 1)] -= 1.0 / (2.0 * h);
        row_plus[plus(i)] -= coefficients.g1;
        row_plus[minus(i)] += coefficients.g2;
        rhs[plus(i)] = -coefficients.g3 * scale;

        auto& row_minus = matrix[minus(i)];
        row_minus[minus(i + 1)] += 1.0 / (2.0 * h);
        row_minus[minus(i - 1)] -= 1.0 / (2.0 * h);
        row_minus[minus(i)] += coefficients.g1;
        row_minus[plus(i)] -= coefficients.g2;
        rhs[minus(i)] = coefficients.g4 * scale;
    }

    // Bottom of the layer, same upwind form for F-.
    {
        const double scale = source_scale * source(tau_star);
        auto& row = matrix[minus(n)];
        row[minus(n)] += 1.0 / h + coefficients.g1;
        row[minus(n - 1)] -= 1.0 / h;
        row[plus(n)] -= coefficients.g2;
        rhs[minus(n)] = coefficients.g4 * scale;
    }

    // Surface reflection: F+(tau*) = albedo (F-(tau*) + direct beam).
    {
        auto& row = matrix[plus(n)];
        row[plus(n)] += 1.0;
        row[minus(n)] -= surface_albedo;
        rhs[plus(n)] = surface_albedo * mu0 * result.direct_transmittance;
    }

    const std::vector<double> solution = solve_dense(matrix, rhs);
    std::vector<double> f_plus(n + 1), f_minus(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        f_plus[i] = solution[plus(i)];
        f_minus[i] = solution[minus(i)];
    }

    const double incident = mu0;  // F0 = 1
    result.reflectance = f_plus[0] / incident;
    result.diffuse_transmittance = f_minus[n] / incident;
    result.total_transmittance =
        (1.0 - surface_albedo) * (f_minus[n] + mu0 * result.direct_transmittance) / incident;
    result.absorptance = 1.0 - result.reflectance - result.total_transmittance;
    return result;
}

double shallow_water_speed(double gravity, double depth) {
    if (gravity <= 0.0) throw std::runtime_error("shallow_water_speed: gravity must be positive");
    if (depth < 0.0) throw std::runtime_error("shallow_water_speed: depth must be non-negative");
    return std::sqrt(gravity * depth);
}

EquatorialWave kelvin_wave(double zonal_wavenumber, double gravity, double depth) {
    EquatorialWave wave;
    const double speed = shallow_water_speed(gravity, depth);
    wave.phase_speed = speed;
    wave.frequency = speed * zonal_wavenumber;
    return wave;
}

std::vector<double> equatorial_wave_frequencies(double k, double beta, double c,
                                                int meridional_mode) {
    if (c <= 0.0) throw std::runtime_error("equatorial_wave_frequencies: c must be positive");
    if (meridional_mode < 0)
        throw std::runtime_error("equatorial_wave_frequencies: mode must be >= 0");

    // omega^3/c^2 + 0 omega^2 - (k^2 + (2n+1) beta/c) omega - beta k = 0
    const double a = 1.0 / (c * c);
    const double b = 0.0;
    const double cc = -(k * k + (2.0 * meridional_mode + 1.0) * beta / c);
    const double d = -beta * k;

    // Depressed cubic t^3 + p t + q = 0 with omega = t - b/(3a).
    const double shift = b / (3.0 * a);
    const double p = (3.0 * a * cc - b * b) / (3.0 * a * a);
    const double q = (2.0 * b * b * b - 9.0 * a * b * cc + 27.0 * a * a * d) / (27.0 * a * a * a);

    std::vector<double> roots;
    const double discriminant = 0.25 * q * q + p * p * p / 27.0;
    if (discriminant > 0.0) {
        const double sqrt_disc = std::sqrt(discriminant);
        const double u = std::cbrt(-0.5 * q + sqrt_disc);
        const double v = std::cbrt(-0.5 * q - sqrt_disc);
        roots.push_back(u + v - shift);
    } else {
        const double radius = 2.0 * std::sqrt(std::max(0.0, -p / 3.0));
        double argument = 0.0;
        if (radius > 0.0) {
            argument = 3.0 * q / (p * radius);
            argument = clamp_value(argument, -1.0, 1.0);
        }
        const double angle = std::acos(argument) / 3.0;
        for (int j = 0; j < 3; ++j) {
            roots.push_back(radius * std::cos(angle - 2.0 * PI * static_cast<double>(j) / 3.0) -
                            shift);
        }
    }

    // Keep physically meaningful (positive) frequencies, sorted ascending.
    std::vector<double> physical;
    for (const double root : roots) {
        if (root > 1e-14) physical.push_back(root);
    }
    std::sort(physical.begin(), physical.end());
    return physical;
}

std::vector<EquatorialWave> equatorial_wave_modes(double k, double beta, double c,
                                                  int meridional_mode) {
    std::vector<EquatorialWave> modes;
    for (const double frequency : equatorial_wave_frequencies(k, beta, c, meridional_mode)) {
        EquatorialWave wave;
        wave.frequency = frequency;
        wave.phase_speed = k != 0.0 ? frequency / k : 0.0;
        modes.push_back(wave);
    }
    return modes;
}

double refractive_index_spectrum(double cn2, double wavenumber) {
    if (cn2 < 0.0) throw std::runtime_error("refractive_index_spectrum: Cn2 must be >= 0");
    if (!(wavenumber > 0.0))
        throw std::runtime_error("refractive_index_spectrum: wavenumber must be positive");
    return 0.033 * cn2 * std::pow(wavenumber, -11.0 / 3.0);
}

double dissipation_rate_from_spectrum(double frequency, double velocity_spectrum,
                                      double wind_speed, double kolmogorov_constant) {
    if (!(frequency > 0.0)) throw std::runtime_error("dissipation_rate: frequency must be positive");
    if (velocity_spectrum < 0.0)
        throw std::runtime_error("dissipation_rate: spectrum must be non-negative");
    if (!(wind_speed > 0.0)) throw std::runtime_error("dissipation_rate: wind speed must be positive");
    if (!(kolmogorov_constant > 0.0))
        throw std::runtime_error("dissipation_rate: Kolmogorov constant must be positive");
    return 2.0 * PI * std::pow(frequency, 2.5) * std::pow(velocity_spectrum, 1.5) /
           (wind_speed * std::pow(kolmogorov_constant, 1.5));
}

double kolmogorov_microscale(double kinematic_viscosity, double dissipation_rate, double a2) {
    if (kinematic_viscosity < 0.0)
        throw std::runtime_error("kolmogorov_microscale: viscosity must be >= 0");
    if (!(dissipation_rate > 0.0))
        throw std::runtime_error("kolmogorov_microscale: dissipation rate must be positive");
    return a2 * std::pow(kinematic_viscosity * kinematic_viscosity * kinematic_viscosity /
                             dissipation_rate,
                         0.25);
}

double kolmogorov_spectrum(double wavenumber, double dissipation_rate, double outer_scale,
                           double kolmogorov_constant) {
    if (!(wavenumber > 0.0)) throw std::runtime_error("kolmogorov_spectrum: wavenumber must be positive");
    if (!(dissipation_rate > 0.0))
        throw std::runtime_error("kolmogorov_spectrum: dissipation rate must be positive");
    if (!(outer_scale > 0.0))
        throw std::runtime_error("kolmogorov_spectrum: outer scale must be positive");
    const double kL = wavenumber * outer_scale;
    return kolmogorov_constant * std::pow(dissipation_rate, 2.0 / 3.0) *
           std::pow(outer_scale, 5.0 / 3.0) * std::pow(kL, 4.0) /
           std::pow(1.0 + kL * kL, 17.0 / 6.0);
}

double stable_boundary_layer_ratio(double stability_parameter, double c_s) {
    if (!(stability_parameter > 0.0))
        throw std::runtime_error("stable_boundary_layer_ratio: mu must be positive");
    if (!(c_s > 0.0)) throw std::runtime_error("stable_boundary_layer_ratio: c_s must be positive");
    return 1.0 / (c_s * std::sqrt(stability_parameter));
}

}  // namespace climate
