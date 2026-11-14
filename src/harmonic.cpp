#include "harmonic.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "common.hpp"
#include "linalg.hpp"
#include "stats.hpp"

namespace climate {
namespace {

std::vector<double> design_row(const HarmonicModel& model, double t) {
    const double u = (t - model.time_mean) / model.time_scale;
    std::vector<double> row;
    row.reserve(static_cast<std::size_t>(model.polynomial_order + 1 + 2 * model.harmonic_count));

    double power = 1.0;
    for (int k = 0; k <= model.polynomial_order; ++k) {
        row.push_back(power);
        power *= u;
    }

    const double angular_scale = 2.0 * PI * model.time_scale / model.period;
    for (int k = 1; k <= model.harmonic_count; ++k) {
        const double angle = static_cast<double>(k) * angular_scale * u;
        row.push_back(std::cos(angle));
        row.push_back(std::sin(angle));
    }
    return row;
}

}  // namespace

HarmonicModel fit_harmonic(const std::vector<double>& time, const std::vector<double>& values,
                           int polynomial_order, int harmonic_count, double period,
                           int reject_iterations, double reject_threshold) {
    if (time.size() != values.size())
        throw std::runtime_error("fit_harmonic: time and values must have equal length");
    if (time.empty()) throw std::runtime_error("fit_harmonic: empty input");
    if (polynomial_order < 0) throw std::runtime_error("fit_harmonic: polynomial order must be >= 0");
    if (harmonic_count < 0) throw std::runtime_error("fit_harmonic: harmonic count must be >= 0");
    if (period <= 0.0) throw std::runtime_error("fit_harmonic: period must be positive");

    HarmonicModel model;
    model.polynomial_order = polynomial_order;
    model.harmonic_count = harmonic_count;
    model.period = period;
    model.time_mean = mean_value(time);

    double spread = 0.0;
    for (double t : time) spread = std::max(spread, std::fabs(t - model.time_mean));
    model.time_scale = spread > 0.0 ? spread : 1.0;

    const std::size_t terms = static_cast<std::size_t>(polynomial_order + 1 + 2 * harmonic_count);

    std::vector<bool> keep(time.size(), true);
    std::size_t rejected_total = 0;

    for (int iteration = 0; iteration <= reject_iterations; ++iteration) {
        std::vector<std::vector<double>> a;
        std::vector<double> b;
        a.reserve(time.size());
        b.reserve(time.size());
        for (std::size_t i = 0; i < time.size(); ++i) {
            if (!keep[i] || !std::isfinite(values[i])) continue;
            a.push_back(design_row(model, time[i]));
            b.push_back(values[i]);
        }

        if (a.size() < terms) {
            // Not enough usable samples for the requested model complexity.
            throw std::runtime_error(
                "fit_harmonic: fewer usable observations than model terms (" +
                std::to_string(a.size()) + " < " + std::to_string(terms) + ")");
        }

        std::vector<double> coefficients;
        if (!least_squares(a, b, coefficients))
            throw std::runtime_error("fit_harmonic: design matrix is rank deficient");
        model.coefficients = coefficients;

        double sum_squared = 0.0;
        for (std::size_t i = 0; i < b.size(); ++i) {
            double fitted = 0.0;
            for (std::size_t j = 0; j < terms; ++j) fitted += coefficients[j] * a[i][j];
            const double residual = b[i] - fitted;
            sum_squared += residual * residual;
        }
        const double dof = static_cast<double>(b.size()) - static_cast<double>(terms);
        model.residual_sd = dof > 0.0 ? std::sqrt(sum_squared / dof) : 0.0;
        model.used_points = b.size();

        // Coefficient of determination over the retained points.
        double mean_used = 0.0;
        for (double v : b) mean_used += v;
        mean_used /= static_cast<double>(b.size());
        double total_squared = 0.0;
        for (double v : b) total_squared += (v - mean_used) * (v - mean_used);
        model.r_squared = total_squared > 0.0 ? 1.0 - sum_squared / total_squared : 0.0;

        if (iteration == reject_iterations || model.residual_sd <= 0.0) break;

        const double limit = reject_threshold * model.residual_sd;
        bool rejected_any = false;
        for (std::size_t i = 0; i < time.size(); ++i) {
            if (!keep[i] || !std::isfinite(values[i])) continue;
            const double fitted = harmonic_evaluate(model, time[i]);
            if (std::fabs(values[i] - fitted) > limit) {
                keep[i] = false;
                ++rejected_total;
                rejected_any = true;
            }
        }
        if (!rejected_any) break;
        (void)rejected_total;
    }

    model.rejected_points = rejected_total;
    return model;
}

double harmonic_evaluate(const HarmonicModel& model, double t) {
    const std::vector<double> row = design_row(model, t);
    double value = 0.0;
    for (std::size_t i = 0; i < row.size() && i < model.coefficients.size(); ++i)
        value += model.coefficients[i] * row[i];
    return value;
}

double harmonic_derivative(const HarmonicModel& model, double t) {
    const std::size_t base = static_cast<std::size_t>(model.polynomial_order + 1);
    if (model.coefficients.size() < base) return 0.0;

    const double u = (t - model.time_mean) / model.time_scale;
    double derivative = 0.0;

    // d/du of the polynomial part, then scaled by du/dt = 1/time_scale.
    double power = 1.0;
    for (int k = 1; k <= model.polynomial_order; ++k) {
        derivative += static_cast<double>(k) * model.coefficients[static_cast<std::size_t>(k)] * power;
        power *= u;
    }

    const double angular_scale = 2.0 * PI * model.time_scale / model.period;
    for (int k = 1; k <= model.harmonic_count; ++k) {
        const double angle = static_cast<double>(k) * angular_scale * u;
        const double amp = model.coefficients[base + static_cast<std::size_t>(2 * (k - 1))];
        const double phase = model.coefficients[base + static_cast<std::size_t>(2 * (k - 1) + 1)];
        const double omega = static_cast<double>(k) * angular_scale;
        derivative += omega * (-amp * std::sin(angle) + phase * std::cos(angle));
    }

    return derivative / model.time_scale;
}

double harmonic_trend_derivative(const HarmonicModel& model, double t) {
    const double u = (t - model.time_mean) / model.time_scale;
    double derivative = 0.0;
    double power = 1.0;
    for (int k = 1; k <= model.polynomial_order; ++k) {
        if (static_cast<std::size_t>(k) >= model.coefficients.size()) break;
        derivative += static_cast<double>(k) * model.coefficients[static_cast<std::size_t>(k)] * power;
        power *= u;
    }
    return derivative / model.time_scale;
}

double harmonic_seasonal_value(const HarmonicModel& model, double t) {
    const std::size_t base = static_cast<std::size_t>(model.polynomial_order + 1);
    const double u = (t - model.time_mean) / model.time_scale;
    const double angular_scale = 2.0 * PI * model.time_scale / model.period;
    double value = 0.0;
    for (int k = 1; k <= model.harmonic_count; ++k) {
        const double angle = static_cast<double>(k) * angular_scale * u;
        const double amp = model.coefficients[base + static_cast<std::size_t>(2 * (k - 1))];
        const double phase = model.coefficients[base + static_cast<std::size_t>(2 * (k - 1) + 1)];
        value += amp * std::cos(angle) + phase * std::sin(angle);
    }
    return value;
}

double harmonic_seasonal_amplitude(const HarmonicModel& model) {
    // Seasonal amplitude is the peak-to-trough range of the harmonic part over
    // one full period, sampled densely enough to resolve the highest harmonic.
    const double start = model.time_mean;
    const int samples = std::max(720, 180 * (model.harmonic_count + 1));
    double lowest = harmonic_seasonal_value(model, start);
    double highest = lowest;
    for (int i = 0; i <= samples; ++i) {
        const double t = start + model.period * static_cast<double>(i) / static_cast<double>(samples);
        const double value = harmonic_seasonal_value(model, t);
        lowest = std::min(lowest, value);
        highest = std::max(highest, value);
    }
    return highest - lowest;
}

std::vector<HarmonicComponent> harmonic_components(const HarmonicModel& model) {
    std::vector<HarmonicComponent> components;
    const std::size_t base = static_cast<std::size_t>(model.polynomial_order + 1);
    for (int k = 1; k <= model.harmonic_count; ++k) {
        if (base + static_cast<std::size_t>(2 * k) > model.coefficients.size()) break;
        const double a = model.coefficients[base + static_cast<std::size_t>(2 * (k - 1))];
        const double b = model.coefficients[base + static_cast<std::size_t>(2 * (k - 1) + 1)];
        HarmonicComponent component;
        component.number = k;
        component.amplitude = std::sqrt(a * a + b * b);
        component.phase_radians = std::atan2(b, a);
        component.phase_degrees = rad_to_deg(component.phase_radians);
        components.push_back(component);
    }
    return components;
}

double harmonic_mean_derivative(const HarmonicModel& model, double start, double end, int samples) {
    if (samples < 2) samples = 2;
    if (end == start) return harmonic_derivative(model, start);
    double sum = 0.0;
    for (int i = 0; i < samples; ++i) {
        const double t = start + (end - start) * (static_cast<double>(i) + 0.5) /
                                     static_cast<double>(samples);
        sum += harmonic_derivative(model, t);
    }
    return sum / static_cast<double>(samples);
}

double harmonic_mean_trend_derivative(const HarmonicModel& model, double start, double end,
                                      int samples) {
    if (samples < 2) samples = 2;
    if (end == start) return harmonic_trend_derivative(model, start);
    double sum = 0.0;
    for (int i = 0; i < samples; ++i) {
        const double t = start + (end - start) * (static_cast<double>(i) + 0.5) /
                                     static_cast<double>(samples);
        sum += harmonic_trend_derivative(model, t);
    }
    return sum / static_cast<double>(samples);
}

}  // namespace climate
