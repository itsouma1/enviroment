#include "stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "loess.hpp"
#include "stats.hpp"

namespace climate {
namespace {

int make_odd(int value) {
    if (value < 1) return 1;
    return (value % 2 == 0) ? value + 1 : value;
}

}  // namespace

STLResult stl_decompose(const std::vector<double>& values, int period, const STLParams& params) {
    const std::size_t n = values.size();
    if (period < 2) throw std::runtime_error("stl_decompose: period must be at least 2");
    if (n < static_cast<std::size_t>(2 * period))
        throw std::runtime_error("stl_decompose: need at least two full seasonal cycles");

    for (double value : values) {
        if (!std::isfinite(value))
            throw std::runtime_error(
                "stl_decompose: series contains missing values; fill them first");
    }

    STLResult result;
    result.seasonal.assign(n, 0.0);
    result.trend.assign(n, 0.0);
    result.remainder.assign(n, 0.0);
    result.robustness.assign(n, 1.0);

    const int seasonal_window = std::max(make_odd(params.seasonal_window), params.seasonal_degree + 2);
    const int trend_window = std::max(make_odd(params.trend_window), params.trend_degree + 2);
    const int low_pass_window = std::max(make_odd(params.low_pass_window), 3);
    const int inner_iterations = std::max(1, params.inner_iterations);
    const int outer_iterations = std::max(0, params.outer_iterations);

    std::vector<double> seasonal(n, 0.0);
    std::vector<double> trend(n, 0.0);
    std::vector<double> robustness(n, 1.0);
    std::vector<double> detrended(n, 0.0);
    std::vector<double> cycle(n, 0.0);
    std::vector<double> deseasonalized(n, 0.0);

    for (int outer = 0; outer <= outer_iterations; ++outer) {
        for (int inner = 0; inner < inner_iterations; ++inner) {
            for (std::size_t i = 0; i < n; ++i) detrended[i] = values[i] - trend[i];

            // Cycle-subseries smoothing: each seasonal position is smoothed
            // independently along its own subsequence.
            for (int phase = 0; phase < period; ++phase) {
                std::vector<std::size_t> indices;
                std::vector<double> subseries;
                for (std::size_t i = static_cast<std::size_t>(phase); i < n;
                     i += static_cast<std::size_t>(period)) {
                    indices.push_back(i);
                    subseries.push_back(detrended[i]);
                }
                if (subseries.empty()) continue;
                const std::vector<double> smoothed =
                    loess_smooth(subseries, seasonal_window, params.seasonal_degree, robustness);
                for (std::size_t k = 0; k < indices.size(); ++k) cycle[indices[k]] = smoothed[k];
            }

            // Low-pass filter of the cycle-subseries component, then remove it.
            std::vector<double> low = moving_average(cycle, period);
            low = moving_average(low, period);
            low = moving_average(low, 3);
            low = loess_smooth(low, low_pass_window, 1, robustness);
            for (std::size_t i = 0; i < n; ++i) seasonal[i] = cycle[i] - low[i];

            for (std::size_t i = 0; i < n; ++i) deseasonalized[i] = values[i] - seasonal[i];
            trend = loess_smooth(deseasonalized, trend_window, params.trend_degree, robustness);
        }

        for (std::size_t i = 0; i < n; ++i) result.remainder[i] = values[i] - seasonal[i] - trend[i];

        if (outer == outer_iterations) break;

        std::vector<double> magnitudes(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) magnitudes[i] = std::fabs(result.remainder[i]);
        const double scale = 6.0 * median_value(magnitudes);
        if (scale <= 1e-12) break;

        bool changed = false;
        for (std::size_t i = 0; i < n; ++i) {
            const double u = magnitudes[i] / scale;
            double weight = 0.0;
            if (u < 1.0) {
                const double t = 1.0 - u * u;
                weight = t * t;
            }
            if (std::fabs(weight - robustness[i]) > 1e-12) changed = true;
            robustness[i] = weight;
        }
        if (!changed) break;
    }

    result.seasonal = seasonal;
    result.trend = trend;
    result.robustness = robustness;
    for (std::size_t i = 0; i < n; ++i) result.remainder[i] = values[i] - seasonal[i] - trend[i];
    return result;
}

STLResult stl_decompose(const std::vector<double>& values, int period) {
    STLParams params;
    params.period = period;
    params.seasonal_window = std::max(7, period + 1);
    params.trend_window = std::max(21, 2 * period + 1);
    params.low_pass_window = std::max(7, period + 1);
    return stl_decompose(values, period, params);
}

}  // namespace climate
