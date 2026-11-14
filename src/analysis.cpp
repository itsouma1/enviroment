#include "analysis.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include "common.hpp"
#include "series.hpp"
#include "stats.hpp"

namespace climate {
namespace {

void sort_by_time(std::vector<double>& time, std::vector<double>& values) {
    std::vector<std::size_t> order(time.size());
    std::iota(order.begin(), order.end(), 0u);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
    std::vector<double> sorted_time(time.size());
    std::vector<double> sorted_values(values.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        sorted_time[i] = time[order[i]];
        sorted_values[i] = values[order[i]];
    }
    time.swap(sorted_time);
    values.swap(sorted_values);
}

double variance_of(const std::vector<double>& values) { return variance_value(values); }

}  // namespace

SeriesAnalysis analyze_series(const std::vector<double>& time_in, const std::vector<double>& values_in,
                              const SeriesOptions& options) {
    if (time_in.size() != values_in.size())
        throw std::runtime_error("analyze_series: time and values must have equal length");
    if (time_in.size() < 3)
        throw std::runtime_error("analyze_series: at least three observations are required");

    std::vector<double> time(time_in);
    std::vector<double> values(values_in);
    sort_by_time(time, values);

    SeriesAnalysis result;
    result.n = values.size();
    result.period = options.period;
    result.alpha = options.alpha;
    result.hamed_rao = options.hamed_rao;
    result.prewhiten = options.prewhiten;
    result.start_time = time.front();
    result.end_time = time.back();
    result.spacing = mean_spacing(time);

    std::size_t filled = 0;
    values = interpolate_missing(values, &filled);
    result.missing_filled = filled;
    if (filled > 0)
        result.warnings.push_back("filled " + std::to_string(filled) +
                                  " missing observations by linear interpolation");

    result.summary = summarize(values);

    const int period = options.period;
    const bool can_decompose =
        period >= 2 && values.size() >= static_cast<std::size_t>(2 * period);

    std::vector<double> seasonal(values.size(), 0.0);

    if (can_decompose) {
        STLParams params;
        params.period = period;
        params.seasonal_window = options.seasonal_window;
        params.trend_window = options.trend_window;
        params.low_pass_window = options.low_pass_window;
        params.seasonal_degree = options.seasonal_degree;
        params.trend_degree = options.trend_degree;
        params.inner_iterations = options.inner_iterations;
        params.outer_iterations = options.outer_iterations;
        try {
            result.decomposition = stl_decompose(values, period, params);
            result.decomposed = true;
            seasonal = result.decomposition.seasonal;

            const double seasonal_plus_remainder = variance_of([&]() {
                std::vector<double> combined(values.size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i)
                    combined[i] = result.decomposition.seasonal[i] + result.decomposition.remainder[i];
                return combined;
            }());
            const double trend_plus_remainder = variance_of([&]() {
                std::vector<double> combined(values.size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i)
                    combined[i] = result.decomposition.trend[i] + result.decomposition.remainder[i];
                return combined;
            }());
            const double remainder_variance = variance_of(result.decomposition.remainder);
            result.seasonal_strength =
                seasonal_plus_remainder > 0.0 ? std::max(0.0, 1.0 - remainder_variance / seasonal_plus_remainder)
                                              : 0.0;
            result.trend_strength =
                trend_plus_remainder > 0.0 ? std::max(0.0, 1.0 - remainder_variance / trend_plus_remainder)
                                           : 0.0;
            double lowest = seasonal.front();
            double highest = seasonal.front();
            for (double value : seasonal) {
                lowest = std::min(lowest, value);
                highest = std::max(highest, value);
            }
            result.seasonal_amplitude = highest - lowest;
        } catch (const std::exception& exception) {
            result.warnings.push_back(std::string("STL decomposition skipped: ") + exception.what());
        }
    } else {
        result.warnings.push_back(
            "series is shorter than two seasonal cycles; using a calendar-month climatology instead "
            "of STL");
    }

    if (!result.decomposed && period >= 2) {
        const std::vector<double> anomaly_series = remove_seasonal_cycle(values, static_cast<std::size_t>(period));
        for (std::size_t i = 0; i < values.size(); ++i) seasonal[i] = values[i] - anomaly_series[i];
        double lowest = seasonal.front();
        double highest = seasonal.front();
        for (double value : seasonal) {
            lowest = std::min(lowest, value);
            highest = std::max(highest, value);
        }
        result.seasonal_amplitude = highest - lowest;
    }

    result.deseasonalised.assign(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
        result.deseasonalised[i] = values[i] - seasonal[i];

    if (period >= 2) {
        const Climatology climatology = phase_climatology(values, static_cast<std::size_t>(period));
        result.standardised_anomalies = standardized_anomalies(values, climatology);
        if (!result.standardised_anomalies.empty()) {
            result.latest_anomaly = result.standardised_anomalies.back();
            result.minimum_anomaly = *std::min_element(result.standardised_anomalies.begin(),
                                                       result.standardised_anomalies.end());
            result.maximum_anomaly = *std::max_element(result.standardised_anomalies.begin(),
                                                       result.standardised_anomalies.end());
        }
    }

    MKSettings settings;
    settings.hamed_rao = options.hamed_rao;
    settings.alpha = options.alpha;

    std::vector<double> test_series = options.test_deseasonalised ? result.deseasonalised : values;
    if (options.prewhiten && period >= 2 && test_series.size() > 3)
        test_series = trend_free_prewhitening(time, test_series);

    result.mk = mann_kendall(test_series, settings);
    result.sen = sens_slope(time, test_series, nullptr, options.alpha);
    result.sen_per_decade = result.sen.computed
                                ? result.sen.slope * 10.0 * options.time_units_per_year
                                : 0.0;

    result.mk_raw = mann_kendall(values, settings);
    result.sen_raw = sens_slope(time, values, nullptr, options.alpha);
    result.sen_raw_per_decade = result.sen_raw.computed
                                    ? result.sen_raw.slope * 10.0 * options.time_units_per_year
                                    : 0.0;

    if (period >= 2) result.mk_seasonal = seasonal_mann_kendall(values, period, settings);

    if (options.hamed_rao && result.mk.correction_factor > 1.5)
        result.warnings.push_back(
            "serial correlation inflated the Mann-Kendall variance by a factor of " +
            std::to_string(result.mk.correction_factor) + "; the uncorrected test would overstate "
                                                        "significance");

    return result;
}

CO2Analysis analyze_co2(const std::vector<double>& time_in, const std::vector<double>& values_in,
                        const CO2Options& options) {
    if (time_in.size() != values_in.size())
        throw std::runtime_error("analyze_co2: time and values must have equal length");
    if (time_in.size() < 12)
        throw std::runtime_error("analyze_co2: at least twelve observations are required");

    std::vector<double> time(time_in);
    std::vector<double> values(values_in);
    sort_by_time(time, values);

    CO2Analysis result;
    result.n = values.size();
    result.start_time = time.front();
    result.end_time = time.back();

    std::size_t filled = 0;
    values = interpolate_missing(values, &filled);
    result.missing_filled = filled;
    if (filled > 0)
        result.warnings.push_back("filled " + std::to_string(filled) +
                                  " missing observations by linear interpolation");

    result.model = fit_harmonic(time, values, options.polynomial_order, options.harmonic_count,
                                options.period, options.reject_iterations,
                                options.reject_threshold);
    if (result.model.rejected_points > 0)
        result.warnings.push_back("rejected " + std::to_string(result.model.rejected_points) +
                                  " outliers during the iterative curve fit");

    result.fitted.resize(values.size());
    result.residual.resize(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        result.fitted[i] = harmonic_evaluate(result.model, time[i]);
        result.residual[i] = values[i] - result.fitted[i];
    }

    result.components = harmonic_components(result.model);
    result.seasonal_amplitude = harmonic_seasonal_amplitude(result.model);

    result.linear_coefficient = result.model.coefficients.size() > 1
                                    ? result.model.coefficients[1] / result.model.time_scale
                                    : 0.0;
    result.quadratic_coefficient =
        result.model.coefficients.size() > 2
            ? 2.0 * result.model.coefficients[2] / (result.model.time_scale * result.model.time_scale)
            : 0.0;

    const double start = time.front();
    const double end = time.back();
    result.growth_rate_start = harmonic_trend_derivative(result.model, start);
    result.growth_rate_end = harmonic_trend_derivative(result.model, end);
    result.growth_rate_mean = harmonic_mean_trend_derivative(result.model, start, end);
    const double recent_start = std::max(start, end - 5.0 * options.time_units_per_year);
    result.growth_rate_last_five_years =
        harmonic_mean_trend_derivative(result.model, recent_start, end);
    result.growth_rate_end_pgc_per_year =
        CARBON_PPM_TO_PGC * result.growth_rate_end;

    std::vector<double> deseasonalised(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
        deseasonalised[i] = values[i] - harmonic_seasonal_value(result.model, time[i]);

    MKSettings settings;
    settings.alpha = options.alpha;
    settings.hamed_rao = true;
    result.mk = mann_kendall(deseasonalised, settings);
    result.sen = sens_slope(time, deseasonalised, nullptr, options.alpha);

    return result;
}

}  // namespace climate
