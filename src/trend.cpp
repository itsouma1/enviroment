#include "trend.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include "stats.hpp"

namespace climate {
namespace {

double tie_corrected_variance(const std::vector<double>& values, double s_unused) {
    (void)s_unused;
    const std::size_t n = values.size();
    if (n < 2) return 0.0;

    double variance = static_cast<double>(n) * static_cast<double>(n - 1) *
                      static_cast<double>(2 * n + 5) / 18.0;

    std::vector<double> sorted(values);
    std::sort(sorted.begin(), sorted.end());
    std::size_t i = 0;
    while (i < sorted.size()) {
        std::size_t j = i + 1;
        while (j < sorted.size() && sorted[j] == sorted[i]) ++j;
        const double tie_length = static_cast<double>(j - i);
        if (tie_length > 1.0) {
            variance -= tie_length * (tie_length - 1.0) * (2.0 * tie_length + 5.0) / 18.0;
        }
        i = j;
    }
    return variance;
}

std::vector<double> ranks(const std::vector<double>& values) {
    const std::size_t n = values.size();
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });

    std::vector<double> result(n, 0.0);
    std::size_t i = 0;
    while (i < n) {
        std::size_t j = i + 1;
        while (j < n && values[order[j]] == values[order[i]]) ++j;
        const double average_rank = (static_cast<double>(i) + static_cast<double>(j - 1)) / 2.0 + 1.0;
        for (std::size_t k = i; k < j; ++k) result[order[k]] = average_rank;
        i = j;
    }
    return result;
}

double hamed_rao_factor(const std::vector<double>& values) {
    const std::size_t n = values.size();
    if (n < 4) return 1.0;

    // Estimate the autocorrelation from the detrended series. A deterministic
    // trend would otherwise produce near-unit autocorrelation at every lag and
    // an absurd variance correction.
    std::vector<double> index(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) index[i] = static_cast<double>(i);
    const double slope = ols_slope(index, values);
    std::vector<double> detrended(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        detrended[i] = values[i] - slope * static_cast<double>(i);

    const std::vector<double> ranked = ranks(detrended);
    const double limit = 1.96 / std::sqrt(static_cast<double>(n));
    double sum = 0.0;
    for (std::size_t k = 0; k + 1 < n; ++k) {
        const int lag = static_cast<int>(k) + 1;
        const double rho = autocorrelation(ranked, lag);
        if (std::fabs(rho) <= limit) continue;
        const double nk = static_cast<double>(n - k);
        sum += (nk) * (nk - 1.0) * (nk - 2.0) * rho;
    }
    const double denominator = static_cast<double>(n) * static_cast<double>(n - 1) *
                               static_cast<double>(n - 2);
    if (denominator <= 0.0) return 1.0;
    double factor = 1.0 + 2.0 * sum / denominator;
    if (factor < 1.0) factor = 1.0;
    return factor;
}

}  // namespace

MKResult mann_kendall(const std::vector<double>& values, const MKSettings& settings) {
    MKResult result;
    result.n = values.size();
    result.critical_value = normal_quantile(1.0 - settings.alpha / 2.0);

    if (values.size() < 3) {
        result.variance = 0.0;
        result.variance_corrected = 0.0;
        result.p_value = 1.0;
        result.direction = "insufficient data";
        return result;
    }

    double s = 0.0;
    const std::size_t n = values.size();
    for (std::size_t k = 0; k + 1 < n; ++k) {
        for (std::size_t j = k + 1; j < n; ++j) {
            if (values[j] > values[k])
                s += 1.0;
            else if (values[j] < values[k])
                s -= 1.0;
        }
    }
    result.s = s;
    result.variance = tie_corrected_variance(values, s);

    double variance = result.variance;
    result.correction_factor = 1.0;
    if (settings.hamed_rao) {
        result.correction_factor = hamed_rao_factor(values);
        variance *= result.correction_factor;
    }
    result.variance_corrected = variance;

    if (variance <= 0.0) {
        result.z = 0.0;
        result.p_value = 1.0;
        result.tau = 0.0;
        result.direction = "no trend";
        return result;
    }

    if (s > 0.0)
        result.z = (s - 1.0) / std::sqrt(variance);
    else if (s < 0.0)
        result.z = (s + 1.0) / std::sqrt(variance);
    else
        result.z = 0.0;

    result.p_value = 2.0 * (1.0 - normal_cdf(std::fabs(result.z)));
    if (result.p_value < 0.0) result.p_value = 0.0;
    if (result.p_value > 1.0) result.p_value = 1.0;
    result.tau = s / (0.5 * static_cast<double>(n) * static_cast<double>(n - 1));
    result.significant = std::fabs(result.z) > result.critical_value;
    if (!result.significant)
        result.direction = "no significant trend";
    else
        result.direction = result.z > 0.0 ? "increasing" : "decreasing";
    return result;
}

SenResult sens_slope(const std::vector<double>& time, const std::vector<double>& values, const SenResult*,
                     double alpha) {
    if (time.size() != values.size())
        throw std::runtime_error("sens_slope: time and values must have equal length");
    SenResult result;
    result.alpha = alpha;
    const std::size_t n = time.size();
    if (n < 3) return result;

    std::vector<double> slopes;
    slopes.reserve(n * (n - 1) / 2);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double dt = time[j] - time[i];
            if (dt == 0.0) continue;
            slopes.push_back((values[j] - values[i]) / dt);
        }
    }
    if (slopes.empty()) return result;

    result.slope = median_value(slopes);
    result.pairs = slopes.size();
    result.computed = true;

    // Intercept at the mean time, which is well conditioned.
    const double time_mean = mean_value(time);
    double value_mean = 0.0;
    for (double v : values) value_mean += v;
    value_mean /= static_cast<double>(n);
    result.intercept = value_mean - result.slope * time_mean;

    const MKResult mk = mann_kendall(values, MKSettings{false, false, alpha});
    const double critical = normal_quantile(1.0 - alpha / 2.0);
    const double rank_offset = critical * std::sqrt(std::max(0.0, mk.variance));
    const double pair_count = static_cast<double>(slopes.size());

    std::sort(slopes.begin(), slopes.end());
    long lower_index = static_cast<long>(std::floor((pair_count - rank_offset) / 2.0));
    long upper_index = static_cast<long>(std::ceil((pair_count + rank_offset) / 2.0)) - 1;

    if (lower_index < 0) lower_index = 0;
    if (upper_index >= static_cast<long>(slopes.size()))
        upper_index = static_cast<long>(slopes.size()) - 1;
    if (upper_index < lower_index) upper_index = lower_index;

    result.lower = slopes[static_cast<std::size_t>(lower_index)];
    result.upper = slopes[static_cast<std::size_t>(upper_index)];
    return result;
}

std::vector<double> trend_free_prewhitening(const std::vector<double>& time,
                                           const std::vector<double>& values) {
    const std::size_t n = values.size();
    if (n < 3) return values;

    const SenResult sen = sens_slope(time, values, nullptr, 0.05);
    if (!sen.computed) return values;

    std::vector<double> detrended(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) detrended[i] = values[i] - sen.slope * time[i];

    const double rho = autocorrelation(detrended, 1);
    if (!std::isfinite(rho) || std::fabs(rho) >= 1.0) return values;

    std::vector<double> whitened(n, 0.0);
    whitened[0] = detrended[0];
    for (std::size_t i = 1; i < n; ++i) whitened[i] = detrended[i] - rho * detrended[i - 1];

    std::vector<double> restored(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) restored[i] = whitened[i] + sen.slope * time[i];
    return restored;
}

MKResult seasonal_mann_kendall(const std::vector<double>& values, int period,
                               const MKSettings& settings) {
    if (period < 2) throw std::runtime_error("seasonal_mann_kendall: period must be >= 2");
    const std::size_t n = values.size();
    if (n < static_cast<std::size_t>(2 * period)) return mann_kendall(values, settings);

    double s = 0.0;
    double variance = 0.0;
    std::size_t total_used = 0;

    for (int phase = 0; phase < period; ++phase) {
        std::vector<double> block;
        for (std::size_t i = static_cast<std::size_t>(phase); i < n;
             i += static_cast<std::size_t>(period))
            block.push_back(values[i]);
        const std::size_t m = block.size();
        total_used += m;
        for (std::size_t k = 0; k + 1 < m; ++k) {
            for (std::size_t j = k + 1; j < m; ++j) {
                if (block[j] > block[k])
                    s += 1.0;
                else if (block[j] < block[k])
                    s -= 1.0;
            }
        }
        variance += tie_corrected_variance(block, 0.0);
    }

    MKResult result;
    result.n = total_used;
    result.s = s;
    result.variance = variance;
    result.variance_corrected = variance;
    result.critical_value = normal_quantile(1.0 - settings.alpha / 2.0);
    if (variance <= 0.0) {
        result.p_value = 1.0;
        result.direction = "no trend";
        return result;
    }
    if (s > 0.0)
        result.z = (s - 1.0) / std::sqrt(variance);
    else if (s < 0.0)
        result.z = (s + 1.0) / std::sqrt(variance);
    else
        result.z = 0.0;
    result.p_value = 2.0 * (1.0 - normal_cdf(std::fabs(result.z)));
    result.tau = 0.0;
    result.significant = std::fabs(result.z) > result.critical_value;
    result.direction = result.significant ? (result.z > 0.0 ? "increasing" : "decreasing")
                                          : "no significant trend";
    return result;
}

}  // namespace climate
