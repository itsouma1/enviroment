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

namespace {

double rank_correlation(const std::vector<double>& first, const std::vector<double>& second) {
    const std::size_t n = first.size();
    if (n < 3) return 0.0;
    auto ranks = [](const std::vector<double>& values) {
        std::vector<std::size_t> order(values.size());
        std::iota(order.begin(), order.end(), 0u);
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });
        std::vector<double> rank(values.size());
        std::size_t i = 0;
        while (i < order.size()) {
            std::size_t j = i + 1;
            while (j < order.size() && values[order[j]] == values[order[i]]) ++j;
            const double average = 0.5 * static_cast<double>(i + j - 1) + 1.0;
            for (std::size_t k = i; k < j; ++k) rank[order[k]] = average;
            i = j;
        }
        return rank;
    };
    const std::vector<double> a = ranks(first);
    const std::vector<double> b = ranks(second);
    const double mean_a = mean_value(a);
    const double mean_b = mean_value(b);
    double covariance = 0.0;
    double variance_a = 0.0;
    double variance_b = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        covariance += (a[i] - mean_a) * (b[i] - mean_b);
        variance_a += (a[i] - mean_a) * (a[i] - mean_a);
        variance_b += (b[i] - mean_b) * (b[i] - mean_b);
    }
    if (variance_a <= 0.0 || variance_b <= 0.0) return 0.0;
    return covariance / std::sqrt(variance_a * variance_b);
}

double line_rmse(const std::vector<double>& time, const std::vector<double>& values, double slope,
                 double intercept) {
    double sse = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double residual = values[i] - (intercept + slope * time[i]);
        sse += residual * residual;
    }
    return std::sqrt(sse / static_cast<double>(values.size()));
}

double two_sided_normal(double z) { return 2.0 * (1.0 - normal_cdf(std::fabs(z))); }

}  // namespace

TrendEnsemble trend_ensemble(const std::vector<double>& time, const std::vector<double>& values,
                             double alpha) {
    TrendEnsemble ensemble;
    if (time.size() != values.size() || time.size() < 4) return ensemble;
    const std::size_t n = values.size();
    const double scale = std::max(stdev_value(values), 1e-12);
    const double mean_time = mean_value(time);
    const double mean_values = mean_value(values);

    const SenResult sen = sens_slope(time, values, nullptr, alpha);
    const double sen_slope = sen.computed ? sen.slope : 0.0;
    const double sen_intercept = sen.computed ? sen.intercept : mean_values;

    auto add = [&](const std::string& name, const std::string& formula, double slope,
                   double intercept, double p_value, bool computed) {
        TrendEstimator member;
        member.name = name;
        member.formula = formula;
        member.slope = slope;
        member.intercept = intercept;
        member.p_value = computed ? p_value : 1.0;
        member.computed = computed;
        member.significant = computed && p_value <= alpha;
        member.rmse = computed ? line_rmse(time, values, slope, intercept) : 0.0;
        const double accuracy = 1.0 / (1.0 + member.rmse / scale);
        const double confidence = 1.0 - std::min(1.0, member.p_value);
        const double sign_agreement = std::signbit(slope) == std::signbit(sen_slope) ? 1.0 : 0.0;
        member.score = computed ? 0.40 * accuracy + 0.30 * confidence + 0.20 * sign_agreement +
                                      0.10 * (member.significant ? 1.0 : 0.0)
                                : 0.0;
        ensemble.members.push_back(member);
    };

    const MKResult kendall = mann_kendall(values);
    add("Sen slope", "median pairwise slope", sen_slope, sen_intercept, kendall.p_value, sen.computed);

    const double ols = ols_slope(time, values);
    const double ols_intercept = mean_values - ols * mean_time;
    double sxx = 0.0;
    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sxx += (time[i] - mean_time) * (time[i] - mean_time);
        const double residual = values[i] - (ols_intercept + ols * time[i]);
        sse += residual * residual;
    }
    const double ols_se = sxx > 0.0 && n > 2 ? std::sqrt(sse / static_cast<double>(n - 2) / sxx) : 0.0;
    const double ols_p = ols_se > 0.0 ? two_sided_normal(ols / ols_se) : 1.0;
    add("OLS slope", "least squares", ols, ols_intercept, ols_p, sxx > 0.0);

    add("Kendall tau", "rank concordance, Sen slope", sen_slope, sen_intercept, kendall.p_value,
        sen.computed);

    const double spearman = rank_correlation(time, values);
    const double spearman_z = spearman * std::sqrt(static_cast<double>(n - 1));
    add("Spearman rho", "rank correlation, Sen slope", sen_slope, sen_intercept,
        two_sided_normal(spearman_z), sen.computed);

    const std::size_t half = n / 2;
    std::size_t positive = 0;
    std::size_t negative = 0;
    for (std::size_t i = 0; i < half; ++i) {
        if (values[i + n - half] > values[i]) ++positive;
        if (values[i + n - half] < values[i]) ++negative;
    }
    const double cox_s = static_cast<double>(positive) - static_cast<double>(negative);
    const double cox_z =
        half > 0 ? (std::fabs(cox_s) - 1.0) / std::sqrt(static_cast<double>(half)) : 0.0;
    const double cox_slope = half > 0
                                 ? (median_value(std::vector<double>(values.end() - half, values.end())) -
                                    median_value(std::vector<double>(values.begin(), values.begin() + half))) /
                                       (mean_value(std::vector<double>(time.end() - half, time.end())) -
                                        mean_value(std::vector<double>(time.begin(), time.begin() + half)))
                                 : 0.0;
    add("Cox-Stuart", "split-half sign test", cox_slope,
        mean_values - cox_slope * mean_time, two_sided_normal(cox_z), half >= 2);

    const std::size_t third = n / 3;
    double first_mean = 0.0;
    double last_mean = 0.0;
    double first_var = 0.0;
    double last_var = 0.0;
    if (third >= 2) {
        for (std::size_t i = 0; i < third; ++i) {
            first_mean += values[i];
            last_mean += values[n - third + i];
        }
        first_mean /= static_cast<double>(third);
        last_mean /= static_cast<double>(third);
        for (std::size_t i = 0; i < third; ++i) {
            first_var += (values[i] - first_mean) * (values[i] - first_mean);
            last_var += (values[n - third + i] - last_mean) * (values[n - third + i] - last_mean);
        }
        first_var /= static_cast<double>(third - 1);
        last_var /= static_cast<double>(third - 1);
    }
    const double endpoint_se = third >= 2 ? std::sqrt((first_var + last_var) / static_cast<double>(third)) : 0.0;
    const double endpoint_slope =
        third >= 2 ? (last_mean - first_mean) /
                         (mean_value(std::vector<double>(time.end() - third, time.end())) -
                          mean_value(std::vector<double>(time.begin(), time.begin() + third)))
                   : 0.0;
    add("Endpoint difference", "first third vs last third", endpoint_slope,
        mean_values - endpoint_slope * mean_time,
        endpoint_se > 0.0 ? two_sided_normal((last_mean - first_mean) / endpoint_se) : 1.0,
        third >= 2);

    double best = -1.0;
    for (std::size_t i = 0; i < ensemble.members.size(); ++i) {
        if (!ensemble.members[i].computed) continue;
        if (ensemble.members[i].score > best + 0.01) {
            best = ensemble.members[i].score;
            ensemble.selected = i;
        }
    }
    if (best < 0.0) return ensemble;
    ensemble.members[ensemble.selected].selected = true;
    ensemble.robust_slope = ensemble.members[ensemble.selected].slope;
    ensemble.robust_p_value = ensemble.members[ensemble.selected].p_value;
    std::size_t agree = 0;
    std::size_t total = 0;
    for (const auto& member : ensemble.members) {
        if (!member.computed) continue;
        ++total;
        if (std::signbit(member.slope) == std::signbit(ensemble.robust_slope)) ++agree;
    }
    ensemble.agreement = total > 0 ? static_cast<double>(agree) / static_cast<double>(total) : 0.0;
    ensemble.consensus = ensemble.members[ensemble.selected].significant
                             ? (ensemble.robust_slope > 0.0 ? "increasing" : "decreasing")
                             : "no significant trend";
    return ensemble;
}

}  // namespace climate
