#include "spatial.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace climate {
namespace {

struct Moments {
    double mean = 0.0;
    double variance = 0.0;  // population variance (divide by n)
    std::size_t count = 0;
};

Moments moments(const std::vector<double>& values) {
    Moments result;
    result.count = values.size();
    if (values.empty()) return result;
    double sum = 0.0;
    for (const double value : values) sum += value;
    result.mean = sum / static_cast<double>(values.size());
    double squares = 0.0;
    for (const double value : values) {
        const double deviation = value - result.mean;
        squares += deviation * deviation;
    }
    result.variance = squares / static_cast<double>(values.size());
    return result;
}

void check_shape(const std::vector<double>& values, const Weights& weights) {
    if (values.empty()) throw std::runtime_error("spatial statistic: no values");
    if (weights.size() != values.size())
        throw std::runtime_error("spatial statistic: weight matrix must be n x n");
    for (const auto& row : weights) {
        if (row.size() != values.size())
            throw std::runtime_error("spatial statistic: weight matrix must be n x n");
    }
}

}  // namespace

double getis_ord_gi(const std::vector<double>& values, const Weights& weights,
                    std::size_t index) {
    check_shape(values, weights);
    if (index >= values.size()) throw std::runtime_error("getis_ord_gi: index out of range");
    const Moments stats = moments(values);
    if (stats.count < 2) throw std::runtime_error("getis_ord_gi: need at least two values");

    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    double weight_squares = 0.0;
    for (std::size_t j = 0; j < values.size(); ++j) {
        const double w = weights[index][j];
        weighted_sum += w * values[j];
        weight_sum += w;
        weight_squares += w * w;
    }

    const double numerator = weighted_sum - stats.mean * weight_sum;
    const double denominator =
        std::sqrt(stats.variance) *
        std::sqrt((stats.count * weight_squares - weight_sum * weight_sum) /
                  static_cast<double>(stats.count - 1));
    if (denominator == 0.0 || !std::isfinite(numerator)) return std::nan("");
    return numerator / denominator;
}

double local_moran_i(const std::vector<double>& values, const Weights& weights,
                     std::size_t index) {
    check_shape(values, weights);
    if (index >= values.size()) throw std::runtime_error("local_moran_i: index out of range");
    const Moments stats = moments(values);
    if (stats.variance == 0.0) return std::nan("");
    double lagged = 0.0;
    for (std::size_t j = 0; j < values.size(); ++j) {
        lagged += weights[index][j] * (values[j] - stats.mean);
    }
    return (values[index] - stats.mean) / stats.variance * lagged;
}

double geodetector_q(const std::vector<double>& values, const std::vector<int>& strata) {
    if (values.size() != strata.size())
        throw std::runtime_error("geodetector_q: values and strata must have equal length");
    if (values.empty()) throw std::runtime_error("geodetector_q: no values");
    const Moments total = moments(values);
    if (total.variance == 0.0) return std::nan("");

    std::vector<int> labels = strata;
    std::sort(labels.begin(), labels.end());
    labels.erase(std::unique(labels.begin(), labels.end()), labels.end());

    double within = 0.0;
    for (const int label : labels) {
        std::vector<double> subset;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (strata[i] == label) subset.push_back(values[i]);
        }
        const Moments stratum = moments(subset);
        within += static_cast<double>(stratum.count) * stratum.variance;
    }
    const double total_scatter = static_cast<double>(values.size()) * total.variance;
    if (total_scatter == 0.0) return std::nan("");
    return 1.0 - within / total_scatter;
}

double walker_threshold(double alpha_global, int hypothesis_count) {
    if (!(alpha_global > 0.0 && alpha_global < 1.0))
        throw std::runtime_error("walker_threshold: alpha must be in (0,1)");
    if (hypothesis_count <= 0)
        throw std::runtime_error("walker_threshold: hypothesis count must be positive");
    return 1.0 - std::pow(1.0 - alpha_global, 1.0 / static_cast<double>(hypothesis_count));
}

std::vector<double> benjamini_hochberg(const std::vector<double>& p_values) {
    const std::size_t m = p_values.size();
    std::vector<std::size_t> order(m);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return p_values[a] < p_values[b]; });

    std::vector<double> adjusted(m, 1.0);
    double running_min = 1.0;
    for (std::size_t rank = m; rank-- > 0;) {
        const double candidate =
            p_values[order[rank]] * static_cast<double>(m) / static_cast<double>(rank + 1);
        running_min = std::min(running_min, candidate);
        adjusted[order[rank]] = std::min(1.0, running_min);
    }
    return adjusted;
}

std::vector<std::size_t> benjamini_hochberg_rejections(const std::vector<double>& p_values,
                                                       double false_discovery_rate) {
    if (!(false_discovery_rate > 0.0 && false_discovery_rate < 1.0))
        throw std::runtime_error("benjamini_hochberg: FDR must be in (0,1)");
    const std::vector<double> adjusted = benjamini_hochberg(p_values);
    std::vector<std::size_t> rejected;
    for (std::size_t i = 0; i < adjusted.size(); ++i) {
        if (adjusted[i] <= false_discovery_rate) rejected.push_back(i);
    }
    return rejected;
}

double poisson_scan_log_likelihood(int inside_count, double inside_expected, int total_count,
                                   double total_expected) {
    if (total_count <= 0) throw std::runtime_error("poisson_scan: total count must be positive");
    if (!(inside_expected > 0.0)) throw std::runtime_error("poisson_scan: inside expected must be > 0");
    if (!(total_expected > inside_expected))
        throw std::runtime_error("poisson_scan: total expected must exceed the inside expected");
    if (inside_count < 0 || inside_count > total_count)
        throw std::runtime_error("poisson_scan: inside count must be within [0, total]");

    const int outside_count = total_count - inside_count;
    const double outside_expected = total_expected - inside_expected;
    if (inside_count == 0 || outside_count == 0) return 0.0;

    const double inside = static_cast<double>(inside_count) *
                          std::log(static_cast<double>(inside_count) / inside_expected);
    const double outside = static_cast<double>(outside_count) *
                           std::log(static_cast<double>(outside_count) / outside_expected);
    return inside + outside;
}

}  // namespace climate
