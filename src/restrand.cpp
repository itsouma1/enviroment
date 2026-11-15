#include "restrand.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

#include "linalg.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<double> project_simplex(std::vector<double> values) {
    const std::size_t n = values.size();
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end(), std::greater<double>());
    double cumulative = 0.0;
    std::size_t rho = 0;
    for (std::size_t i = 0; i < n; ++i) {
        cumulative += sorted[i];
        if (sorted[i] - (cumulative - 1.0) / static_cast<double>(i + 1) > 0.0) rho = i + 1;
    }
    double threshold = 0.0;
    for (std::size_t i = 0; i < rho; ++i) threshold += sorted[i];
    threshold = (threshold - 1.0) / static_cast<double>(rho);
    for (double& value : values) value = std::max(0.0, value - threshold);
    return values;
}

double correlation(const std::vector<double>& x, const std::vector<double>& y) {
    const double n = static_cast<double>(x.size());
    double mean_x = 0.0, mean_y = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }
    mean_x /= n;
    mean_y /= n;
    double covariance = 0.0, variance_x = 0.0, variance_y = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        covariance += (x[i] - mean_x) * (y[i] - mean_y);
        variance_x += (x[i] - mean_x) * (x[i] - mean_x);
        variance_y += (y[i] - mean_y) * (y[i] - mean_y);
    }
    return covariance / std::sqrt(std::max(1e-300, variance_x * variance_y));
}

}  // namespace

std::vector<double> restrand_residuals(const std::vector<double>& response,
                                       const Matrix& predictors) {
    require(!response.empty() && predictors.size() == response.size(),
            "restrand_residuals: size");
    const std::size_t n = response.size(), columns = predictors[0].size() + 1;
    Matrix design(n, std::vector<double>(columns, 1.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < predictors[i].size(); ++j)
            design[i][j + 1] = predictors[i][j];
    std::vector<double> coefficients;
    require(least_squares(design, response, coefficients),
            "restrand_residuals: rank deficient");
    std::vector<double> residuals(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        residuals[i] = response[i];
        for (std::size_t j = 0; j < coefficients.size(); ++j)
            residuals[i] -= coefficients[j] * design[i][j];
    }
    return residuals;
}

double residual_trend(const std::vector<double>& residuals) {
    require(residuals.size() >= 2, "residual_trend: size");
    const double n = static_cast<double>(residuals.size());
    const double mean_t = 0.5 * (n - 1.0);
    double mean_y = 0.0;
    for (double value : residuals) mean_y += value;
    mean_y /= n;
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t t = 0; t < residuals.size(); ++t) {
        const double dt = static_cast<double>(t) - mean_t;
        numerator += dt * (residuals[t] - mean_y);
        denominator += dt * dt;
    }
    return numerator / denominator;
}

SegmentedTrend segmented_trend(const std::vector<double>& times,
                               const std::vector<double>& values, double breakpoint) {
    require(times.size() == values.size() && times.size() >= 4,
            "segmented_trend: size");
    const std::size_t n = times.size();
    Matrix design(n, std::vector<double>(4, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        const double post = times[i] > breakpoint ? 1.0 : 0.0;
        design[i][0] = 1.0;
        design[i][1] = times[i];
        design[i][2] = post;
        design[i][3] = times[i] * post;
    }
    std::vector<double> coefficients;
    require(least_squares(design, values, coefficients), "segmented_trend: rank deficient");
    return {coefficients[0], coefficients[1], coefficients[2], coefficients[3]};
}

std::vector<double> sdid_unit_weights(const Matrix& control_outcomes,
                                      const std::vector<double>& treated_mean,
                                      double ridge) {
    require(!control_outcomes.empty() &&
                control_outcomes.size() == treated_mean.size() && ridge >= 0.0,
            "sdid_unit_weights: size");
    const std::size_t times = control_outcomes.size();
    const std::size_t units = control_outcomes[0].size();
    require(units > 0, "sdid_unit_weights: units");
    std::vector<double> weights(units, 1.0 / static_cast<double>(units));
    const double step = 1.0 / (static_cast<double>(times) * 10.0 + ridge + 1.0);
    for (int iteration = 0; iteration < 2000; ++iteration) {
        std::vector<double> gradient(units, 0.0);
        for (std::size_t t = 0; t < times; ++t) {
            double synthetic = 0.0;
            for (std::size_t i = 0; i < units; ++i) synthetic += weights[i] * control_outcomes[t][i];
            const double residual = synthetic - treated_mean[t];
            for (std::size_t i = 0; i < units; ++i)
                gradient[i] += 2.0 * residual * control_outcomes[t][i];
        }
        for (std::size_t i = 0; i < units; ++i) gradient[i] += 2.0 * ridge * weights[i];
        for (std::size_t i = 0; i < units; ++i) weights[i] -= step * gradient[i];
        weights = project_simplex(std::move(weights));
    }
    return weights;
}

double sdid_effect(const Matrix& outcomes, const std::vector<int>& treated_units,
                   std::size_t treatment_time, const std::vector<double>& unit_weights,
                   const std::vector<double>& time_weights) {
    require(!outcomes.empty() && outcomes[0].size() == treated_units.size(),
            "sdid_effect: size");
    const std::size_t times = outcomes.size(), units = treated_units.size();
    require(treatment_time > 0 && treatment_time < times, "sdid_effect: treatment time");
    std::vector<std::size_t> controls;
    for (std::size_t i = 0; i < units; ++i)
        if (!treated_units[i]) controls.push_back(i);
    require(!controls.empty() && unit_weights.size() == controls.size(),
            "sdid_effect: controls");
    require(time_weights.size() == treatment_time, "sdid_effect: time weights");
    double weight_sum = 0.0;
    for (double weight : time_weights) weight_sum += weight;
    require(weight_sum > 0.0, "sdid_effect: zero time weight");
    auto treated_mean = [&](std::size_t start, std::size_t end) {
        double sum = 0.0;
        std::size_t count = 0;
        for (std::size_t t = start; t < end; ++t)
            for (std::size_t i = 0; i < units; ++i)
                if (treated_units[i]) {
                    sum += outcomes[t][i];
                    ++count;
                }
        return sum / static_cast<double>(count);
    };
    auto synthetic_mean = [&](std::size_t start, std::size_t end) {
        double sum = 0.0;
        std::size_t count = 0;
        for (std::size_t t = start; t < end; ++t) {
            double value = 0.0;
            for (std::size_t k = 0; k < controls.size(); ++k)
                value += unit_weights[k] * outcomes[t][controls[k]];
            sum += value;
            ++count;
        }
        return sum / static_cast<double>(count);
    };
    double treated_pre = 0.0, synthetic_pre = 0.0;
    for (std::size_t t = 0; t < treatment_time; ++t) {
        double treated_t = 0.0;
        std::size_t count = 0;
        for (std::size_t i = 0; i < units; ++i)
            if (treated_units[i]) {
                treated_t += outcomes[t][i];
                ++count;
            }
        treated_t /= static_cast<double>(count);
        double synthetic_t = 0.0;
        for (std::size_t k = 0; k < controls.size(); ++k)
            synthetic_t += unit_weights[k] * outcomes[t][controls[k]];
        treated_pre += time_weights[t] * treated_t;
        synthetic_pre += time_weights[t] * synthetic_t;
    }
    treated_pre /= weight_sum;
    synthetic_pre /= weight_sum;
    const double treated_post = treated_mean(treatment_time, times);
    const double synthetic_post = synthetic_mean(treatment_time, times);
    return (treated_post - treated_pre) - (synthetic_post - synthetic_pre);
}

double placebo_p_value(std::vector<double> placebo_effects, double observed_effect) {
    require(!placebo_effects.empty(), "placebo_p_value: empty");
    std::size_t exceedances = 0;
    for (double effect : placebo_effects)
        if (effect >= observed_effect) ++exceedances;
    return static_cast<double>(exceedances + 1) /
           static_cast<double>(placebo_effects.size() + 1);
}

std::vector<double> stability_selection(const Matrix& predictors,
                                        const std::vector<double>& response,
                                        std::size_t subsamples, double threshold,
                                        unsigned seed) {
    require(!predictors.empty() && predictors.size() == response.size() && subsamples > 0,
            "stability_selection: size");
    const std::size_t n = response.size(), variables = predictors[0].size();
    std::mt19937 generator(seed);
    std::vector<double> counts(variables, 0.0);
    const std::size_t subset = std::max<std::size_t>(5, n / 2);
    std::vector<std::size_t> index(n);
    std::iota(index.begin(), index.end(), 0);
    for (std::size_t replicate = 0; replicate < subsamples; ++replicate) {
        std::shuffle(index.begin(), index.end(), generator);
        std::vector<double> scores(variables, 0.0);
        for (std::size_t j = 0; j < variables; ++j) {
            std::vector<double> x(subset, 0.0), y(subset, 0.0);
            for (std::size_t i = 0; i < subset; ++i) {
                x[i] = predictors[index[i]][j];
                y[i] = response[index[i]];
            }
            scores[j] = std::fabs(correlation(x, y));
        }
        std::vector<double> sorted = scores;
        std::sort(sorted.begin(), sorted.end(), std::greater<double>());
        const double cutoff = sorted[std::min(variables - 1, variables / 2)];
        for (std::size_t j = 0; j < variables; ++j)
            if (scores[j] >= cutoff && scores[j] >= threshold) counts[j] += 1.0;
    }
    for (double& count : counts) count /= static_cast<double>(subsamples);
    return counts;
}

std::size_t knockoff_threshold(std::vector<double> statistics, double q) {
    require(!statistics.empty() && q > 0.0 && q < 1.0, "knockoff_threshold: inputs");
    std::sort(statistics.begin(), statistics.end());
    std::size_t best = statistics.size();
    for (std::size_t k = 0; k < statistics.size(); ++k) {
        const double t = statistics[k];
        if (t <= 0.0) continue;
        std::size_t positive = 0, negative = 0;
        for (double value : statistics) {
            if (value >= t) ++positive;
            if (value <= -t) ++negative;
        }
        if (positive > 0 &&
            static_cast<double>(1 + negative) / static_cast<double>(positive) <= q) {
            best = k;
            break;
        }
    }
    return best;
}

std::vector<std::size_t> knapsack_select(const std::vector<double>& weights,
                                         const std::vector<int>& costs, int budget) {
    require(weights.size() == costs.size() && budget >= 0, "knapsack_select: inputs");
    const std::size_t n = weights.size();
    Matrix best(n + 1, std::vector<double>(static_cast<std::size_t>(budget) + 1, 0.0));
    for (std::size_t i = 1; i <= n; ++i)
        for (int capacity = 0; capacity <= budget; ++capacity) {
            best[i][static_cast<std::size_t>(capacity)] = best[i - 1][static_cast<std::size_t>(capacity)];
            if (costs[i - 1] <= capacity)
                best[i][static_cast<std::size_t>(capacity)] = std::max(
                    best[i][static_cast<std::size_t>(capacity)],
                    best[i - 1][static_cast<std::size_t>(capacity - costs[i - 1])] +
                        weights[i - 1]);
        }
    std::vector<std::size_t> chosen;
    int capacity = budget;
    for (std::size_t i = n; i > 0; --i)
        if (best[i][static_cast<std::size_t>(capacity)] !=
            best[i - 1][static_cast<std::size_t>(capacity)]) {
            chosen.push_back(i - 1);
            capacity -= costs[i - 1];
        }
    std::sort(chosen.begin(), chosen.end());
    return chosen;
}

double theil_sen_slope(const std::vector<double>& times, const std::vector<double>& values) {
    require(times.size() == values.size() && times.size() >= 2, "theil_sen_slope: size");
    std::vector<double> slopes;
    for (std::size_t i = 0; i < times.size(); ++i)
        for (std::size_t j = i + 1; j < times.size(); ++j)
            if (times[j] != times[i])
                slopes.push_back((values[j] - values[i]) / (times[j] - times[i]));
    require(!slopes.empty(), "theil_sen_slope: no pairs");
    const std::size_t middle = slopes.size() / 2;
    std::nth_element(slopes.begin(), slopes.begin() + static_cast<long>(middle), slopes.end());
    return slopes[middle];
}

double pettitt_statistic(const std::vector<double>& values, std::size_t& change_point) {
    require(values.size() >= 3, "pettitt_statistic: size");
    double maximum = 0.0;
    change_point = 0;
    for (std::size_t t = 0; t + 1 < values.size(); ++t) {
        double sum = 0.0;
        for (std::size_t i = 0; i <= t; ++i)
            for (std::size_t j = t + 1; j < values.size(); ++j)
                sum += values[j] > values[i] ? 1.0 : (values[j] < values[i] ? -1.0 : 0.0);
        if (std::fabs(sum) > std::fabs(maximum)) {
            maximum = sum;
            change_point = t;
        }
    }
    return maximum;
}

std::vector<std::size_t> benjamini_hochberg(const std::vector<double>& p_values, double q) {
    require(q > 0.0 && q < 1.0, "benjamini_hochberg: q");
    std::vector<std::size_t> order(p_values.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return p_values[a] < p_values[b]; });
    std::size_t largest = 0;
    for (std::size_t k = 0; k < order.size(); ++k)
        if (p_values[order[k]] <=
            q * static_cast<double>(k + 1) / static_cast<double>(order.size()))
            largest = k + 1;
    std::vector<std::size_t> accepted(order.begin(), order.begin() + static_cast<long>(largest));
    std::sort(accepted.begin(), accepted.end());
    return accepted;
}

}  // namespace climate
