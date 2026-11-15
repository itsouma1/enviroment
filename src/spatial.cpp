#include "spatial.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "linalg.hpp"

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

LisaClass lisa_classification(const std::vector<double>& values, const Weights& weights,
                              std::size_t index) {
    check_shape(values, weights);
    if (index >= values.size()) throw std::runtime_error("lisa_classification: index");
    const Moments stats = moments(values);
    double lagged = 0.0;
    for (std::size_t j = 0; j < values.size(); ++j)
        lagged += weights[index][j] * (values[j] - stats.mean);
    const bool centre_high = values[index] > stats.mean;
    const bool neighbours_high = lagged > 0.0;
    if (centre_high && neighbours_high) return LisaClass::HighHigh;
    if (!centre_high && !neighbours_high) return LisaClass::LowLow;
    return centre_high ? LisaClass::HighLow : LisaClass::LowHigh;
}

double gaussian_kernel(double distance, double bandwidth) {
    if (!(bandwidth > 0.0)) throw std::runtime_error("gaussian_kernel: bandwidth");
    return std::exp(-0.5 * distance * distance / (bandwidth * bandwidth));
}

std::vector<double> gwr_coefficients(const Matrix& predictors, const std::vector<double>& response,
                                     const std::vector<double>& distances, double bandwidth) {
    if (predictors.empty() || predictors.size() != response.size() ||
        predictors.size() != distances.size())
        throw std::runtime_error("gwr_coefficients: size mismatch");
    const std::size_t parameters = predictors[0].size();
    Matrix normal(parameters, std::vector<double>(parameters, 0.0));
    std::vector<double> rhs(parameters, 0.0);
    for (std::size_t row = 0; row < predictors.size(); ++row) {
        if (predictors[row].size() != parameters)
            throw std::runtime_error("gwr_coefficients: ragged predictors");
        const double weight = gaussian_kernel(distances[row], bandwidth);
        for (std::size_t i = 0; i < parameters; ++i) {
            rhs[i] += weight * predictors[row][i] * response[row];
            for (std::size_t j = 0; j < parameters; ++j)
                normal[i][j] += weight * predictors[row][i] * predictors[row][j];
        }
    }
    std::vector<double> coefficients;
    if (!solve_linear_system(normal, rhs, coefficients))
        throw std::runtime_error("gwr_coefficients: singular local design");
    return coefficients;
}

double matern32_covariance(double distance, double variance, double length_scale) {
    if (!(variance >= 0.0) || !(length_scale > 0.0))
        throw std::runtime_error("matern32_covariance: invalid parameters");
    const double scaled = std::sqrt(3.0) * distance / length_scale;
    return variance * (1.0 + scaled) * std::exp(-scaled);
}

double probit_category_probability(double latent_mean, double lower_threshold,
                                   double upper_threshold, double sigma) {
    if (!(sigma > 0.0)) throw std::runtime_error("probit_category_probability: sigma");
    const auto normal_cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
    const double lower =
        std::isinf(lower_threshold) ? (lower_threshold < 0.0 ? 0.0 : 1.0)
                                    : normal_cdf((lower_threshold - latent_mean) / sigma);
    const double upper =
        std::isinf(upper_threshold) ? (upper_threshold < 0.0 ? 0.0 : 1.0)
                                    : normal_cdf((upper_threshold - latent_mean) / sigma);
    return std::max(0.0, upper - lower);
}

double hotelling_functional_scan(int inside_count, int outside_count,
                                 const std::vector<double>& inside_mean,
                                 const std::vector<double>& outside_mean,
                                 const Matrix& pooled_covariance) {
    if (inside_count <= 0 || outside_count <= 0)
        throw std::runtime_error("hotelling_functional_scan: counts must be positive");
    const std::size_t dimensions = inside_mean.size();
    if (outside_mean.size() != dimensions || pooled_covariance.size() != dimensions)
        throw std::runtime_error("hotelling_functional_scan: size mismatch");
    std::vector<double> difference(dimensions, 0.0);
    for (std::size_t i = 0; i < dimensions; ++i) difference[i] = inside_mean[i] - outside_mean[i];
    std::vector<double> solved;
    if (!solve_linear_system(pooled_covariance, difference, solved))
        throw std::runtime_error("hotelling_functional_scan: singular covariance");
    double quadratic = 0.0;
    for (std::size_t i = 0; i < dimensions; ++i) quadratic += difference[i] * solved[i];
    return static_cast<double>(inside_count) * static_cast<double>(outside_count) /
           static_cast<double>(inside_count + outside_count) * quadratic;
}

double crisp_cusum_statistic(double inside_sum, double outside_sum, std::size_t inside_count,
                             std::size_t outside_count) {
    if (inside_count == 0 || outside_count == 0) return 0.0;
    const double n = static_cast<double>(inside_count + outside_count);
    const double contrast = inside_sum / static_cast<double>(inside_count) -
                            outside_sum / static_cast<double>(outside_count);
    return std::sqrt(static_cast<double>(inside_count) *
                     static_cast<double>(outside_count) / n) *
           contrast;
}

std::size_t crisp_argmax(const std::vector<double>& statistics) {
    if (statistics.empty()) throw std::runtime_error("crisp_argmax: empty");
    std::size_t best = 0;
    for (std::size_t i = 1; i < statistics.size(); ++i)
        if (std::fabs(statistics[i]) > std::fabs(statistics[best])) best = i;
    return best;
}

double crisp_criterion(const std::vector<double>& values, const std::vector<int>& region_mask,
                       double penalty, double lambda) {
    if (values.size() != region_mask.size())
        throw std::runtime_error("crisp_criterion: size mismatch");
    double inside_sum = 0.0, outside_sum = 0.0;
    int inside_count = 0, outside_count = 0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (region_mask[i] != 0) {
            inside_sum += values[i];
            ++inside_count;
        } else {
            outside_sum += values[i];
            ++outside_count;
        }
    }
    const double inside_mean = inside_count > 0 ? inside_sum / inside_count : 0.0;
    const double outside_mean = outside_count > 0 ? outside_sum / outside_count : 0.0;
    double criterion = lambda * penalty;
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double mean = region_mask[i] != 0 ? inside_mean : outside_mean;
        criterion += (values[i] - mean) * (values[i] - mean);
    }
    return criterion;
}

double log_gaussian_cox_intensity(double intercept, const std::vector<double>& covariates,
                                  const std::vector<double>& coefficients,
                                  double spatial_effect) {
    if (covariates.size() != coefficients.size())
        throw std::runtime_error("log_gaussian_cox_intensity: size mismatch");
    double linear = intercept + spatial_effect;
    for (std::size_t i = 0; i < covariates.size(); ++i) linear += covariates[i] * coefficients[i];
    return linear;
}

double deviance_information_criterion(double mean_deviance, double deviance_at_mean) {
    const double effective_parameters = mean_deviance - deviance_at_mean;
    return effective_parameters + mean_deviance;
}

}  // namespace climate
