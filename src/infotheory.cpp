#include "infotheory.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::size_t bin_index(double value, double minimum, double maximum, std::size_t bins) {
    if (maximum <= minimum) return 0;
    const double scaled = (value - minimum) / (maximum - minimum) * static_cast<double>(bins);
    const long index = static_cast<long>(std::floor(scaled));
    return static_cast<std::size_t>(std::max(0L, std::min<long>(index, static_cast<long>(bins) - 1)));
}

std::vector<std::size_t> discretize(const std::vector<double>& values, std::size_t bins) {
    require(!values.empty() && bins > 0, "discretize: inputs");
    const auto range = std::minmax_element(values.begin(), values.end());
    std::vector<std::size_t> result(values.size(), 0);
    for (std::size_t i = 0; i < values.size(); ++i)
        result[i] = bin_index(values[i], *range.first, *range.second, bins);
    return result;
}

double entropy_from_counts(const std::map<std::size_t, std::size_t>& counts,
                           std::size_t samples) {
    double entropy = 0.0;
    for (const auto& entry : counts) {
        const double probability = static_cast<double>(entry.second) /
                                   static_cast<double>(samples);
        entropy -= probability * std::log(probability);
    }
    return entropy;
}

double plugin_entropy(const std::vector<std::size_t>& values) {
    std::map<std::size_t, std::size_t> counts;
    for (std::size_t value : values) ++counts[value];
    return entropy_from_counts(counts, values.size());
}

}  // namespace

double digamma(double value) {
    require(value > 0.0, "digamma: positive argument");
    double result = 0.0;
    while (value < 6.0) {
        result -= 1.0 / value;
        value += 1.0;
    }
    const double inverse = 1.0 / value;
    const double inverse2 = inverse * inverse;
    return result + std::log(value) - 0.5 * inverse -
           inverse2 * (1.0 / 12.0 - inverse2 * (1.0 / 120.0 - inverse2 / 252.0));
}

double mutual_information_plugin(const std::vector<double>& x, const std::vector<double>& y,
                                 std::size_t bins) {
    require(x.size() == y.size() && x.size() >= 2 && bins > 0,
            "mutual_information_plugin: inputs");
    const std::vector<std::size_t> xb = discretize(x, bins);
    const std::vector<std::size_t> yb = discretize(y, bins);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> joint;
    std::map<std::size_t, std::size_t> xcounts, ycounts;
    for (std::size_t i = 0; i < xb.size(); ++i) {
        ++joint[{xb[i], yb[i]}];
        ++xcounts[xb[i]];
        ++ycounts[yb[i]];
    }
    const double n = static_cast<double>(xb.size());
    double information = 0.0;
    for (const auto& entry : joint) {
        const double probability = static_cast<double>(entry.second) / n;
        const double marginal_product =
            static_cast<double>(xcounts[entry.first.first]) / n *
            static_cast<double>(ycounts[entry.first.second]) / n;
        information += probability * std::log(probability / marginal_product);
    }
    return information;
}

double miller_madow_correction(std::size_t x_bins, std::size_t y_bins, std::size_t samples) {
    require(samples > 0, "miller_madow_correction: samples");
    return static_cast<double>((x_bins - 1) * (y_bins - 1)) /
           (2.0 * static_cast<double>(samples) * std::log(2.0));
}

double partial_mutual_information_plugin(const std::vector<double>& x,
                                         const std::vector<double>& y,
                                         const Matrix& conditioning, std::size_t bins) {
    require(x.size() == y.size() && conditioning.size() == x.size() && bins > 0,
            "partial_mutual_information_plugin: size");
    if (conditioning.empty() || conditioning[0].empty())
        return mutual_information_plugin(x, y, bins);
    std::vector<std::vector<std::size_t>> strata(conditioning[0].size());
    for (std::size_t column = 0; column < conditioning[0].size(); ++column) {
        std::vector<double> values(conditioning.size(), 0.0);
        for (std::size_t i = 0; i < conditioning.size(); ++i)
            values[i] = conditioning[i][column];
        strata[column] = discretize(values, bins);
    }
    std::map<std::vector<std::size_t>, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::vector<std::size_t> key(strata.size(), 0);
        for (std::size_t column = 0; column < strata.size(); ++column)
            key[column] = strata[column][i];
        groups[key].push_back(i);
    }
    const double n = static_cast<double>(x.size());
    double information = 0.0;
    for (const auto& group : groups) {
        if (group.second.size() < 2) continue;
        std::vector<double> group_x(group.second.size()), group_y(group.second.size());
        for (std::size_t i = 0; i < group.second.size(); ++i) {
            group_x[i] = x[group.second[i]];
            group_y[i] = y[group.second[i]];
        }
        information += static_cast<double>(group.second.size()) / n *
                       mutual_information_plugin(group_x, group_y, bins);
    }
    return information;
}

double ksg_mutual_information(const std::vector<double>& x, const std::vector<double>& y,
                              std::size_t neighbours) {
    require(x.size() == y.size() && x.size() > neighbours && neighbours > 0,
            "ksg_mutual_information: inputs");
    const std::size_t n = x.size();
    double average = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<double> distances(n, 0.0);
        for (std::size_t j = 0; j < n; ++j)
            distances[j] = std::max(std::fabs(x[i] - x[j]), std::fabs(y[i] - y[j]));
        std::sort(distances.begin(), distances.end());
        const double epsilon = distances[neighbours];
        std::size_t nx = 0, ny = 0;
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i && std::fabs(x[i] - x[j]) < epsilon) ++nx;
            if (j != i && std::fabs(y[i] - y[j]) < epsilon) ++ny;
        }
        average += digamma(static_cast<double>(nx) + 1.0) +
                   digamma(static_cast<double>(ny) + 1.0);
    }
    average /= static_cast<double>(n);
    return digamma(static_cast<double>(neighbours)) + digamma(static_cast<double>(n)) -
           average;
}

double chi_square_quantile(double probability, double degrees_of_freedom) {
    require(probability > 0.0 && probability < 1.0 && degrees_of_freedom > 0.0,
            "chi_square_quantile: inputs");
    const double z = normal_quantile(0.5 + 0.5 * probability);
    if (degrees_of_freedom == 1.0) return z * z;
    if (degrees_of_freedom == 2.0) return -2.0 * std::log(1.0 - probability);
    const double term = 1.0 - 2.0 / (9.0 * degrees_of_freedom) +
                        z * std::sqrt(2.0 / (9.0 * degrees_of_freedom));
    return degrees_of_freedom * term * term * term;
}

double pmi_threshold(double alpha, std::size_t conditioning_states, std::size_t x_bins,
                     std::size_t y_bins, std::size_t samples) {
    require(alpha > 0.0 && alpha < 1.0 && samples > 0,
            "pmi_threshold: inputs");
    const double degrees = static_cast<double>(conditioning_states) *
                           static_cast<double>(x_bins > 0 ? x_bins - 1 : 1) *
                           static_cast<double>(y_bins > 0 ? y_bins - 1 : 1);
    return chi_square_quantile(1.0 - alpha, degrees) /
           (2.0 * static_cast<double>(samples) * std::log(2.0));
}

double optimal_bin_width(std::size_t samples, std::size_t dimensions, double c1, double c2) {
    require(samples > 0 && dimensions > 0 && c1 > 0.0 && c2 > 0.0,
            "optimal_bin_width: inputs");
    return std::pow(static_cast<double>(dimensions) * c2 /
                        (2.0 * c1 * static_cast<double>(samples)),
                    1.0 / (static_cast<double>(dimensions) + 2.0));
}

double mutual_information_influence(double log_density_ratio, double mutual_information) {
    return log_density_ratio - mutual_information;
}

double mutual_information_influence_variance(const std::vector<double>& log_density_ratios) {
    require(log_density_ratios.size() >= 2, "mi_influence_variance: size");
    double mean = 0.0;
    for (double value : log_density_ratios) mean += value;
    mean /= static_cast<double>(log_density_ratios.size());
    double variance = 0.0;
    for (double value : log_density_ratios) variance += (value - mean) * (value - mean);
    return variance / (static_cast<double>(log_density_ratios.size()) *
                       static_cast<double>(log_density_ratios.size() - 1));
}

double gaussian_partial_mutual_information(double delta_r_squared) {
    require(delta_r_squared >= 0.0 && delta_r_squared < 1.0,
            "gaussian_partial_mutual_information: range");
    return -0.5 * std::log(1.0 - delta_r_squared);
}

double information_redundancy(const std::vector<double>& marginal_information,
                              double joint_information) {
    double sum = 0.0;
    for (double value : marginal_information) sum += value;
    return sum - joint_information;
}

std::vector<std::size_t> pmi_forward_selection(const std::vector<double>& response,
                                               const Matrix& predictors, double threshold,
                                               std::size_t bins) {
    require(!predictors.empty() && predictors.size() == response.size() && bins > 0,
            "pmi_forward_selection: size");
    const std::size_t variables = predictors[0].size();
    std::vector<std::size_t> selected;
    std::vector<double> scores(variables, 0.0);
    while (true) {
        double best = threshold;
        int best_variable = -1;
        for (std::size_t j = 0; j < variables; ++j) {
            if (std::find(selected.begin(), selected.end(), j) != selected.end()) continue;
            std::vector<double> x(predictors.size(), 0.0);
            for (std::size_t i = 0; i < predictors.size(); ++i) x[i] = predictors[i][j];
            Matrix conditioning(predictors.size(), std::vector<double>(selected.size(), 0.0));
            for (std::size_t s = 0; s < selected.size(); ++s)
                for (std::size_t i = 0; i < predictors.size(); ++i)
                    conditioning[i][s] = predictors[i][selected[s]];
            const double pmi =
                partial_mutual_information_plugin(x, response, conditioning, bins);
            scores[j] = pmi;
            if (pmi > best) {
                best = pmi;
                best_variable = static_cast<int>(j);
            }
        }
        if (best_variable < 0) break;
        selected.push_back(static_cast<std::size_t>(best_variable));
    }
    return selected;
}

}  // namespace climate
