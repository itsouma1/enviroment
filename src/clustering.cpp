#include "clustering.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double squared_distance(const std::vector<double>& a, const std::vector<double>& b) {
    require(a.size() == b.size(), "squared_distance: dimension mismatch");
    double sum = 0.0;
    for (std::size_t d = 0; d < a.size(); ++d) {
        const double difference = a[d] - b[d];
        sum += difference * difference;
    }
    return sum;
}

double ols_weight(std::size_t index, std::size_t n) {
    const double mean = 0.5 * (static_cast<double>(n) - 1.0);
    double variance = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dt = static_cast<double>(i) - mean;
        variance += dt * dt;
    }
    return (static_cast<double>(index) - mean) / variance;
}

}  // namespace

double fuzzy_objective(const Matrix& points, const Matrix& membership,
                       const Matrix& centroids, double fuzzifier) {
    require(points.size() == membership.size() && !centroids.empty(),
            "fuzzy_objective: size");
    double objective = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        for (std::size_t k = 0; k < centroids.size(); ++k) {
            const double weight = std::pow(membership[i][k], fuzzifier);
            objective += weight * squared_distance(points[i], centroids[k]);
        }
    }
    return objective;
}

FuzzyCMeansResult fuzzy_c_means(const Matrix& points, std::size_t clusters,
                                double fuzzifier, double tolerance, int max_iterations) {
    require(!points.empty() && !points[0].empty() && clusters > 0,
            "fuzzy_c_means: inputs");
    require(fuzzifier > 1.0 && tolerance > 0.0 && max_iterations > 0,
            "fuzzy_c_means: parameters");
    const std::size_t n = points.size(), dimensions = points[0].size();
    FuzzyCMeansResult result;
    result.centroids.assign(clusters, std::vector<double>(dimensions, 0.0));
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return points[a][0] < points[b][0]; });
    for (std::size_t k = 0; k < clusters; ++k) {
        const std::size_t pick = std::min(n - 1, (k * n) / clusters);
        result.centroids[k] = points[order[pick]];
    }
    result.membership.assign(n, std::vector<double>(clusters, 0.0));
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        for (std::size_t i = 0; i < n; ++i) {
            bool singular = false;
            std::size_t singular_cluster = 0;
            for (std::size_t k = 0; k < clusters; ++k)
                if (squared_distance(points[i], result.centroids[k]) < 1e-24) {
                    singular = true;
                    singular_cluster = k;
                    break;
                }
            if (singular) {
                for (std::size_t k = 0; k < clusters; ++k)
                    result.membership[i][k] = k == singular_cluster ? 1.0 : 0.0;
                continue;
            }
            double total = 0.0;
            for (std::size_t k = 0; k < clusters; ++k) {
                const double numerator = std::sqrt(squared_distance(points[i], result.centroids[k]));
                double sum = 0.0;
                for (std::size_t j = 0; j < clusters; ++j) {
                    const double denominator =
                        std::sqrt(squared_distance(points[i], result.centroids[j]));
                    sum += std::pow(numerator / denominator, 2.0 / (fuzzifier - 1.0));
                }
                result.membership[i][k] = 1.0 / sum;
                total += result.membership[i][k];
            }
            for (std::size_t k = 0; k < clusters; ++k) result.membership[i][k] /= total;
        }
        Matrix updated(clusters, std::vector<double>(dimensions, 0.0));
        for (std::size_t k = 0; k < clusters; ++k) {
            double weight_sum = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                const double weight = std::pow(result.membership[i][k], fuzzifier);
                weight_sum += weight;
                for (std::size_t d = 0; d < dimensions; ++d)
                    updated[k][d] += weight * points[i][d];
            }
            require(weight_sum > 0.0, "fuzzy_c_means: empty cluster");
            for (double& value : updated[k]) value /= weight_sum;
        }
        double change = 0.0;
        for (std::size_t k = 0; k < clusters; ++k)
            for (std::size_t d = 0; d < dimensions; ++d)
                change = std::max(change, std::fabs(updated[k][d] - result.centroids[k][d]));
        result.centroids = updated;
        result.iterations = iteration + 1;
        if (change < tolerance) break;
    }
    result.objective =
        fuzzy_objective(points, result.membership, result.centroids, fuzzifier);
    return result;
}

std::vector<double> fuzzy_membership_sensitivity(const std::vector<double>& point,
                                                 const Matrix& centroids,
                                                 const std::vector<double>& membership) {
    require(!centroids.empty() && centroids.size() == membership.size() &&
                point.size() == centroids[0].size(),
            "fuzzy_membership_sensitivity: size");
    const std::size_t clusters = centroids.size();
    std::vector<double> weighted_mean(point.size(), 0.0);
    std::vector<double> distances(clusters, 0.0);
    for (std::size_t k = 0; k < clusters; ++k) {
        distances[k] = squared_distance(point, centroids[k]);
        require(distances[k] > 0.0, "fuzzy_membership_sensitivity: point on centroid");
        for (std::size_t d = 0; d < point.size(); ++d)
            weighted_mean[d] += membership[k] * (point[d] - centroids[k][d]) / distances[k];
    }
    std::vector<double> sensitivity(clusters, 0.0);
    for (std::size_t k = 0; k < clusters; ++k)
        sensitivity[k] = -membership[k] *
                         ((point[0] - centroids[k][0]) / distances[k] - weighted_mean[0]);
    return sensitivity;
}

double fuzzy_centroid_trend(const std::vector<double>& trends,
                            const std::vector<double>& memberships, double fuzzifier) {
    require(!trends.empty() && trends.size() == memberships.size(),
            "fuzzy_centroid_trend: size");
    double weighted = 0.0, total = 0.0;
    for (std::size_t i = 0; i < trends.size(); ++i) {
        const double weight = std::pow(memberships[i], fuzzifier);
        weighted += weight * trends[i];
        total += weight;
    }
    require(total > 0.0, "fuzzy_centroid_trend: zero membership");
    return weighted / total;
}

double attribute_distance(const std::vector<double>& first,
                          const std::vector<double>& second, double attribute_weight) {
    require(first.size() == second.size() && first.size() >= 3 && attribute_weight >= 0.0,
            "attribute_distance: inputs");
    double sum = 0.0;
    for (std::size_t d = 0; d + 1 < first.size(); ++d) {
        const double difference = first[d] - second[d];
        sum += difference * difference;
    }
    const double attribute = first.back() - second.back();
    return std::sqrt(sum + attribute_weight * attribute * attribute);
}

std::vector<int> dbscan(const Matrix& points, double epsilon, std::size_t min_points) {
    require(!points.empty() && epsilon > 0.0 && min_points > 0, "dbscan: inputs");
    const std::size_t n = points.size();
    std::vector<int> labels(n, -1);
    std::vector<int> visited(n, 0);
    int cluster = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (visited[i]) continue;
        visited[i] = 1;
        std::vector<std::size_t> neighbours;
        for (std::size_t j = 0; j < n; ++j)
            if (std::sqrt(squared_distance(points[i], points[j])) <= epsilon)
                neighbours.push_back(j);
        if (neighbours.size() < min_points) continue;
        labels[i] = cluster;
        for (std::size_t index = 0; index < neighbours.size(); ++index) {
            const std::size_t candidate = neighbours[index];
            if (!visited[candidate]) {
                visited[candidate] = 1;
                std::vector<std::size_t> candidate_neighbours;
                for (std::size_t j = 0; j < n; ++j)
                    if (std::sqrt(squared_distance(points[candidate], points[j])) <= epsilon)
                        candidate_neighbours.push_back(j);
                if (candidate_neighbours.size() >= min_points)
                    neighbours.insert(neighbours.end(), candidate_neighbours.begin(),
                                      candidate_neighbours.end());
            }
            if (labels[candidate] < 0) labels[candidate] = cluster;
        }
        ++cluster;
    }
    return labels;
}

double dbscan_critical_weight(double epsilon, double spatial_distance, double attribute_gap) {
    require(attribute_gap != 0.0, "dbscan_critical_weight: zero attribute gap");
    return (epsilon * epsilon - spatial_distance * spatial_distance) /
           (attribute_gap * attribute_gap);
}

std::vector<int> dbscan_attribute(const Matrix& points, double attribute_weight,
                                  double epsilon, std::size_t min_points) {
    require(!points.empty() && !points[0].empty() && attribute_weight >= 0.0,
            "dbscan_attribute: inputs");
    Matrix augmented(points.size(), std::vector<double>(points[0].size(), 0.0));
    for (std::size_t i = 0; i < points.size(); ++i)
        for (std::size_t d = 0; d < points[i].size(); ++d)
            augmented[i][d] = d + 1 == points[i].size()
                                  ? std::sqrt(attribute_weight) * points[i][d]
                                  : points[i][d];
    return dbscan(augmented, epsilon, min_points);
}

double cluster_mean_slope(const std::vector<double>& slopes, const std::vector<int>& labels,
                          int cluster) {
    require(slopes.size() == labels.size(), "cluster_mean_slope: size");
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < slopes.size(); ++i)
        if (labels[i] == cluster) {
            sum += slopes[i];
            ++count;
        }
    require(count > 0, "cluster_mean_slope: empty cluster");
    return sum / static_cast<double>(count);
}

double varimax_criterion(const Matrix& loadings) {
    require(!loadings.empty() && !loadings[0].empty(), "varimax_criterion: loadings");
    const std::size_t sites = loadings.size(), modes = loadings[0].size();
    double criterion = 0.0;
    for (std::size_t k = 0; k < modes; ++k) {
        double fourth = 0.0, second = 0.0;
        for (std::size_t s = 0; s < sites; ++s) {
            const double squared = loadings[s][k] * loadings[s][k];
            fourth += squared * squared;
            second += squared;
        }
        fourth /= static_cast<double>(sites);
        second /= static_cast<double>(sites);
        criterion += fourth - second * second;
    }
    return criterion;
}

std::vector<int> reof_partition(const Matrix& loadings, double threshold) {
    require(!loadings.empty() && !loadings[0].empty(), "reof_partition: loadings");
    std::vector<int> partition(loadings.size(), -1);
    for (std::size_t s = 0; s < loadings.size(); ++s) {
        double best = threshold;
        for (std::size_t k = 0; k < loadings[s].size(); ++k)
            if (std::fabs(loadings[s][k]) > best) {
                best = std::fabs(loadings[s][k]);
                partition[s] = static_cast<int>(k);
            }
    }
    return partition;
}

Matrix rpc_trend_sensitivity(const std::vector<double>& times, const Matrix& loadings,
                             const Matrix& rotation, std::size_t mode) {
    require(!times.empty() && !loadings.empty() && rotation.size() == loadings[0].size() &&
                mode < rotation.size(),
            "rpc_trend_sensitivity: size");
    const std::size_t sites = loadings.size();
    Matrix sensitivity(times.size(), std::vector<double>(sites, 0.0));
    for (std::size_t t = 0; t < times.size(); ++t) {
        const double weight = ols_weight(t, times.size());
        for (std::size_t s = 0; s < sites; ++s)
            for (std::size_t j = 0; j < rotation.size(); ++j)
                sensitivity[t][s] += weight * loadings[s][j] * rotation[j][mode];
    }
    return sensitivity;
}

double trend_sign_flip_probability(double trend, double standard_error, double alpha) {
    require(standard_error > 0.0 && alpha > 0.0 && alpha < 1.0,
            "trend_sign_flip_probability: inputs");
    const double critical = normal_quantile(1.0 - 0.5 * alpha);
    const double ratio = trend / standard_error;
    return normal_cdf(critical - ratio) - normal_cdf(-critical - ratio);
}

std::vector<double> categorical_assignment_probabilities(
    const std::vector<double>& log_likelihoods) {
    require(!log_likelihoods.empty(), "categorical_assignment_probabilities: empty");
    const double maximum = *std::max_element(log_likelihoods.begin(), log_likelihoods.end());
    std::vector<double> probabilities(log_likelihoods.size());
    double total = 0.0;
    for (std::size_t k = 0; k < log_likelihoods.size(); ++k) {
        probabilities[k] = std::exp(log_likelihoods[k] - maximum);
        total += probabilities[k];
    }
    for (double& value : probabilities) value /= total;
    return probabilities;
}

}  // namespace climate
