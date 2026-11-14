#include "ml.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double mean_of(const std::vector<double>& values, const std::vector<std::size_t>& indices) {
    if (indices.empty()) return 0.0;
    double total = 0.0;
    for (std::size_t index : indices) total += values[index];
    return total / static_cast<double>(indices.size());
}

}  // namespace

void RegressionTree::fit(const Matrix& predictors, const std::vector<double>& response,
                         int max_depth, int min_leaf, int features_per_split, unsigned seed) {
    require(!predictors.empty() && predictors.size() == response.size(),
            "RegressionTree::fit: size mismatch");
    const std::size_t feature_count = predictors[0].size();
    require(feature_count > 0, "RegressionTree::fit: no predictors");
    for (const auto& row : predictors)
        require(row.size() == feature_count, "RegressionTree::fit: ragged predictors");
    require(min_leaf >= 1, "RegressionTree::fit: min_leaf must be positive");

    nodes.clear();
    importance.assign(feature_count, 0.0);
    if (features_per_split <= 0)
        features_per_split = static_cast<int>(std::sqrt(static_cast<double>(feature_count)));
    features_per_split = std::min(features_per_split, static_cast<int>(feature_count));

    std::mt19937 generator(seed);
    std::vector<std::size_t> indices(predictors.size());
    std::iota(indices.begin(), indices.end(), 0);

    struct Builder {
        const Matrix& predictors;
        const std::vector<double>& response;
        std::vector<double>& importance;
        std::vector<Node>& nodes;
        int max_depth;
        int min_leaf;
        int features_per_split;
        std::mt19937& generator;

        int build(const std::vector<std::size_t>& rows, int depth) {
            const int node_index = static_cast<int>(nodes.size());
            nodes.push_back(Node{});
            nodes[node_index].value = mean_of(response, rows);
            if (depth >= max_depth || static_cast<int>(rows.size()) < 2 * min_leaf) {
                return node_index;
            }

            const std::size_t feature_count = predictors[0].size();
            std::vector<std::size_t> features(feature_count);
            std::iota(features.begin(), features.end(), 0);
            std::shuffle(features.begin(), features.end(), generator);
            if (features_per_split < static_cast<int>(feature_count))
                features.resize(static_cast<std::size_t>(features_per_split));

            double total_sum = 0.0, total_square = 0.0;
            for (std::size_t index : rows) {
                total_sum += response[index];
                total_square += response[index] * response[index];
            }
            const double parent_sse =
                total_square - total_sum * total_sum / static_cast<double>(rows.size());

            double best_gain = -1.0;
            int best_feature = -1;
            double best_threshold = 0.0;

            for (std::size_t feature : features) {
                std::vector<double> candidates;
                candidates.reserve(rows.size());
                for (std::size_t index : rows) candidates.push_back(predictors[index][feature]);
                std::sort(candidates.begin(), candidates.end());
                candidates.erase(std::unique(candidates.begin(), candidates.end()),
                                 candidates.end());
                if (candidates.size() < 2) continue;

                for (std::size_t split = 0; split + 1 < candidates.size(); ++split) {
                    const double threshold = 0.5 * (candidates[split] + candidates[split + 1]);
                    double left_sum = 0.0, left_square = 0.0, right_sum = 0.0, right_square = 0.0;
                    std::size_t left_count = 0, right_count = 0;
                    for (std::size_t index : rows) {
                        if (predictors[index][feature] <= threshold) {
                            left_sum += response[index];
                            left_square += response[index] * response[index];
                            ++left_count;
                        } else {
                            right_sum += response[index];
                            right_square += response[index] * response[index];
                            ++right_count;
                        }
                    }
                    if (static_cast<int>(left_count) < min_leaf ||
                        static_cast<int>(right_count) < min_leaf)
                        continue;
                    const double child_sse =
                        (left_square - left_sum * left_sum / static_cast<double>(left_count)) +
                        (right_square - right_sum * right_sum / static_cast<double>(right_count));
                    const double gain = parent_sse - child_sse;
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_feature = static_cast<int>(feature);
                        best_threshold = threshold;
                    }
                }
            }

            if (best_feature < 0 || best_gain <= 1e-15) return node_index;

            std::vector<std::size_t> left_rows, right_rows;
            for (std::size_t index : rows) {
                if (predictors[index][static_cast<std::size_t>(best_feature)] <= best_threshold)
                    left_rows.push_back(index);
                else
                    right_rows.push_back(index);
            }
            nodes[node_index].feature = best_feature;
            nodes[node_index].threshold = best_threshold;
            nodes[node_index].left = build(left_rows, depth + 1);
            nodes[node_index].right = build(right_rows, depth + 1);
            importance[static_cast<std::size_t>(best_feature)] += best_gain;
            return node_index;
        }
    };

    Builder builder{predictors, response, importance, nodes, max_depth, min_leaf,
                    features_per_split, generator};
    builder.build(indices, 0);
    const double total = std::accumulate(importance.begin(), importance.end(), 0.0);
    if (total > 0.0)
        for (double& value : importance) value /= total;
}

double RegressionTree::predict(const std::vector<double>& predictors) const {
    require(!nodes.empty(), "RegressionTree::predict: tree is empty");
    int index = 0;
    while (nodes[static_cast<std::size_t>(index)].feature >= 0) {
        const Node& node = nodes[static_cast<std::size_t>(index)];
        require(static_cast<std::size_t>(node.feature) < predictors.size(),
                "RegressionTree::predict: predictor too short");
        index = predictors[static_cast<std::size_t>(node.feature)] <= node.threshold ? node.left
                                                                                     : node.right;
    }
    return nodes[static_cast<std::size_t>(index)].value;
}

void RandomForest::fit(const Matrix& predictors, const std::vector<double>& response,
                       int tree_count, int max_depth, int min_leaf, int features_per_split,
                       unsigned seed) {
    require(!predictors.empty() && predictors.size() == response.size(),
            "RandomForest::fit: size mismatch");
    require(tree_count > 0, "RandomForest::fit: tree_count must be positive");
    trees.clear();
    const std::size_t sample_count = predictors.size();
    const std::size_t feature_count = predictors[0].size();
    importance.assign(feature_count, 0.0);

    std::mt19937 generator(seed);
    std::uniform_int_distribution<std::size_t> pick(0, sample_count - 1);
    for (int tree = 0; tree < tree_count; ++tree) {
        Matrix bootstrap(sample_count);
        std::vector<double> bootstrap_response(sample_count);
        for (std::size_t i = 0; i < sample_count; ++i) {
            const std::size_t source = pick(generator);
            bootstrap[i] = predictors[source];
            bootstrap_response[i] = response[source];
        }
        RegressionTree model;
        model.fit(bootstrap, bootstrap_response, max_depth, min_leaf, features_per_split,
                  generator());
        for (std::size_t i = 0; i < feature_count; ++i) importance[i] += model.importance[i];
        trees.push_back(std::move(model));
    }
    for (double& value : importance) value /= static_cast<double>(tree_count);
}

double RandomForest::predict(const std::vector<double>& predictors) const {
    require(!trees.empty(), "RandomForest::predict: forest is empty");
    double total = 0.0;
    for (const auto& tree : trees) total += tree.predict(predictors);
    return total / static_cast<double>(trees.size());
}

double RandomForest::deweathered_prediction(const std::vector<double>& fixed_state) const {
    return predict(fixed_state);
}

double causal_forest_split_criterion(double left_count, double right_count, double left_effect,
                                     double right_effect) {
    const double total = left_count + right_count;
    require(total > 0.0, "causal_forest_split_criterion: zero sample count");
    return left_count * right_count / (total * total) *
           (left_effect - right_effect) * (left_effect - right_effect);
}

double honest_emse_criterion(double treatment_effect_squared_sum, double splitting_count,
                             double estimation_count, double leaf_variance_sum) {
    require(splitting_count > 0.0 && estimation_count > 0.0,
            "honest_emse_criterion: non-positive sample count");
    return treatment_effect_squared_sum / splitting_count -
           (1.0 / splitting_count + 1.0 / estimation_count) * leaf_variance_sum;
}

double cate_estimate(const std::vector<double>& treated_outcomes,
                     const std::vector<double>& control_outcomes) {
    require(!treated_outcomes.empty() && !control_outcomes.empty(),
            "cate_estimate: empty group");
    double treated = 0.0, control = 0.0;
    for (double value : treated_outcomes) treated += value;
    for (double value : control_outcomes) control += value;
    return treated / static_cast<double>(treated_outcomes.size()) -
           control / static_cast<double>(control_outcomes.size());
}

}  // namespace climate
