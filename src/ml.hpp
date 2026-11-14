#ifndef CLIMATE_ML_HPP
#define CLIMATE_ML_HPP

#include <vector>

namespace climate {

using Matrix = std::vector<std::vector<double>>;

// ---------------------------------------------------------------------------
// CART regression tree
// ---------------------------------------------------------------------------
struct RegressionTree {
    struct Node {
        int feature = -1;
        double threshold = 0.0;
        double value = 0.0;
        int left = -1;
        int right = -1;
    };

    std::vector<Node> nodes;
    std::vector<double> importance;

    void fit(const Matrix& predictors, const std::vector<double>& response,
             int max_depth = 8, int min_leaf = 2, int features_per_split = 0,
             unsigned seed = 0);

    double predict(const std::vector<double>& predictors) const;
};

// ---------------------------------------------------------------------------
// Random forest (deweathered-ML learner)
// ---------------------------------------------------------------------------
struct RandomForest {
    std::vector<RegressionTree> trees;
    std::vector<double> importance;

    void fit(const Matrix& predictors, const std::vector<double>& response, int tree_count = 100,
             int max_depth = 8, int min_leaf = 2, int features_per_split = 0,
             unsigned seed = 0);

    double predict(const std::vector<double>& predictors) const;

    // Apply the trained forest to a fixed meteorological state.
    double deweathered_prediction(const std::vector<double>& fixed_state) const;
};

// ---------------------------------------------------------------------------
// Causal forests
// ---------------------------------------------------------------------------
// Honest split criterion: n1 n2 / (n1 + n2)^2 * (tau1 - tau2)^2.
double causal_forest_split_criterion(double left_count, double right_count, double left_effect,
                                     double right_effect);

// CATE estimand: mean(Y | treatment = 1) - mean(Y | treatment = 0).
double cate_estimate(const std::vector<double>& treated_outcomes,
                     const std::vector<double>& control_outcomes);

}  // namespace climate

#endif
