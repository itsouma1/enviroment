#ifndef CLIMATE_ASSIM_HPP
#define CLIMATE_ASSIM_HPP

#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Hybrid 4D-Var - EnKF gain
// ---------------------------------------------------------------------------
// x_hyb = gamma x_EnKF + (1 - gamma) x_EnVar.
std::vector<double> hybrid_analysis(const std::vector<double>& enkf_analysis,
                                     const std::vector<double>& enviar_analysis, double gamma);

// d x_hyb / d gamma = x_EnKF - x_EnVar.
std::vector<double> hybrid_analysis_sensitivity(const std::vector<double>& enkf_analysis,
                                                const std::vector<double>& enviar_analysis);

// K_hyb = gamma K_EnKF + (1 - gamma) K_EnVar, applied entry-wise.
std::vector<std::vector<double>> hybrid_gain(
    const std::vector<std::vector<double>>& enkf_gain,
    const std::vector<std::vector<double>>& enviar_gain, double gamma);

// ---------------------------------------------------------------------------
// Ensemble Adjustment Kalman Filter (Anderson 2001)
// ---------------------------------------------------------------------------
struct EnsembleAnalysis {
    std::vector<std::vector<double>> members;  // members x state
    std::vector<double> mean;
    std::vector<std::vector<double>> covariance;
};

// Deterministic square-root update:
//   Pa = [Pb^-1 + H^T R^-1 H]^-1
//   xa = Pa [Pb^-1 xb + H^T R^-1 y]
//   X_i^a = La^T Lb^-1 (X_i^b - xb) + xa, with Pb = Lb Lb^T, Pa = La La^T.
EnsembleAnalysis ensemble_adjustment_filter(
    const std::vector<std::vector<double>>& members, const std::vector<double>& observation,
    const std::vector<std::vector<double>>& observation_operator,
    const std::vector<std::vector<double>>& observation_covariance);

// ---------------------------------------------------------------------------
// Feed-forward network (Hybrid-JSBACH4 architecture)
// ---------------------------------------------------------------------------
// The published network is a 3 x 128 ReLU FNN with three outputs. The paper
// does not release a universal weight set, so this provides inference and an
// exact Jacobian for weights supplied by the caller.
struct FeedForwardNetwork {
    // weights[l][o][i], biases[l][o]; activation[l] is 0 (linear) or 1 (ReLU).
    std::vector<std::vector<std::vector<double>>> weights;
    std::vector<std::vector<double>> biases;
    std::vector<int> activation;
};

std::vector<double> network_forward(const FeedForwardNetwork& network,
                                    const std::vector<double>& input);

std::vector<std::vector<double>> network_jacobian(const FeedForwardNetwork& network,
                                                  const std::vector<double>& input);

// Ridge regression: beta = (X^T X + lambda I)^-1 X^T y.
std::vector<double> ridge_regression(const std::vector<std::vector<double>>& design,
                                     const std::vector<double>& response, double lambda);

// Prediction w^T x + b.
double linear_prediction(const std::vector<double>& coefficients,
                         const std::vector<double>& predictors, double intercept = 0.0);

}  // namespace climate

#endif
