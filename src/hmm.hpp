#ifndef CLIMATE_HMM_HPP
#define CLIMATE_HMM_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// Gaussian-emission hidden Markov model.
struct GaussianHmm {
    std::vector<double> initial;   // pi, length K
    Matrix transition;             // A, K x K
    std::vector<double> means;     // mu, length K
    std::vector<double> variances; // sigma^2, length K
};

struct HmmPosterior {
    Matrix gamma;              // T x K posterior state probabilities
    std::vector<Matrix> xi;    // T-1 transition posteriors, each K x K
    double log_likelihood = 0.0;
};

HmmPosterior hmm_forward_backward(const GaussianHmm& model,
                                  const std::vector<double>& observations);

// Baum-Welch EM re-estimation of pi, A, mu and sigma^2.
GaussianHmm baum_welch(const GaussianHmm& model, const std::vector<double>& observations,
                       int max_iterations = 100, double tolerance = 1e-6);

}  // namespace climate

#endif
