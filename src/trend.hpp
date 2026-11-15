#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace climate {

struct MKSettings {
    bool hamed_rao = true;               // autocorrelation-corrected variance
    bool trend_free_prewhitening = false; // remove AR(1) from the detrended series first
    double alpha = 0.05;
};

struct MKResult {
    std::size_t n = 0;
    double s = 0.0;
    double variance = 0.0;            // tie-corrected variance of S
    double variance_corrected = 0.0;  // variance after Hamed-Rao
    double correction_factor = 1.0;   // variance_corrected / variance
    double z = 0.0;
    double p_value = 1.0;
    double tau = 0.0;
    double critical_value = 1.96;
    bool significant = false;
    std::string direction = "no trend";
};

struct SenResult {
    double slope = 0.0;
    double intercept = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    double alpha = 0.05;
    std::size_t pairs = 0;
    bool computed = false;
};

// Mann-Kendall trend test with tie correction and optional Hamed-Rao variance
// correction for serial correlation.
MKResult mann_kendall(const std::vector<double>& values, const MKSettings& settings = {});

// Sen's non-parametric slope with the rank-based confidence interval.
SenResult sens_slope(const std::vector<double>& time, const std::vector<double>& values,
                     const SenResult* unused = nullptr, double alpha = 0.05);

// Trend-free prewhitening (Yue et al., 2002): estimate the Sen slope, remove it,
// fit and remove lag-1 autocorrelation, then restore the trend.
std::vector<double> trend_free_prewhitening(const std::vector<double>& time,
                                            const std::vector<double>& values);

// Seasonal (within-block) Mann-Kendall test: S is accumulated per season across
// years, which removes the seasonal cycle before testing.
MKResult seasonal_mann_kendall(const std::vector<double>& values, int period,
                               const MKSettings& settings = {});

// A single trend estimator in the accuracy-weighted ensemble. Different
// estimators answer the same question with different assumptions; showing all
// of them makes the choice auditable instead of hidden.
struct TrendEstimator {
    std::string name;
    std::string formula;
    double slope = 0.0;
    double intercept = 0.0;
    double rmse = 0.0;
    double p_value = 1.0;
    double score = 0.0;
    bool significant = false;
    bool computed = false;
    bool selected = false;
};

struct TrendEnsemble {
    std::vector<TrendEstimator> members;
    std::size_t selected = 0;
    double agreement = 0.0;  // fraction of members sharing the selected sign
    double robust_slope = 0.0;
    double robust_p_value = 1.0;
    std::string consensus = "no trend";
};

// Runs Sen, OLS, Kendall, Spearman, Cox-Stuart and endpoint-difference
// estimators. Scores each on residual accuracy, significance and sign
// agreement; the highest score is primary and the rest are fallbacks.
TrendEnsemble trend_ensemble(const std::vector<double>& time, const std::vector<double>& values,
                             double alpha = 0.05);

}  // namespace climate
