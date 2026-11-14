#pragma once

#include <vector>

namespace climate {

struct STLParams {
    int period = 12;              // number of observations in one seasonal cycle
    int seasonal_window = 13;     // loess span for the cycle-subseries smoother
    int trend_window = 25;        // loess span for the trend smoother
    int low_pass_window = 13;     // loess span applied after the low-pass moving averages
    int seasonal_degree = 1;
    int trend_degree = 1;
    int inner_iterations = 2;
    int outer_iterations = 5;     // robustness iterations
};

struct STLResult {
    std::vector<double> seasonal;
    std::vector<double> trend;
    std::vector<double> remainder;
    std::vector<double> robustness;
};

// Seasonal-Trend decomposition using LOESS (Cleveland et al., 1990):
//   y = trend + seasonal + remainder
// The series must be regularly spaced and free of missing values (fill gaps
// first with interpolate_missing); a runtime_error is thrown otherwise.
STLResult stl_decompose(const std::vector<double>& values, int period, const STLParams& params);

// Convenience overload using default parameters for the given period.
STLResult stl_decompose(const std::vector<double>& values, int period);

}  // namespace climate
