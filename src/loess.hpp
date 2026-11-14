#pragma once

#include <vector>

namespace climate {

// Local polynomial regression (LOESS). The estimate at each index is the value
// at the centre of a weighted local polynomial fit over the "span" nearest
// points, using tricube weights. robustness may be empty (all weights 1) or
// contain one non-negative robustness weight per observation, as used by the
// STL outer loop.
std::vector<double> loess_smooth(const std::vector<double>& values, int span, int degree,
                                 const std::vector<double>& robustness = {});

// Centred moving average. Even window lengths use half weights on the two end
// points so the filter stays symmetric. Missing values must be filled first.
std::vector<double> moving_average(const std::vector<double>& values, int window);

// Median absolute deviation (unscaled median of |x - median(x)|).
double median_absolute_deviation(const std::vector<double>& values);

}  // namespace climate
