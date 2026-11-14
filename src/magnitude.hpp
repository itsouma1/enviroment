#ifndef CLIMATE_MAGNITUDE_HPP
#define CLIMATE_MAGNITUDE_HPP

#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Polynomial trend derivatives
// ---------------------------------------------------------------------------
double quadratic_trend_level(double intercept, double linear, double quadratic, double time);

double quadratic_trend_rate(double linear, double quadratic, double time);

double quadratic_trend_acceleration(double quadratic);

double linearithmic_trend_level(double intercept, double coefficient, double time);

double linearithmic_trend_rate(double coefficient, double time);

double linearithmic_trend_acceleration(double coefficient, double time);

// ---------------------------------------------------------------------------
// Abrupt and integrated change
// ---------------------------------------------------------------------------
double step_magnitude(double level_after, double level_before);

double integrated_change(const std::vector<double>& values, double time_step);

// Shift in return level between two non-stationary fits.
double return_level_shift(double future_return_level, double past_return_level);

// Newey-West long-run variance of a residual series with a Bartlett kernel.
double newey_west_variance(const std::vector<double>& residuals, int lag);

}  // namespace climate

#endif
