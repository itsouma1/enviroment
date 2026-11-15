#ifndef CLIMATE_RESTRAND_HPP
#define CLIMATE_RESTRAND_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

std::vector<double> restrand_residuals(const std::vector<double>& response,
                                       const Matrix& predictors);

double residual_trend(const std::vector<double>& residuals);

struct SegmentedTrend {
    double intercept = 0.0;
    double slope = 0.0;
    double level_shift = 0.0;
    double slope_change = 0.0;
};

SegmentedTrend segmented_trend(const std::vector<double>& times,
                               const std::vector<double>& values, double breakpoint);

std::vector<double> sdid_unit_weights(const Matrix& control_outcomes,
                                      const std::vector<double>& treated_mean,
                                      double ridge);

double sdid_effect(const Matrix& outcomes, const std::vector<int>& treated_units,
                   std::size_t treatment_time, const std::vector<double>& unit_weights,
                   const std::vector<double>& time_weights);

double placebo_p_value(std::vector<double> placebo_effects, double observed_effect);

std::vector<double> stability_selection(const Matrix& predictors,
                                        const std::vector<double>& response,
                                        std::size_t subsamples, double threshold,
                                        unsigned seed);

std::size_t knockoff_threshold(std::vector<double> statistics, double q);

// 0-1 knapsack sensor selection under an integer budget.
std::vector<std::size_t> knapsack_select(const std::vector<double>& weights,
                                         const std::vector<int>& costs, int budget);

double theil_sen_slope(const std::vector<double>& times, const std::vector<double>& values);

double pettitt_statistic(const std::vector<double>& values, std::size_t& change_point);

std::vector<std::size_t> benjamini_hochberg(const std::vector<double>& p_values, double q);

}  // namespace climate

#endif
