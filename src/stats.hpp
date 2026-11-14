#pragma once

#include <cstddef>
#include <vector>

namespace climate {

struct Summary {
    std::size_t n = 0;
    double mean = 0.0;
    double sd = 0.0;
    double minimum = 0.0;
    double maximum = 0.0;
    double median = 0.0;
    double q25 = 0.0;
    double q75 = 0.0;
};

Summary summarize(const std::vector<double>& values);

double mean_value(const std::vector<double>& values);
double variance_value(const std::vector<double>& values);  // sample variance (n-1)
double stdev_value(const std::vector<double>& values);
double median_value(std::vector<double> values);           // takes by value, sorts in place
double quantile_value(std::vector<double> values, double probability);

double normal_cdf(double z);
double normal_quantile(double probability);

// Sample autocorrelation of the series at the given lag, returned in [-1, 1].
double autocorrelation(const std::vector<double>& values, int lag);

// Ordinary least squares slope of y against x (used for simple diagnostics).
double ols_slope(const std::vector<double>& x, const std::vector<double>& y);

}  // namespace climate
