#include "stats.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "common.hpp"

namespace climate {

double mean_value(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    double sum = 0.0;
    for (double v : values) sum += v;
    return sum / static_cast<double>(values.size());
}

double variance_value(const std::vector<double>& values) {
    const std::size_t n = values.size();
    if (n < 2) return 0.0;
    const double m = mean_value(values);
    double sum = 0.0;
    for (double v : values) {
        const double d = v - m;
        sum += d * d;
    }
    return sum / static_cast<double>(n - 1);
}

double stdev_value(const std::vector<double>& values) { return std::sqrt(variance_value(values)); }

double median_value(std::vector<double> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if (n % 2 == 1) return values[n / 2];
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

double quantile_value(std::vector<double> values, double probability) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const double p = std::min(1.0, std::max(0.0, probability));
    const double position = p * static_cast<double>(values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    if (lower == upper) return values[lower];
    const double fraction = position - static_cast<double>(lower);
    return values[lower] * (1.0 - fraction) + values[upper] * fraction;
}

Summary summarize(const std::vector<double>& values) {
    Summary result;
    result.n = values.size();
    if (values.empty()) return result;
    result.minimum = values[0];
    result.maximum = values[0];
    for (double v : values) {
        if (v < result.minimum) result.minimum = v;
        if (v > result.maximum) result.maximum = v;
    }
    result.mean = mean_value(values);
    result.sd = stdev_value(values);
    result.median = median_value(values);
    result.q25 = quantile_value(values, 0.25);
    result.q75 = quantile_value(values, 0.75);
    return result;
}

double normal_cdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

// Acklam's rational approximation to the inverse normal CDF, refined with one
// Newton step so the result is accurate to close to machine precision.
double normal_quantile(double probability) {
    static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285104469687e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                                2.445134137142996e+00, 3.754408661907416e+00};

    const double p_low = 0.02425;
    const double p_high = 1.0 - p_low;
    double p = std::min(1.0 - 1e-15, std::max(1e-15, probability));
    double x = 0.0;

    if (p < p_low) {
        const double q = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else if (p <= p_high) {
        const double q = p - 0.5;
        const double r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    } else {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }

    // One Halley refinement step against the exact CDF.
    const double e = normal_cdf(x) - p;
    const double u = e * std::sqrt(2.0 * PI) * std::exp(0.5 * x * x);
    x = x - u / (1.0 + 0.5 * x * u);
    return x;
}

double autocorrelation(const std::vector<double>& values, int lag) {
    const std::size_t n = values.size();
    if (lag < 0 || static_cast<std::size_t>(lag) >= n) return 0.0;
    const double m = mean_value(values);
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = values[i] - m;
        denominator += d * d;
    }
    if (denominator <= 0.0) return 0.0;
    for (std::size_t i = 0; i + static_cast<std::size_t>(lag) < n; ++i) {
        numerator += (values[i] - m) * (values[i + static_cast<std::size_t>(lag)] - m);
    }
    // The classic estimator divides by (n - lag) * variance, which keeps the
    // estimate unbiased for short series.
    const double lag_count = static_cast<double>(n - static_cast<std::size_t>(lag));
    return (numerator / lag_count) / (denominator / static_cast<double>(n));
}

double ols_slope(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = std::min(x.size(), y.size());
    if (n < 2) return 0.0;
    const double mx = mean_value(std::vector<double>(x.begin(), x.begin() + n));
    const double my = mean_value(std::vector<double>(y.begin(), y.begin() + n));
    double covariance = 0.0;
    double variance = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        covariance += (x[i] - mx) * (y[i] - my);
        variance += (x[i] - mx) * (x[i] - mx);
    }
    if (variance <= 0.0) return 0.0;
    return covariance / variance;
}

}  // namespace climate
