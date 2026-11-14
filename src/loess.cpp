#include "loess.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "linalg.hpp"
#include "stats.hpp"

namespace climate {

std::vector<double> loess_smooth(const std::vector<double>& values, int span, int degree,
                                 const std::vector<double>& robustness) {
    const std::size_t n = values.size();
    std::vector<double> result(n, 0.0);
    if (n == 0) return result;
    if (degree < 0) throw std::runtime_error("loess_smooth: degree must be non-negative");

    const bool have_robustness = robustness.size() == n;
    int window = std::max(span, degree + 1);
    window = std::max(window, 1);
    if (static_cast<std::size_t>(window) > n) window = static_cast<int>(n);

    for (std::size_t i = 0; i < n; ++i) {
        int left = static_cast<int>(i) - window / 2;
        int right = left + window - 1;
        if (left < 0) {
            left = 0;
            right = window - 1;
        }
        if (right >= static_cast<int>(n)) {
            right = static_cast<int>(n) - 1;
            left = right - window + 1;
        }

        double max_distance = 0.0;
        for (int j = left; j <= right; ++j)
            max_distance = std::max(max_distance, std::fabs(static_cast<double>(j) - static_cast<double>(i)));
        if (max_distance <= 0.0) {
            result[i] = values[i];
            continue;
        }

        const int terms = degree + 1;
        std::vector<std::vector<double>> normal(
            static_cast<std::size_t>(terms), std::vector<double>(static_cast<std::size_t>(terms), 0.0));
        std::vector<double> rhs(static_cast<std::size_t>(terms), 0.0);

        double total_weight = 0.0;
        for (int j = left; j <= right; ++j) {
            if (!std::isfinite(values[static_cast<std::size_t>(j)])) continue;
            const double distance =
                std::fabs(static_cast<double>(j) - static_cast<double>(i)) / max_distance;
            if (distance >= 1.0) continue;
            const double u = 1.0 - distance * distance * distance;
            double weight = u * u * u;
            if (have_robustness) weight *= robustness[static_cast<std::size_t>(j)];
            if (weight <= 0.0) continue;
            total_weight += weight;

            const double x = static_cast<double>(j) - static_cast<double>(i);
            std::vector<double> powers(static_cast<std::size_t>(2 * terms - 1), 1.0);
            for (std::size_t p = 1; p < powers.size(); ++p) powers[p] = powers[p - 1] * x;

            for (int row = 0; row < terms; ++row) {
                for (int col = 0; col < terms; ++col) {
                    normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] +=
                        weight * powers[static_cast<std::size_t>(row + col)];
                }
                rhs[static_cast<std::size_t>(row)] +=
                    weight * powers[static_cast<std::size_t>(row)] *
                    values[static_cast<std::size_t>(j)];
            }
        }

        if (total_weight <= 0.0) {
            result[i] = values[i];
            continue;
        }

        std::vector<double> coefficients;
        if (solve_linear_system(normal, rhs, coefficients))
            result[i] = coefficients[0];
        else
            result[i] = values[i];
    }
    return result;
}

std::vector<double> moving_average(const std::vector<double>& values, int window) {
    const std::size_t n = values.size();
    std::vector<double> result(n, 0.0);
    if (n == 0 || window <= 0) return result;
    if (window == 1) return values;

    if (window % 2 == 1) {
        const int half = window / 2;
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            int count = 0;
            for (int k = -half; k <= half; ++k) {
                long index = static_cast<long>(i) + k;
                if (index < 0) index = 0;
                if (index >= static_cast<long>(n)) index = static_cast<long>(n) - 1;
                sum += values[static_cast<std::size_t>(index)];
                ++count;
            }
            result[i] = sum / static_cast<double>(count);
        }
    } else {
        const int half = window / 2;
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            double weight_sum = 0.0;
            for (int k = -half; k <= half; ++k) {
                const double weight = (k == -half || k == half) ? 0.5 : 1.0;
                long index = static_cast<long>(i) + k;
                if (index < 0) index = 0;
                if (index >= static_cast<long>(n)) index = static_cast<long>(n) - 1;
                sum += weight * values[static_cast<std::size_t>(index)];
                weight_sum += weight;
            }
            result[i] = sum / weight_sum;
        }
    }
    return result;
}

double median_absolute_deviation(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    std::vector<double> magnitudes;
    magnitudes.reserve(values.size());
    const double centre = median_value(values);
    for (double v : values) magnitudes.push_back(std::fabs(v - centre));
    return median_value(magnitudes);
}

}  // namespace climate
