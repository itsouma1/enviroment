#include "lagrangian.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

std::pair<double, double> hessian_eigenvalues(double xx, double xy, double yy) {
    const double mean = 0.5 * (xx + yy);
    const double radius = std::sqrt(0.25 * (xx - yy) * (xx - yy) + xy * xy);
    return {mean + radius, mean - radius};
}

double shape_index(double k1, double k2) {
    if (k1 < k2) std::swap(k1, k2);
    if (std::fabs(k2 - k1) < 1e-15) {
        if (k1 < 0.0) return 1.0;
        if (k1 > 0.0) return -1.0;
        return 0.0;
    }
    const double ratio = (k2 + k1) / (k2 - k1);
    return 2.0 / kPi * std::atan(ratio);
}

double shape_index_k1_derivative(double k1, double k2) {
    require(std::fabs(k2 - k1) > 1e-15, "shape_index_k1_derivative: isotropic curvature");
    const double u = (k2 + k1) / (k2 - k1);
    return 2.0 / kPi / (1.0 + u * u) * (-2.0 * k2) / ((k2 - k1) * (k2 - k1));
}

double shape_index_k2_derivative(double k1, double k2) {
    require(std::fabs(k2 - k1) > 1e-15, "shape_index_k2_derivative: isotropic curvature");
    const double u = (k2 + k1) / (k2 - k1);
    return 2.0 / kPi / (1.0 + u * u) * (2.0 * k1) / ((k2 - k1) * (k2 - k1));
}

Matrix gaussian_smooth(const Matrix& field, double sigma_cells) {
    require(!field.empty() && !field[0].empty(), "gaussian_smooth: empty field");
    require(sigma_cells > 0.0, "gaussian_smooth: sigma");
    const std::size_t rows = field.size(), columns = field[0].size();
    const int radius = std::max(1, static_cast<int>(std::ceil(3.0 * sigma_cells)));
    std::vector<double> kernel(static_cast<std::size_t>(2 * radius + 1));
    double normalisation = 0.0;
    for (int i = -radius; i <= radius; ++i) {
        kernel[static_cast<std::size_t>(i + radius)] =
            std::exp(-0.5 * (i / sigma_cells) * (i / sigma_cells));
        normalisation += kernel[static_cast<std::size_t>(i + radius)];
    }
    for (double& value : kernel) value /= normalisation;
    Matrix horizontal(rows, std::vector<double>(columns, 0.0));
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < columns; ++c) {
            double sum = 0.0;
            for (int k = -radius; k <= radius; ++k) {
                const int index = std::min<int>(static_cast<int>(columns) - 1,
                                                std::max(0, static_cast<int>(c) + k));
                sum += field[r][static_cast<std::size_t>(index)] *
                       kernel[static_cast<std::size_t>(k + radius)];
            }
            horizontal[r][c] = sum;
        }
    Matrix smoothed(rows, std::vector<double>(columns, 0.0));
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < columns; ++c) {
            double sum = 0.0;
            for (int k = -radius; k <= radius; ++k) {
                const int index = std::min<int>(static_cast<int>(rows) - 1,
                                                std::max(0, static_cast<int>(r) + k));
                sum += horizontal[static_cast<std::size_t>(index)][c] *
                       kernel[static_cast<std::size_t>(k + radius)];
            }
            smoothed[r][c] = sum;
        }
    return smoothed;
}

double linear_trend(const std::vector<double>& values) {
    require(values.size() > 1, "linear_trend: too few values");
    const double n = static_cast<double>(values.size());
    const double mean_t = 0.5 * (n - 1.0);
    double mean_y = 0.0;
    for (double value : values) mean_y += value;
    mean_y /= n;
    double covariance = 0.0, variance = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double dt = static_cast<double>(i) - mean_t;
        covariance += dt * (values[i] - mean_y);
        variance += dt * dt;
    }
    return covariance / variance;
}

double object_frequency(const std::vector<int>& present) {
    double count = 0.0;
    for (int flag : present) count += flag != 0 ? 1.0 : 0.0;
    return count;
}

double mean_object_intensity(const std::vector<double>& intensities,
                             const std::vector<int>& present) {
    require(intensities.size() == present.size(), "mean_object_intensity: size");
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < intensities.size(); ++i)
        if (present[i] != 0) {
            sum += intensities[i];
            ++count;
        }
    require(count > 0, "mean_object_intensity: no objects");
    return sum / static_cast<double>(count);
}

double nearest_link_distance(double x, double y, double next_x, double next_y, double radius) {
    require(radius > 0.0, "nearest_link_distance: radius");
    const double distance = std::hypot(next_x - x, next_y - y);
    return distance < radius ? distance : -1.0;
}

double object_population_trend(const std::vector<double>& counts) {
    return linear_trend(counts);
}

double along_track_trend(const std::vector<double>& intensities) {
    return linear_trend(intensities);
}

}  // namespace climate
