#include "flow.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::pair<double, double> advance(double x, double y, double t, double dt,
                                  const VelocityField& velocity) {
    const auto k1 = velocity(x, y, t);
    const auto k2 = velocity(x + 0.5 * dt * k1.first, y + 0.5 * dt * k1.second, t + 0.5 * dt);
    const auto k3 = velocity(x + 0.5 * dt * k2.first, y + 0.5 * dt * k2.second, t + 0.5 * dt);
    const auto k4 = velocity(x + dt * k3.first, y + dt * k3.second, t + dt);
    return {x + dt / 6.0 * (k1.first + 2.0 * k2.first + 2.0 * k3.first + k4.first),
            y + dt / 6.0 * (k1.second + 2.0 * k2.second + 2.0 * k3.second + k4.second)};
}

std::pair<double, double> flow_map(double x0, double y0, double t0, double integration_time,
                                   const VelocityField& velocity, int steps) {
    const double dt = integration_time / static_cast<double>(steps);
    double x = x0, y = y0, t = t0;
    for (int i = 0; i < steps; ++i) {
        const auto next = advance(x, y, t, dt, velocity);
        x = next.first;
        y = next.second;
        t += dt;
    }
    return {x, y};
}

}  // namespace

double ftle(double x0, double y0, double t0, double integration_time,
            const VelocityField& velocity, double jacobian_step, int steps) {
    require(std::fabs(integration_time) > 0.0 && jacobian_step > 0.0 && steps > 0,
            "ftle: inputs");
    const auto centre = flow_map(x0, y0, t0, integration_time, velocity, steps);
    const auto east = flow_map(x0 + jacobian_step, y0, t0, integration_time, velocity, steps);
    const auto north = flow_map(x0, y0 + jacobian_step, t0, integration_time, velocity, steps);
    const double a = (east.first - centre.first) / jacobian_step;
    const double c = (east.second - centre.second) / jacobian_step;
    const double b = (north.first - centre.first) / jacobian_step;
    const double d = (north.second - centre.second) / jacobian_step;
    const double xx = a * a + c * c;
    const double xy = a * b + c * d;
    const double yy = b * b + d * d;
    const double mean = 0.5 * (xx + yy);
    const double radius = std::sqrt(0.25 * (xx - yy) * (xx - yy) + xy * xy);
    const double largest = std::max(1e-15, mean + radius);
    return 0.5 * std::log(largest) / std::fabs(integration_time);
}

Matrix ftle_field(double t0, double integration_time, double x_min, double x_max,
                  double y_min, double y_max, std::size_t nx, std::size_t ny,
                  const VelocityField& velocity, double jacobian_step, int steps) {
    require(nx > 1 && ny > 1 && x_max > x_min && y_max > y_min,
            "ftle_field: grid");
    Matrix field(ny, std::vector<double>(nx, 0.0));
    for (std::size_t j = 0; j < ny; ++j) {
        const double y = y_min + (y_max - y_min) * static_cast<double>(j) /
                                     static_cast<double>(ny - 1);
        for (std::size_t i = 0; i < nx; ++i) {
            const double x = x_min + (x_max - x_min) * static_cast<double>(i) /
                                         static_cast<double>(nx - 1);
            field[j][i] = ftle(x, y, t0, integration_time, velocity, jacobian_step, steps);
        }
    }
    return field;
}

double lagrangian_descriptor(double x0, double y0, double t0, double tau,
                             const VelocityField& velocity, int steps) {
    require(tau >= 0.0 && steps > 0, "lagrangian_descriptor: inputs");
    double length = 0.0;
    for (int direction = -1; direction <= 1; direction += 2) {
        double x = x0, y = y0, t = t0;
        const double dt = direction * tau / static_cast<double>(steps);
        for (int i = 0; i < steps; ++i) {
            const auto next = advance(x, y, t, dt, velocity);
            length += std::hypot(next.first - x, next.second - y);
            x = next.first;
            y = next.second;
            t += dt;
        }
    }
    return length;
}

std::pair<double, double> jet_axis(double t0, double tau, double x_min, double x_max,
                                   double y_min, double y_max, std::size_t nx,
                                   std::size_t ny, const VelocityField& velocity, int steps) {
    require(nx > 1 && ny > 1, "jet_axis: grid");
    double best = -std::numeric_limits<double>::infinity();
    std::pair<double, double> axis{x_min, y_min};
    for (std::size_t j = 0; j < ny; ++j) {
        const double y = y_min + (y_max - y_min) * static_cast<double>(j) /
                                     static_cast<double>(ny - 1);
        for (std::size_t i = 0; i < nx; ++i) {
            const double x = x_min + (x_max - x_min) * static_cast<double>(i) /
                                         static_cast<double>(nx - 1);
            const double value = lagrangian_descriptor(x, y, t0, tau, velocity, steps);
            if (value > best) {
                best = value;
                axis = {x, y};
            }
        }
    }
    return axis;
}

double occurrence_percent(std::size_t strong_days, std::size_t sampled_days) {
    require(sampled_days > 0, "occurrence_percent: no samples");
    return 100.0 * static_cast<double>(strong_days) / static_cast<double>(sampled_days);
}

double year_to_year_variability(const std::vector<double>& values) {
    require(values.size() > 1, "year_to_year_variability: size");
    double mean = 0.0;
    for (double value : values) mean += value;
    mean /= static_cast<double>(values.size());
    double sum = 0.0;
    for (double value : values) sum += (value - mean) * (value - mean);
    return std::sqrt(sum / static_cast<double>(values.size() - 1));
}

double regime_shift_difference(double post_mean, double pre_mean) {
    return post_mean - pre_mean;
}

std::vector<double> cyclic_moving_average(const std::vector<double>& values, std::size_t window) {
    require(!values.empty() && window > 0 && window % 2 == 1, "cyclic_moving_average: inputs");
    const std::size_t n = values.size();
    const std::size_t half = window / 2;
    std::vector<double> smoothed(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (std::size_t offset = 0; offset < window; ++offset) {
            const long index = static_cast<long>(i) + static_cast<long>(offset) -
                               static_cast<long>(half);
            const long wrapped = ((index % static_cast<long>(n)) + static_cast<long>(n)) %
                                 static_cast<long>(n);
            sum += values[static_cast<std::size_t>(wrapped)];
        }
        smoothed[i] = sum / static_cast<double>(window);
    }
    return smoothed;
}

bool can_stitch(double distance, double time_gap, double radius, double max_gap) {
    return distance <= radius && time_gap <= max_gap;
}

std::vector<std::vector<std::pair<std::size_t, std::size_t>>> stitch_tracks(
    const std::vector<std::vector<std::pair<double, double>>>& frames, double radius,
    double max_gap) {
    require(radius > 0.0 && max_gap >= 0.0, "stitch_tracks: inputs");
    struct Active {
        std::size_t track = 0;
        std::size_t frame = 0;
        std::pair<double, double> position;
    };
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> tracks;
    std::vector<Active> active;
    for (std::size_t f = 0; f < frames.size(); ++f) {
        std::vector<Active> next_active;
        std::vector<int> used(frames[f].size(), 0);
        for (const Active& track : active) {
            if (static_cast<double>(f - track.frame) > max_gap) continue;
            double best = radius;
            int best_index = -1;
            for (std::size_t c = 0; c < frames[f].size(); ++c) {
                if (used[c]) continue;
                const double distance = std::hypot(frames[f][c].first - track.position.first,
                                                   frames[f][c].second - track.position.second);
                if (distance <= best) {
                    best = distance;
                    best_index = static_cast<int>(c);
                }
            }
            if (best_index >= 0) {
                used[static_cast<std::size_t>(best_index)] = 1;
                tracks[track.track].push_back({f, static_cast<std::size_t>(best_index)});
                next_active.push_back({track.track, f, frames[f][static_cast<std::size_t>(best_index)]});
            } else if (static_cast<double>(f - track.frame) <= max_gap) {
                next_active.push_back({track.track, track.frame, track.position});
            }
        }
        for (std::size_t c = 0; c < frames[f].size(); ++c) {
            if (used[c]) continue;
            tracks.push_back({{f, c}});
            next_active.push_back({tracks.size() - 1, f, frames[f][c]});
        }
        active = next_active;
    }
    return tracks;
}

double jet_trend_correction(double dynamical_trend, double latitude_sensitivity,
                            double isentrope_trend) {
    return dynamical_trend + latitude_sensitivity * isentrope_trend;
}

Matrix sobel_gradient(const Matrix& field, int axis) {
    require(!field.empty() && !field[0].empty() && axis >= 0 && axis <= 1,
            "sobel_gradient: inputs");
    const std::size_t rows = field.size(), columns = field[0].size();
    Matrix gradient(rows, std::vector<double>(columns, 0.0));
    if (rows < 3 || columns < 3) return gradient;
    const double kernel_x[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    const double kernel_y[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
    for (std::size_t r = 1; r + 1 < rows; ++r)
        for (std::size_t c = 1; c + 1 < columns; ++c) {
            double sum = 0.0;
            for (int i = -1; i <= 1; ++i)
                for (int j = -1; j <= 1; ++j) {
                    const double weight =
                        axis == 0 ? kernel_x[i + 1][j + 1] : kernel_y[i + 1][j + 1];
                    sum += weight * field[static_cast<std::size_t>(static_cast<long>(r) + i)]
                                          [static_cast<std::size_t>(static_cast<long>(c) + j)];
                }
            gradient[r][c] = 0.125 * sum;
        }
    return gradient;
}

double total_sobel_signal(double spatial_x, double spatial_y, double temporal,
                          double calibration) {
    const double spatial = calibration * (spatial_x + spatial_y);
    return std::sqrt(spatial * spatial + temporal * temporal);
}

double calibration_factor(double spatial_quantile, double temporal_quantile) {
    require(spatial_quantile > 0.0 && temporal_quantile > 0.0,
            "calibration_factor: quantiles");
    return spatial_quantile / temporal_quantile;
}

std::vector<double> non_maximum_suppression_1d(const std::vector<double>& values) {
    std::vector<double> suppressed(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const bool left = i == 0 || values[i] > values[i - 1];
        const bool right = i + 1 == values.size() || values[i] > values[i + 1];
        if (left && right) suppressed[i] = values[i];
    }
    return suppressed;
}

std::vector<int> hysteresis_threshold_1d(const std::vector<double>& values, double high,
                                         double low) {
    require(high >= low, "hysteresis_threshold_1d: thresholds");
    std::vector<int> result(values.size(), 0);
    for (std::size_t i = 0; i < values.size(); ++i)
        if (values[i] >= high) result[i] = 1;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (result[i] == 1) continue;
        if (values[i] < low) continue;
        const bool touches_strong =
            (i > 0 && result[i - 1] == 1) || (i + 1 < values.size() && result[i + 1] == 1);
        if (touches_strong) result[i] = 2;
    }
    return result;
}

double cohens_d(double before_mean, double after_mean, double pooled_sd) {
    require(pooled_sd > 0.0, "cohens_d: pooled sd");
    return (after_mean - before_mean) / pooled_sd;
}

std::size_t connected_component_count(const std::vector<int>& mask, std::size_t nx,
                                      std::size_t ny, std::size_t nz) {
    require(mask.size() == nx * ny * nz, "connected_component_count: size");
    std::vector<int> visited(mask.size(), 0);
    std::size_t components = 0;
    std::vector<std::size_t> stack;
    for (std::size_t start = 0; start < mask.size(); ++start) {
        if (mask[start] == 0 || visited[start]) continue;
        ++components;
        stack.push_back(start);
        visited[start] = 1;
        while (!stack.empty()) {
            const std::size_t index = stack.back();
            stack.pop_back();
            const std::size_t x = index % nx;
            const std::size_t y = (index / nx) % ny;
            const std::size_t z = index / (nx * ny);
            const int dx[6] = {-1, 1, 0, 0, 0, 0};
            const int dy[6] = {0, 0, -1, 1, 0, 0};
            const int dz[6] = {0, 0, 0, 0, -1, 1};
            for (int k = 0; k < 6; ++k) {
                const long nxx = static_cast<long>(x) + dx[k];
                const long nyy = static_cast<long>(y) + dy[k];
                const long nzz = static_cast<long>(z) + dz[k];
                if (nxx < 0 || nyy < 0 || nzz < 0 || nxx >= static_cast<long>(nx) ||
                    nyy >= static_cast<long>(ny) || nzz >= static_cast<long>(nz))
                    continue;
                const std::size_t neighbour = static_cast<std::size_t>(nzz) * nx * ny +
                                              static_cast<std::size_t>(nyy) * nx +
                                              static_cast<std::size_t>(nxx);
                if (mask[neighbour] != 0 && !visited[neighbour]) {
                    visited[neighbour] = 1;
                    stack.push_back(neighbour);
                }
            }
        }
    }
    return components;
}

}  // namespace climate
