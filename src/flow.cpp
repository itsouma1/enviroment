#include "flow.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
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

std::vector<std::pair<std::size_t, std::size_t>> constrained_ridge_path(
    const Matrix& field, std::size_t max_step, double max_difference) {
    require(!field.empty() && !field[0].empty() && max_step > 0,
            "constrained_ridge_path: inputs");
    const std::size_t rows = field.size(), columns = field[0].size();
    std::size_t best_row = 0, best_column = 0;
    double best_value = field[0][0];
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < columns; ++c)
            if (field[r][c] > best_value) {
                best_value = field[r][c];
                best_row = r;
                best_column = c;
            }
    std::vector<std::pair<std::size_t, std::size_t>> path{{best_row, best_column}};
    std::vector<int> visited(rows * columns, 0);
    visited[best_row * columns + best_column] = 1;
    while (true) {
        const std::size_t last_row = path.back().first, last_column = path.back().second;
        int next_row = -1, next_column = -1;
        double next_value = -1e300;
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t c = 0; c < columns; ++c) {
                if (visited[r * columns + c]) continue;
                const std::size_t dr = r > last_row ? r - last_row : last_row - r;
                const std::size_t dc = c > last_column ? c - last_column : last_column - c;
                if (std::max(dr, dc) > max_step) continue;
                if (std::fabs(field[r][c] - field[last_row][last_column]) > max_difference)
                    continue;
                bool local_maximum = true;
                for (int i = -1; i <= 1 && local_maximum; ++i)
                    for (int j = -1; j <= 1; ++j) {
                        const long rr = static_cast<long>(r) + i;
                        const long cc = static_cast<long>(c) + j;
                        if (rr < 0 || cc < 0 || rr >= static_cast<long>(rows) ||
                            cc >= static_cast<long>(columns))
                            continue;
                        if (field[static_cast<std::size_t>(rr)][static_cast<std::size_t>(cc)] >
                            field[r][c])
                            local_maximum = false;
                    }
                if (!local_maximum) continue;
                if (field[r][c] > next_value) {
                    next_value = field[r][c];
                    next_row = static_cast<int>(r);
                    next_column = static_cast<int>(c);
                }
            }
        if (next_row < 0) break;
        path.push_back({static_cast<std::size_t>(next_row),
                        static_cast<std::size_t>(next_column)});
        visited[static_cast<std::size_t>(next_row) * columns +
                static_cast<std::size_t>(next_column)] = 1;
    }
    return path;
}

bool has_closed_contour(const Matrix& field, std::size_t row, std::size_t column,
                        double threshold, std::size_t max_radius) {
    require(!field.empty() && !field[0].empty() && row < field.size() &&
                column < field[0].size(),
            "has_closed_contour: index");
    const double target = field[row][column] + threshold;
    const std::size_t rows = field.size(), columns = field[0].size();
    for (std::size_t radius = 1; radius <= max_radius; ++radius) {
        if (radius > row || radius > column || row + radius >= rows ||
            column + radius >= columns)
            break;
        bool closed = true;
        for (long offset = -static_cast<long>(radius);
             offset <= static_cast<long>(radius) && closed; ++offset) {
            const long r1 = static_cast<long>(row) - static_cast<long>(radius);
            const long r2 = static_cast<long>(row) + static_cast<long>(radius);
            const long c1 = static_cast<long>(column) + offset;
            if (field[static_cast<std::size_t>(r1)][static_cast<std::size_t>(c1)] < target ||
                field[static_cast<std::size_t>(r2)][static_cast<std::size_t>(c1)] < target)
                closed = false;
            const long r3 = static_cast<long>(row) + offset;
            const long c2 = static_cast<long>(column) - static_cast<long>(radius);
            const long c3 = static_cast<long>(column) + static_cast<long>(radius);
            if (field[static_cast<std::size_t>(r3)][static_cast<std::size_t>(c2)] < target ||
                field[static_cast<std::size_t>(r3)][static_cast<std::size_t>(c3)] < target)
                closed = false;
        }
        if (closed) return true;
    }
    return false;
}

std::vector<std::pair<std::size_t, std::size_t>> longest_path_ridge(
    const Matrix& field, std::size_t max_step, double max_difference) {
    require(!field.empty() && !field[0].empty() && max_step > 0,
            "longest_path_ridge: inputs");
    const std::size_t rows = field.size(), columns = field[0].size();
    struct Vertex {
        std::size_t row;
        std::size_t column;
        double value;
    };
    std::vector<Vertex> vertices;
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < columns; ++c) {
            bool maximum = true;
            for (int dr = -1; dr <= 1 && maximum; ++dr)
                for (int dc = -1; dc <= 1; ++dc) {
                    if (dr == 0 && dc == 0) continue;
                    const long rr = static_cast<long>(r) + dr;
                    const long cc = static_cast<long>(c) + dc;
                    if (rr < 0 || cc < 0 || rr >= static_cast<long>(rows) ||
                        cc >= static_cast<long>(columns))
                        continue;
                    if (field[static_cast<std::size_t>(rr)][static_cast<std::size_t>(cc)] >
                        field[r][c])
                        maximum = false;
                }
            if (maximum) vertices.push_back({r, c, field[r][c]});
        }
    require(!vertices.empty(), "longest_path_ridge: no local maxima");
    std::sort(vertices.begin(), vertices.end(),
              [](const Vertex& a, const Vertex& b) { return a.value < b.value; });
    std::vector<double> score(vertices.size(), 0.0);
    std::vector<int> parent(vertices.size(), -1);
    for (std::size_t v = 0; v < vertices.size(); ++v) {
        score[v] = vertices[v].value;
        for (std::size_t u = 0; u < v; ++u) {
            const std::size_t dr = vertices[v].row > vertices[u].row
                                       ? vertices[v].row - vertices[u].row
                                       : vertices[u].row - vertices[v].row;
            const std::size_t dc = vertices[v].column > vertices[u].column
                                       ? vertices[v].column - vertices[u].column
                                       : vertices[u].column - vertices[v].column;
            if (std::max(dr, dc) > max_step) continue;
            if (std::fabs(vertices[v].value - vertices[u].value) > max_difference) continue;
            if (score[u] + vertices[v].value > score[v]) {
                score[v] = score[u] + vertices[v].value;
                parent[v] = static_cast<int>(u);
            }
        }
    }
    std::size_t best = 0;
    for (std::size_t v = 1; v < vertices.size(); ++v)
        if (score[v] > score[best]) best = v;
    std::vector<std::pair<std::size_t, std::size_t>> path;
    for (int v = static_cast<int>(best); v >= 0; v = parent[static_cast<std::size_t>(v)])
        path.push_back({vertices[static_cast<std::size_t>(v)].row,
                        vertices[static_cast<std::size_t>(v)].column});
    std::reverse(path.begin(), path.end());
    return path;
}

namespace {

std::pair<double, double> edge_crossing(const Matrix& field, std::size_t row,
                                        std::size_t column, int edge, double threshold) {
    const double bl = field[row][column], br = field[row][column + 1];
    const double tr = field[row + 1][column + 1], tl = field[row + 1][column];
    auto interpolate = [threshold](double a, double b) {
        const double denominator = b - a;
        if (std::fabs(denominator) < 1e-15) return 0.5;
        return (threshold - a) / denominator;
    };
    switch (edge) {
        case 0: return {static_cast<double>(column) + interpolate(bl, br),
                        static_cast<double>(row)};
        case 1: return {static_cast<double>(column + 1),
                        static_cast<double>(row) + interpolate(br, tr)};
        case 2: return {static_cast<double>(column) + interpolate(tl, tr),
                        static_cast<double>(row + 1)};
        default: return {static_cast<double>(column),
                         static_cast<double>(row) + interpolate(bl, tl)};
    }
}

std::vector<std::pair<int, int>> case_edges(int code, bool centre_high) {
    switch (code) {
        case 0:
        case 15: return {};
        case 1: return {{3, 0}};
        case 2: return {{0, 1}};
        case 3: return {{3, 1}};
        case 4: return {{1, 2}};
        case 5: return centre_high ? std::vector<std::pair<int, int>>{{3, 2}, {0, 1}}
                                   : std::vector<std::pair<int, int>>{{3, 0}, {1, 2}};
        case 6: return {{0, 2}};
        case 7: return {{3, 2}};
        case 8: return {{2, 3}};
        case 9: return {{2, 0}};
        case 10: return centre_high ? std::vector<std::pair<int, int>>{{3, 0}, {2, 1}}
                                    : std::vector<std::pair<int, int>>{{3, 2}, {0, 1}};
        case 11: return {{2, 1}};
        case 12: return {{3, 1}};
        case 13: return {{0, 1}};
        case 14: return {{0, 3}};
        default: return {};
    }
}

}  // namespace

std::vector<ContourSegment> marching_squares_segments(const Matrix& field,
                                                       double threshold) {
    require(field.size() >= 2 && field[0].size() >= 2, "marching_squares_segments: grid");
    std::vector<ContourSegment> segments;
    for (std::size_t row = 0; row + 1 < field.size(); ++row)
        for (std::size_t column = 0; column + 1 < field[0].size(); ++column) {
            int code = 0;
            if (field[row][column] >= threshold) code |= 1;
            if (field[row][column + 1] >= threshold) code |= 2;
            if (field[row + 1][column + 1] >= threshold) code |= 4;
            if (field[row + 1][column] >= threshold) code |= 8;
            const double centre = 0.25 * (field[row][column] + field[row][column + 1] +
                                          field[row + 1][column] + field[row + 1][column + 1]);
            for (const auto& edges : case_edges(code, centre >= threshold)) {
                const auto first = edge_crossing(field, row, column, edges.first, threshold);
                const auto second = edge_crossing(field, row, column, edges.second, threshold);
                segments.push_back({first.first, first.second, second.first, second.second});
            }
        }
    return segments;
}

std::vector<std::vector<std::pair<double, double>>> closed_contours(
    const Matrix& field, double threshold) {
    const std::vector<ContourSegment> segments = marching_squares_segments(field, threshold);
    std::map<std::pair<long long, long long>, std::vector<std::size_t>> adjacency;
    auto key = [](double x, double y) {
        return std::make_pair(static_cast<long long>(std::llround(x * 1e6)),
                              static_cast<long long>(std::llround(y * 1e6)));
    };
    for (std::size_t s = 0; s < segments.size(); ++s) {
        adjacency[key(segments[s].x1, segments[s].y1)].push_back(s);
        adjacency[key(segments[s].x2, segments[s].y2)].push_back(s);
    }
    std::vector<int> used(segments.size(), 0);
    std::vector<std::vector<std::pair<double, double>>> loops;
    for (std::size_t start = 0; start < segments.size(); ++start) {
        if (used[start]) continue;
        std::vector<std::pair<double, double>> loop;
        std::size_t current = start;
        double x = segments[start].x1, y = segments[start].y1;
        const double start_x = x, start_y = y;
        bool closed = false;
        while (true) {
            used[current] = 1;
            loop.push_back({x, y});
            const ContourSegment& segment = segments[current];
            const bool at_first = key(x, y) == key(segment.x1, segment.y1);
            x = at_first ? segment.x2 : segment.x1;
            y = at_first ? segment.y2 : segment.y1;
            if (key(x, y) == key(start_x, start_y)) {
                closed = true;
                break;
            }
            const auto found = adjacency.find(key(x, y));
            if (found == adjacency.end()) break;
            std::size_t next = segments.size();
            for (std::size_t candidate : found->second)
                if (!used[candidate]) {
                    next = candidate;
                    break;
                }
            if (next == segments.size()) break;
            current = next;
        }
        if (closed && loop.size() >= 3) loops.push_back(loop);
    }
    return loops;
}

std::vector<double> non_maximum_suppression_3d(const std::vector<double>& magnitude,
                                               const std::vector<double>& gradient_x,
                                               const std::vector<double>& gradient_y,
                                               const std::vector<double>& gradient_t,
                                               std::size_t nx, std::size_t ny,
                                               std::size_t nt, double step) {
    const std::size_t n = nx * ny * nt;
    require(magnitude.size() == n && gradient_x.size() == n && gradient_y.size() == n &&
                gradient_t.size() == n && step > 0.0,
            "non_maximum_suppression_3d: size");
    auto sample = [&](double x, double y, double z) {
        if (x < 0.0 || y < 0.0 || z < 0.0 || x >= static_cast<double>(nx - 1) ||
            y >= static_cast<double>(ny - 1) || z >= static_cast<double>(nt - 1))
            return std::numeric_limits<double>::quiet_NaN();
        const std::size_t i = static_cast<std::size_t>(x), j = static_cast<std::size_t>(y),
                          k = static_cast<std::size_t>(z);
        const double fx = x - static_cast<double>(i), fy = y - static_cast<double>(j),
                     fz = z - static_cast<double>(k);
        double sum = 0.0;
        for (int dz = 0; dz <= 1; ++dz)
            for (int dy = 0; dy <= 1; ++dy)
                for (int dx = 0; dx <= 1; ++dx) {
                    const double weight = (dx ? fx : 1.0 - fx) * (dy ? fy : 1.0 - fy) *
                                          (dz ? fz : 1.0 - fz);
                    sum += weight * magnitude[(i + dx) + nx * ((j + dy) + ny * (k + dz))];
                }
        return sum;
    };
    std::vector<double> result(n, 0.0);
    for (std::size_t k = 0; k < nt; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                const double gx = gradient_x[index], gy = gradient_y[index],
                             gt = gradient_t[index];
                const double norm = std::sqrt(gx * gx + gy * gy + gt * gt);
                if (norm <= 1e-15) continue;
                const double forward = sample(static_cast<double>(i) + step * gx / norm,
                                              static_cast<double>(j) + step * gy / norm,
                                              static_cast<double>(k) + step * gt / norm);
                const double backward = sample(static_cast<double>(i) - step * gx / norm,
                                               static_cast<double>(j) - step * gy / norm,
                                               static_cast<double>(k) - step * gt / norm);
                const double value = magnitude[index];
                bool keep = true;
                if (!std::isnan(forward) && value <= forward) keep = false;
                if (!std::isnan(backward) && value <= backward) keep = false;
                if (keep) result[index] = value;
            }
    return result;
}

std::vector<std::size_t> penalized_ridge_path(const Matrix& field, double penalty) {
    require(!field.empty() && !field[0].empty() && penalty >= 0.0,
            "penalized_ridge_path: inputs");
    const std::size_t rows = field.size(), columns = field[0].size();
    std::vector<std::vector<double>> cost(columns, std::vector<double>(rows, 0.0));
    std::vector<std::vector<int>> parent(columns, std::vector<int>(rows, -1));
    for (std::size_t row = 0; row < rows; ++row)
        cost[0][row] = -std::log(std::max(field[row][0], 1e-300));
    for (std::size_t column = 1; column < columns; ++column)
        for (std::size_t row = 0; row < rows; ++row) {
            double best = std::numeric_limits<double>::infinity();
            int best_parent = -1;
            for (std::size_t previous = 0; previous < rows; ++previous) {
                const double difference = static_cast<double>(row) -
                                          static_cast<double>(previous);
                const double candidate = cost[column - 1][previous] +
                                         penalty * difference * difference;
                if (candidate < best) {
                    best = candidate;
                    best_parent = static_cast<int>(previous);
                }
            }
            cost[column][row] = -std::log(std::max(field[row][column], 1e-300)) + best;
            parent[column][row] = best_parent;
        }
    std::size_t best_row = 0;
    for (std::size_t row = 1; row < rows; ++row)
        if (cost[columns - 1][row] < cost[columns - 1][best_row]) best_row = row;
    std::vector<std::size_t> path(columns, 0);
    int row = static_cast<int>(best_row);
    for (std::size_t column = columns; column-- > 0;) {
        path[column] = static_cast<std::size_t>(row);
        row = parent[column][static_cast<std::size_t>(row)];
        if (column == 0) break;
    }
    return path;
}

std::vector<std::size_t> cyclic_penalized_ridge(const Matrix& field, double penalty,
                                                int passes) {
    require(!field.empty() && !field[0].empty() && passes > 0,
            "cyclic_penalized_ridge: inputs");
    const std::size_t rows = field.size(), columns = field[0].size();
    std::vector<std::vector<std::size_t>> votes(columns, std::vector<std::size_t>(rows, 0));
    for (int pass = 0; pass < passes; ++pass) {
        const std::size_t offset = static_cast<std::size_t>(
            (static_cast<long>(pass) * static_cast<long>(columns)) / passes);
        Matrix shifted(rows, std::vector<double>(columns, 0.0));
        for (std::size_t row = 0; row < rows; ++row)
            for (std::size_t column = 0; column < columns; ++column)
                shifted[row][column] = field[row][(column + offset) % columns];
        const std::vector<std::size_t> path = penalized_ridge_path(shifted, penalty);
        for (std::size_t column = 0; column < columns; ++column)
            ++votes[(column + offset) % columns][path[column]];
    }
    std::vector<std::size_t> ridge(columns, 0);
    for (std::size_t column = 0; column < columns; ++column) {
        std::size_t best = 0;
        for (std::size_t row = 1; row < rows; ++row)
            if (votes[column][row] > votes[column][best]) best = row;
        ridge[column] = best;
    }
    return ridge;
}

double rossby_integration_time(double zonal_wavenumber, double beta) {
    require(beta > 0.0, "rossby_integration_time: beta");
    return std::fabs(zonal_wavenumber) / (2.0 * beta);
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
