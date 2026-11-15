#include "readers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "netcdf.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::size_t lower_bound_index(const std::vector<double>& axis, double value) {
    require(!axis.empty(), "lower_bound_index: empty axis");
    if (value <= axis.front()) return 0;
    if (value >= axis.back()) return axis.size() - 2;
    std::size_t low = 0, high = axis.size() - 1;
    while (high - low > 1) {
        const std::size_t middle = (low + high) / 2;
        if (axis[middle] <= value)
            low = middle;
        else
            high = middle;
    }
    return low;
}

std::vector<double> read_axis(const NetCDFFile& file, const std::string& first,
                              const std::string& second) {
    const NcVariable* variable = file.find_variable(first);
    if (variable == nullptr && !second.empty()) variable = file.find_variable(second);
    require(variable != nullptr, "read_axis: coordinate variable missing");
    return file.read_variable(variable->name);
}

const NcVariable* find_first(const NetCDFFile& file, const std::vector<std::string>& names) {
    for (const std::string& name : names)
        if (const NcVariable* variable = file.find_variable(name)) return variable;
    return nullptr;
}

}  // namespace

double normalize_longitude(double longitude) {
    double result = std::fmod(longitude + 180.0, 360.0);
    if (result < 0.0) result += 360.0;
    return result - 180.0;
}

std::vector<double> bilinear_weights(double x, double y, double x0, double x1, double y0,
                                     double y1) {
    require(x1 != x0 && y1 != y0, "bilinear_weights: degenerate cell");
    const double tx = (x - x0) / (x1 - x0);
    const double ty = (y - y0) / (y1 - y0);
    return {(1.0 - tx) * (1.0 - ty), tx * (1.0 - ty), (1.0 - tx) * ty, tx * ty};
}

double bilinear_sample(const std::vector<double>& field, std::size_t rows,
                       std::size_t columns, const std::vector<double>& longitudes,
                       const std::vector<double>& latitudes, double longitude,
                       double latitude) {
    require(field.size() == rows * columns && longitudes.size() == columns &&
                latitudes.size() == rows && columns >= 2 && rows >= 2,
            "bilinear_sample: size");
    const std::size_t i = lower_bound_index(longitudes, longitude);
    const std::size_t j = lower_bound_index(latitudes, latitude);
    const std::vector<double> weights = bilinear_weights(
        longitude, latitude, longitudes[i], longitudes[i + 1], latitudes[j], latitudes[j + 1]);
    return weights[0] * field[j * columns + i] + weights[1] * field[j * columns + i + 1] +
           weights[2] * field[(j + 1) * columns + i] +
           weights[3] * field[(j + 1) * columns + i + 1];
}

double bilinear_time_weight(double time, double first_time, double second_time) {
    require(second_time != first_time, "bilinear_time_weight: degenerate interval");
    return (time - first_time) / (second_time - first_time);
}

VelocityField bilinear_velocity(const VelocityGrid& grid) {
    return [grid](double longitude, double latitude, double) {
        return std::make_pair(
            bilinear_sample(grid.u, grid.latitude.size(), grid.longitude.size(),
                            grid.longitude, grid.latitude, longitude, latitude),
            bilinear_sample(grid.v, grid.latitude.size(), grid.longitude.size(),
                            grid.longitude, grid.latitude, longitude, latitude));
    };
}

VelocityGrid read_velocity_netcdf(const std::string& path, const std::string& u_name,
                                  const std::string& v_name, std::size_t time_index) {
    require(!is_hdf5_file(path),
            "read_velocity_netcdf: NetCDF-4/HDF5 must be converted to NetCDF-3");
    const NetCDFFile file = NetCDFFile::open(path);
    const NcVariable* u = u_name.empty()
                              ? find_first(file, {"u", "EWCT", "u10", "water_u"})
                              : file.find_variable(u_name);
    const NcVariable* v = v_name.empty()
                              ? find_first(file, {"v", "NSCT", "v10", "water_v"})
                              : file.find_variable(v_name);
    require(u != nullptr && v != nullptr, "read_velocity_netcdf: u/v variables missing");
    require(u->shape.size() >= 2 && v->shape.size() == u->shape.size(),
            "read_velocity_netcdf: unsupported rank");
    const std::size_t rows = u->shape[u->shape.size() - 2];
    const std::size_t columns = u->shape.back();
    const NetCDFFile& handle = file;
    VelocityGrid grid;
    grid.latitude = read_axis(handle, "lat", "latitude");
    grid.longitude = read_axis(handle, "lon", "longitude");
    if (grid.latitude.size() != rows) grid.latitude.resize(rows);
    if (grid.longitude.size() != columns) grid.longitude.resize(columns);
    const std::vector<double> raw_u = file.read_variable(u->name);
    const std::vector<double> raw_v = file.read_variable(v->name);
    std::size_t offset = 0;
    if (u->shape.size() == 3) {
        require(time_index < u->shape[0], "read_velocity_netcdf: time index");
        offset = time_index * rows * columns;
    }
    require(raw_u.size() >= offset + rows * columns && raw_v.size() >= offset + rows * columns,
            "read_velocity_netcdf: short variable");
    grid.u.assign(raw_u.begin() + static_cast<long>(offset),
                  raw_u.begin() + static_cast<long>(offset + rows * columns));
    grid.v.assign(raw_v.begin() + static_cast<long>(offset),
                  raw_v.begin() + static_cast<long>(offset + rows * columns));
    grid.fill_value = file.attribute_number(*u, "_FillValue", -9999.0);
    for (double& value : grid.u)
        if (value == grid.fill_value) value = std::numeric_limits<double>::quiet_NaN();
    for (double& value : grid.v)
        if (value == grid.fill_value) value = std::numeric_limits<double>::quiet_NaN();
    return grid;
}

VelocityGrid read_era5_wind_netcdf(const std::string& path, std::size_t time_index) {
    VelocityGrid grid = read_velocity_netcdf(path, "u10", "v10", time_index);
    for (double& longitude : grid.longitude) longitude = normalize_longitude(longitude);
    std::vector<std::size_t> order(grid.longitude.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return grid.longitude[a] < grid.longitude[b];
    });
    std::vector<double> sorted_longitude(grid.longitude.size());
    for (std::size_t i = 0; i < order.size(); ++i) sorted_longitude[i] = grid.longitude[order[i]];
    std::vector<double> sorted_u(grid.u.size()), sorted_v(grid.v.size());
    for (std::size_t row = 0; row < grid.latitude.size(); ++row)
        for (std::size_t column = 0; column < order.size(); ++column) {
            sorted_u[row * order.size() + column] = grid.u[row * order.size() + order[column]];
            sorted_v[row * order.size() + column] = grid.v[row * order.size() + order[column]];
        }
    grid.longitude = sorted_longitude;
    grid.u = sorted_u;
    grid.v = sorted_v;
    return grid;
}

double cosine_latitude_mean(const std::vector<double>& field, std::size_t rows,
                            std::size_t columns, const std::vector<double>& latitudes) {
    require(field.size() == rows * columns && latitudes.size() == rows && rows > 0,
            "cosine_latitude_mean: size");
    double weighted = 0.0, total = 0.0;
    for (std::size_t row = 0; row < rows; ++row) {
        const double weight = std::cos(latitudes[row] * 3.14159265358979323846 / 180.0);
        for (std::size_t column = 0; column < columns; ++column) {
            weighted += weight * field[row * columns + column];
            total += weight;
        }
    }
    require(total > 0.0, "cosine_latitude_mean: zero weight");
    return weighted / total;
}

}  // namespace climate
