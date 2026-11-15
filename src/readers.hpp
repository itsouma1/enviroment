#ifndef CLIMATE_READERS_HPP
#define CLIMATE_READERS_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "flow.hpp"
#include "ml.hpp"

namespace climate {

struct VelocityGrid {
    std::vector<double> longitude;
    std::vector<double> latitude;
    std::vector<double> u;  // latitude rows, longitude columns
    std::vector<double> v;
    double fill_value = -9999.0;
};

double normalize_longitude(double longitude);

// Four bilinear weights ordered (x0,y0), (x1,y0), (x0,y1), (x1,y1).
std::vector<double> bilinear_weights(double x, double y, double x0, double x1, double y0,
                                     double y1);

double bilinear_sample(const std::vector<double>& field, std::size_t rows,
                       std::size_t columns, const std::vector<double>& longitudes,
                       const std::vector<double>& latitudes, double longitude,
                       double latitude);

double bilinear_time_weight(double time, double first_time, double second_time);

// Builds a VelocityField from co-located u/v grids; coordinates are ascending.
VelocityField bilinear_velocity(const VelocityGrid& grid);

// Reads the first time slice of u/v (or EWCT/NSCT) from a classic NetCDF file.
// `u_name`/`v_name` may be empty to use the default aliases.
VelocityGrid read_velocity_netcdf(const std::string& path, const std::string& u_name,
                                  const std::string& v_name, std::size_t time_index);

// Reads 10 m wind (u10/v10), converting longitude to [-180, 180).
VelocityGrid read_era5_wind_netcdf(const std::string& path, std::size_t time_index);

// Mass-weighted spatial average over a grid; used for HFR/ERA5 diagnostics.
double cosine_latitude_mean(const std::vector<double>& field, std::size_t rows,
                            std::size_t columns, const std::vector<double>& latitudes);

}  // namespace climate

#endif
