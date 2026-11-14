#pragma once

#include <cmath>
#include <string>

namespace climate {

// Physical and numerical constants used across the retrieval and analysis code.
constexpr double PI = 3.14159265358979323846264338327950288;
constexpr double EARTH_RADIUS_KM = 6371.0088;
constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;   // W m^-2 K^-4
constexpr double CARBON_PPM_TO_PGC = 2.124;           // PgC per ppm of atmospheric CO2
constexpr double GPI_MM_PER_HOUR = 3.0;               // GPI rain rate for full cold-cloud coverage
constexpr double GPI_COLD_CLOUD_THRESHOLD_K = 235.0;  // IR threshold used by GPI
constexpr double SIC_EXTENT_THRESHOLD = 0.15;         // 15% sea ice concentration convention
constexpr double DAYS_PER_YEAR = 365.2425;

inline double deg_to_rad(double degrees) { return degrees * PI / 180.0; }
inline double rad_to_deg(double radians) { return radians * 180.0 / PI; }

inline double clamp_value(double value, double low, double high) {
    if (value < low) return low;
    if (value > high) return high;
    return value;
}

inline bool is_missing(double value) { return !std::isfinite(value); }

}  // namespace climate
