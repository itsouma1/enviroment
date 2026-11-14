// Generates the synthetic demonstration datasets shipped in data/.
// The series are physically plausible but synthetic: they exist so that every
// command in the tool can be exercised end to end without downloading NASA
// products. Replace them with real CDR files for scientific use.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../src/common.hpp"
#include "../src/retrieval.hpp"

namespace {

std::mt19937_64 rng(20251114u);
std::normal_distribution<double> standard_normal(0.0, 1.0);

double noise(double sigma) { return sigma * standard_normal(rng); }

std::string number(double value, int precision = 6) {
    // Fixed-point formatting: significant-digit formatting would collapse all
    // twelve months of a year onto the same timestamp.
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

void write_file(const std::filesystem::path& path, const std::string& content) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot write " + path.string());
    output << content;
    std::cout << "wrote " << path << " (" << content.size() << " bytes)\n";
}

double clamp01(double value) { return value < 0.0 ? 0.0 : (value > 1.0 ? 1.0 : value); }

}  // namespace

int main(int argc, char** argv) {
    const std::filesystem::path directory = argc > 1 ? argv[1] : "data";
    std::filesystem::create_directories(directory);

    // Sea surface temperature, 1982-2024, monthly, with a warming trend and a
    // few missing values.
    {
        std::ostringstream out;
        out << "time,sst\n";
        for (int year = 1982; year <= 2024; ++year) {
            for (int month = 1; month <= 12; ++month) {
                const double t = year + (month - 0.5) / 12.0;
                const double sst = 18.1 + 0.014 * (t - 1982.0) +
                                   0.32 * std::sin(2.0 * climate::PI * (t - 0.12)) + noise(0.05);
                out << number(t, 4) << ",";
                if (month == 7 && year % 13 == 0)
                    out << "NA";
                else
                    out << number(sst, 4);
                out << "\n";
            }
        }
        write_file(directory / "sst_monthly.csv", out.str());
    }

    // Global surface air temperature anomaly, 1950-2024, monthly.
    {
        std::ostringstream out;
        out << "time,value\n";
        for (int year = 1950; year <= 2024; ++year) {
            for (int month = 1; month <= 12; ++month) {
                const double t = year + (month - 0.5) / 12.0;
                const double anomaly = -0.05 + 0.021 * (t - 1980.0) +
                                       0.10 * std::sin(2.0 * climate::PI * t) + noise(0.07);
                out << number(t, 4) << "," << number(anomaly, 4) << "\n";
            }
        }
        write_file(directory / "air_temperature_monthly.csv", out.str());
    }

    // Atmospheric CO2, 1980-2024, monthly: quadratic trend plus four harmonics,
    // a deliberately contaminated point and several gaps.
    {
        std::ostringstream out;
        out << "time,co2\n";
        for (int year = 1980; year <= 2024; ++year) {
            for (int month = 1; month <= 12; ++month) {
                const double t = year + (month - 0.5) / 12.0;
                const double dt = t - 1980.0;
                double co2 = 338.0 + 1.55 * dt + 0.012 * dt * dt +
                             3.0 * std::sin(2.0 * climate::PI * t) +
                             1.1 * std::sin(4.0 * climate::PI * t + 0.7) + noise(0.25);
                const bool missing = (year == 1985 && month == 3) || (year == 1999 && month == 9) ||
                                     (year == 2011 && month == 1);
                const bool outlier = (year == 1991 && month == 6);
                out << number(t, 4) << ",";
                if (missing)
                    out << "NA";
                else if (outlier)
                    out << number(co2 - 9.0, 4);
                else
                    out << number(co2, 4);
                out << "\n";
            }
        }
        write_file(directory / "co2_monthly.csv", out.str());
    }

    // Amazon NDVI, 2000-2024, monthly.
    {
        std::ostringstream out;
        out << "time,value\n";
        for (int year = 2000; year <= 2024; ++year) {
            for (int month = 1; month <= 12; ++month) {
                const double t = year + (month - 0.5) / 12.0;
                const double ndvi = 0.74 + 0.0035 * (t - 2000.0) +
                                    0.05 * std::sin(2.0 * climate::PI * (t - 0.18)) + noise(0.008);
                out << number(t, 4) << "," << number(ndvi, 5) << "\n";
            }
        }
        write_file(directory / "ndvi_monthly.csv", out.str());
    }

    // Tropical precipitation, 1979-2024, monthly, mm/day.
    {
        std::ostringstream out;
        out << "time,value\n";
        for (int year = 1979; year <= 2024; ++year) {
            for (int month = 1; month <= 12; ++month) {
                const double t = year + (month - 0.5) / 12.0;
                const double precipitation = 3.2 + 1.1 * std::sin(2.0 * climate::PI * (t - 0.2)) +
                                             0.0015 * (t - 1980.0) + noise(0.25);
                out << number(t, 4) << "," << number(precipitation, 5) << "\n";
            }
        }
        write_file(directory / "precipitation_monthly.csv", out.str());
    }

    // Arctic September sea ice concentration grid, 1 degree cells. The ice edge
    // sits near 79 N, so the 15% extent is a realistic fraction of the grid.
    {
        std::ostringstream out;
        out << "lat_min,lat_max,lon_min,lon_max,concentration\n";
        for (int lat = 68; lat <= 88; ++lat) {
            for (int lon = -180; lon < 180; lon += 5) {
                double concentration =
                    0.02 + 0.14 * (static_cast<double>(lat) + 1.0 - 78.5) +
                    0.04 * std::sin(climate::deg_to_rad(3.0 * lon));
                concentration = clamp01(concentration + noise(0.03));
                out << lat << "," << lat + 1 << "," << lon << "," << lon + 5 << ","
                    << number(concentration, 4) << "\n";
            }
        }
        write_file(directory / "sea_ice_grid_2024_09.csv", out.str());
    }

    // Arctic sea ice extent, 1979-2024, monthly, million km2: seasonal cycle plus
    // a declining trend of roughly 0.45 million km2 per decade.
    {
        std::ostringstream out;
        out << "time,value\n";
        for (int year = 1979; year <= 2024; ++year) {
            for (int month = 1; month <= 12; ++month) {
                const double t = year + (month - 0.5) / 12.0;
                const double extent = 10.0 - 0.045 * (t - 1979.0) +
                                      4.5 * std::cos(2.0 * climate::PI * (t - 0.12)) +
                                      noise(0.15);
                out << number(t, 4) << "," << number(extent, 4) << "\n";
            }
        }
        write_file(directory / "sea_ice_extent_monthly.csv", out.str());
    }

    // Satellite-buoy matchups for NLSST calibration. Truth is generated from a
    // known coefficient vector, so calibration accuracy can be verified.
    {
        const double truth[7] = {1.2, 0.995, 0.00035, 0.02, 0.01, -0.05, 0.04};
        std::uniform_real_distribution<double> zenith(0.0, 55.0);
        std::uniform_real_distribution<double> azimuth(-180.0, 180.0);
        std::uniform_real_distribution<double> guess(298.5, 301.5);
        std::uniform_int_distribution<int> mirror(0, 1);

        std::ostringstream out;
        out << "bt11,bt12,first_guess,sat_zenith_deg,mirror,theta_star_deg,buoy_sst\n";
        for (int i = 0; i < 2000; ++i) {
            const double tsfc = guess(rng);
            const double bt11 = tsfc - 0.6 + noise(0.4);
            const double bt12 = bt11 - 0.9 + noise(0.3);
            const double zenith_deg = zenith(rng);
            const double azimuth_deg = azimuth(rng);
            const int mirror_side = mirror(rng);

            const std::vector<double> features =
                climate::nlsst_features(bt11, bt12, tsfc, climate::deg_to_rad(zenith_deg),
                                        static_cast<double>(mirror_side),
                                        climate::deg_to_rad(azimuth_deg));
            double buoy = 0.0;
            for (std::size_t j = 0; j < 7; ++j) buoy += truth[j] * features[j];
            buoy += noise(0.05);

            out << number(bt11, 4) << "," << number(bt12, 4) << "," << number(tsfc, 4) << ","
                << number(zenith_deg, 3) << "," << mirror_side << "," << number(azimuth_deg, 3)
                << "," << number(buoy, 4) << "\n";
        }
        write_file(directory / "sst_matchups.csv", out.str());
    }

    // Passive-microwave tie points (representative example values).
    {
        std::ostringstream out;
        out << "constituent,channel_19v,channel_19h,channel_37v\n";
        out << "open_water,170.0,105.0,195.0\n";
        out << "first_year,245.0,235.0,245.0\n";
        out << "multiyear,215.0,200.0,200.0\n";
        write_file(directory / "nasa_team_tiepoints.csv", out.str());
    }

    // Independent precipitation estimates for the GPCP merge.
    {
        std::ostringstream out;
        out << "source,value,error_sd\n";
        out << "gpi_ir,2.80,0.90\n";
        out << "ssmis,3.40,0.60\n";
        out << "gauge,3.10,0.35\n";
        write_file(directory / "precip_estimates.csv", out.str());
    }

    // Meridional velocity section at 26.5 N for the AMOC diagnostic.
    {
        std::ostringstream out;
        out << "latitude,longitude,depth_m,velocity_m_s\n";
        for (int depth = 0; depth <= 1500; depth += 100) {
            for (int lon = -80; lon <= 0; lon += 2) {
                const double western =
                    0.75 * std::exp(-std::pow((lon + 76.0) / 1.6, 2.0)) * std::exp(-depth / 700.0);
                const double interior =
                    -0.02 * std::exp(-std::pow((lon + 40.0) / 30.0, 2.0)) * std::exp(-depth / 1200.0);
                out << "26.5," << lon << "," << depth << "," << number(western + interior, 6)
                    << "\n";
            }
        }
        write_file(directory / "amoc_velocity_26n.csv", out.str());
    }

    // Reflectance pairs for the NDVI retrieval.
    {
        std::ostringstream out;
        out << "red,nir\n";
        for (int i = 0; i < 200; ++i) {
            const double red = 0.06 + 0.10 * std::fabs(noise(0.5));
            const double nir = 0.35 + 0.25 * std::fabs(noise(0.5));
            out << number(red, 4) << "," << number(nir, 4) << "\n";
        }
        write_file(directory / "reflectance_pairs.csv", out.str());
    }

    std::cout << "demo data written to " << directory << "\n";
    return 0;
}
