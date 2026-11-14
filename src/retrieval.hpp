#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace climate {

// ---------------------------------------------------------------------------
// Sea surface temperature: modified non-linear SST (Walton et al., 1998;
// Kilpatrick et al., 2015)
//
//   SST = a0 + a1 BT11 + a2 (BT11 - BT12) T_sfc
//         + a3 (sec(theta) - 1)(BT11 - BT12)
//         + a4 mirror + a5 theta* + a6 theta^2
// ---------------------------------------------------------------------------
struct NLSSTSample {
    double bt11_k = 0.0;
    double bt12_k = 0.0;
    double first_guess_k = 0.0;  // climatological or reference SST
    double sat_zenith_rad = 0.0;
    double mirror = 0.0;         // scan-mirror side indicator (0 or 1)
    double theta_star_rad = 0.0; // relative azimuth term
    double buoy_sst_k = 0.0;     // matchup truth
};

struct NLSSTModel {
    std::vector<double> coefficients;  // a0 .. a6
    double rmse = 0.0;
    double bias = 0.0;
    std::size_t samples = 0;
    bool calibrated = false;
};

std::vector<double> nlsst_features(double bt11_k, double bt12_k, double first_guess_k,
                                   double sat_zenith_rad, double mirror, double theta_star_rad);

NLSSTModel calibrate_nlsst(const std::vector<NLSSTSample>& samples);

double nlsst_retrieve(const NLSSTModel& model, double bt11_k, double bt12_k, double first_guess_k,
                      double sat_zenith_rad, double mirror = 0.0, double theta_star_rad = 0.0);

bool read_nlsst_samples(const std::string& path, std::vector<NLSSTSample>& out, std::string& error);
bool write_nlsst_model(const std::string& path, const NLSSTModel& model);
bool read_nlsst_model(const std::string& path, NLSSTModel& model, std::string& error);

// ---------------------------------------------------------------------------
// Vegetation greenness: NDVI (Rouse et al., 1974)
//   NDVI = (rho_NIR - rho_RED) / (rho_NIR + rho_RED)
// ---------------------------------------------------------------------------
double ndvi(double red_reflectance, double nir_reflectance);

// ---------------------------------------------------------------------------
// Sea ice: NASA Team algorithm (Cavalieri et al., 1984)
//   PR(19)        = (TB19V - TB19H) / (TB19V + TB19H)
//   GR(37V/19V)   = (TB37V - TB19V) / (TB37V + TB19V)
// The physical mixture model with open water, first-year ice and multiyear ice
// tie points is inverted directly; the empirical polarization/gradient form
//   C = (a0 + a1 PR + a2 GR + a3 PR GR) / (1 + d1 PR + d2 GR + d3 PR GR)
// is calibrated from the same tie points.
// ---------------------------------------------------------------------------
struct SeaIceTiePoints {
    // [constituent][channel], constituents open water, first-year, multiyear;
    // channels 19V, 19H, 37V (Kelvin).
    std::array<std::array<double, 3>, 3> tb{{
        {{170.0, 105.0, 195.0}},
        {{245.0, 235.0, 245.0}},
        {{215.0, 200.0, 200.0}},
    }};
};

struct SeaIceConcentrations {
    double open_water = 0.0;
    double first_year = 0.0;
    double multiyear = 0.0;

    double total() const { return first_year + multiyear; }
    double extent_indicator(double threshold = 0.15) const { return total() >= threshold ? 1.0 : 0.0; }
};

double polarization_ratio(double tb19v, double tb19h);
double gradient_ratio(double tb37v, double tb19v);

SeaIceConcentrations nasa_team_mixture_invert(const SeaIceTiePoints& tie_points, double tb19v,
                                              double tb19h, double tb37v);

struct PolarGradientModel {
    std::vector<double> cf;  // 4 coefficients
    std::vector<double> cm;  // 4 coefficients
    std::vector<double> d;   // d0 fixed to 1, d1 .. d3 calibrated
    double rmse_cf = 0.0;
    double rmse_cm = 0.0;
    std::size_t samples = 0;
    bool calibrated = false;
};

PolarGradientModel calibrate_polar_gradient(const SeaIceTiePoints& tie_points, int steps = 20);

SeaIceConcentrations nasa_team_polar_gradient(const PolarGradientModel& model, double tb19v,
                                              double tb19h, double tb37v);

bool read_sea_ice_tie_points(const std::string& path, SeaIceTiePoints& tie_points,
                             std::string& error);
bool write_polar_gradient_model(const std::string& path, const PolarGradientModel& model);
bool read_polar_gradient_model(const std::string& path, PolarGradientModel& model,
                               std::string& error);

struct SeaIceCell {
    double lat_min_deg = 0.0;
    double lat_max_deg = 0.0;
    double lon_min_deg = 0.0;
    double lon_max_deg = 0.0;
    double concentration = 0.0;
};

double spherical_cell_area_km2(double lat_min_deg, double lat_max_deg, double lon_min_deg,
                               double lon_max_deg);

struct SeaIceAreaResult {
    double extent_km2 = 0.0;
    double area_km2 = 0.0;
    double total_grid_km2 = 0.0;
    std::size_t cells_total = 0;
    std::size_t cells_ice_covered = 0;
};

SeaIceAreaResult sea_ice_extent_area(const std::vector<SeaIceCell>& cells,
                                     double threshold = 0.15);

bool read_sea_ice_grid(const std::string& path, std::vector<SeaIceCell>& cells,
                       std::string& error);

// ---------------------------------------------------------------------------
// Precipitation
//   GPI: P = 3 mm/h * fraction of IR pixels colder than 235 K
//   water-vapour adjustment: P = exp(a1 (w - a2)) * GPI
//   GPCP merge: P = sum(w_i P_i) / sum(w_i), w_i = 1 / sigma_i^2
// ---------------------------------------------------------------------------
double gpi(double cold_cloud_fraction);

double gpi_with_water_vapour(double gpi_value, double precipitable_water_cm, double a1, double a2);

struct PrecipEstimate {
    std::string source;
    double value = 0.0;
    double error_sd = 0.0;
};

struct MergedPrecip {
    double value = 0.0;
    std::vector<double> weights;
    std::size_t used = 0;
    bool valid = false;
};

MergedPrecip gpcp_merge(const std::vector<PrecipEstimate>& estimates);

bool read_precip_estimates(const std::string& path, std::vector<PrecipEstimate>& out,
                           std::string& error);

// ---------------------------------------------------------------------------
// Simple coupled-system diagnostics
// ---------------------------------------------------------------------------
// N = F + lambda * deltaT
double toa_energy_imbalance(double effective_forcing, double feedback_parameter,
                            double temperature_change);

// Rn = (1 - albedo) Rs + eps_a Rl - sigma eps_s Ts^4
double net_surface_radiation(double albedo, double shortwave_down, double emissivity_air,
                             double longwave_down, double emissivity_surface,
                             double surface_temperature);

// Residual of Rn = LE + H + G.
double surface_energy_residual(double net_radiation, double latent_heat, double sensible_heat,
                               double ground_heat);

// G_ATM = k * (M2 - M1) / (t2 - t1), k = 2.124 PgC per ppm
double carbon_growth_rate_pgc_per_year(double ppm_start, double ppm_end, double year_start,
                                       double year_end);

struct VelocitySample {
    double latitude_deg = 0.0;
    double longitude_deg = 0.0;
    double depth_m = 0.0;
    double velocity_m_s = 0.0;
};

// Zonally and vertically integrated meridional transport Q(y, z) = int_z^0 int v dx dz,
// expressed in Sverdrup (1 Sv = 1e6 m^3/s).
double amoc_streamfunction_sv(const std::vector<VelocitySample>& samples, double depth_m);

bool read_velocity_samples(const std::string& path, std::vector<VelocitySample>& out,
                           std::string& error);

}  // namespace climate
