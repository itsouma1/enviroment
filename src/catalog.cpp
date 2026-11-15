#include "catalog.hpp"

#include <algorithm>
#include <cctype>

namespace climate {
namespace {

std::string lowered(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (char c : text) out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    return out;
}

const std::vector<CatalogEntry> kEntries = {
    {"sst",
     "Sea surface temperature",
     "degC",
     {"temperature", "ocean"},
     {"global", "arctic", "antarctic", "tropics", "north_atlantic", "midlatitudes", "amazon", "sahel"},
     "Modified non-linear SST (NLSST) split-window retrieval from 11 and 12 micrometre brightness "
     "temperatures, scaled by a first-guess SST, with satellite zenith angle and mirror-side terms; "
     "coefficients regressed from satellite-buoy matchups.",
     "Monthly climatology removal then STL decomposition; Mann-Kendall and Sen's slope on the "
     "deseasonalised series.",
     "Skin sea surface temperature from MODIS/VIIRS thermal infrared; the primary ocean warming "
     "indicator and a boundary condition for sea ice and weather."},

    {"air_temperature",
     "Surface air temperature",
     "degC",
     {"temperature", "warming"},
     {"global", "arctic", "antarctic", "tropics", "north_atlantic", "midlatitudes", "amazon", "sahel"},
     "Station and reanalysis blend; anomalies are computed against a monthly climatology.",
     "Monthly anomalies, STL decomposition, Mann-Kendall and Sen's slope with autocorrelation "
     "correction.",
     "Near-surface air temperature over land and ocean; the headline indicator of surface warming."},

    {"co2",
     "Atmospheric carbon dioxide",
     "ppm",
     {"carbon", "greenhouse"},
     {"global", "arctic", "midlatitudes"},
     "Ground-based or satellite column CO2 mole fraction; the NOAA/CMDL curve fit separates the "
     "long-term trend from the seasonal cycle.",
     "Polynomial plus harmonic regression (Thoning et al., 1989) with iterative outlier rejection; "
     "growth rate from the analytic derivative of the fitted trend.",
     "Well-mixed greenhouse gas driving radiative forcing; Mauna Loa style records show the "
     "Keeling curve and its accelerating growth rate."},

    {"ch4",
     "Atmospheric methane",
     "ppb",
     {"carbon", "greenhouse"},
     {"global", "arctic", "tropics"},
     "Ground-based or satellite column CH4 mole fraction.",
     "Polynomial plus harmonic regression with iterative outlier rejection.",
     "Second most important well-mixed greenhouse gas; wetland and permafrost sources are a "
     "climate feedback."},

    {"ndvi",
     "Normalized Difference Vegetation Index",
     "unitless",
     {"vegetation", "drought"},
     {"global", "tropics", "amazon", "sahel", "midlatitudes", "arctic"},
     "NDVI = (rho_NIR - rho_RED) / (rho_NIR + rho_RED) from surface reflectances, composited to "
     "remove cloud and view-angle artefacts.",
     "Standardised monthly anomalies and STL decomposition; Mann-Kendall and Sen's slope on the "
     "growing-season or deseasonalised series.",
     "Vegetation greenness and canopy density; greening and browning trends track drought, "
     "land-use change and growing-season length."},

    {"sea_ice_extent",
     "Sea ice extent",
     "million km2",
     {"sea_ice", "cryosphere"},
     {"global", "arctic", "antarctic"},
     "Sum of grid-cell areas where NASA Team ice concentration is at least 15%.",
     "Monthly anomalies on the 15% extent series; Mann-Kendall and Sen's slope per decade, most "
     "informative for September (Arctic minimum) and February (Arctic maximum).",
     "The 15% concentration contour area; the standard indicator of polar sea ice loss."},

    {"sea_ice_concentration",
     "Sea ice concentration",
     "fraction",
     {"sea_ice", "cryosphere"},
     {"global", "arctic", "antarctic"},
     "NASA Team algorithm applied to passive-microwave brightness temperatures: polarization ratio "
     "PR(19) and gradient ratio GR(37V/19V) inverted for first-year and multiyear ice fractions.",
     "Area-weighted anomalies and trends; grid-cell aggregation before Mann-Kendall testing.",
     "Fraction of ocean covered by ice; drives albedo feedback and ice-volume change."},

    {"precipitation",
     "Precipitation rate",
     "mm/day",
     {"precipitation", "drought"},
     {"global", "tropics", "amazon", "sahel", "midlatitudes", "arctic"},
     "GPI cold-cloud index from infrared brightness temperatures, adjusted by precipitable water, "
     "then merged with passive-microwave and gauge estimates by inverse-error variance.",
     "Monthly anomalies, STL decomposition and Mann-Kendall/Sen's slope.",
     "Rainfall and snowfall; central to the hydrological cycle, drought and freshwater availability."},

    {"ocean_heat_content",
     "Ocean heat content",
     "ZJ",
     {"ocean", "warming"},
     {"global", "north_atlantic", "arctic", "tropics"},
     "Vertically integrated temperature anomaly from Argo and repeat hydrography.",
     "STL decomposition with a long trend window; Mann-Kendall and Sen's slope per decade.",
     "The dominant reservoir of the Earth energy imbalance, storing roughly 90% of excess heat."},

    {"toa_energy_imbalance",
     "Top-of-atmosphere energy imbalance",
     "W/m2",
     {"energy_balance", "warming"},
     {"global"},
     "Net shortwave minus longwave flux at the top of the atmosphere from CERES-class radiometry; "
     "N = F + lambda dT relates it to forcing and feedback.",
     "Annual-mean trend testing with Mann-Kendall and Sen's slope.",
     "The physical measure of ongoing global heating, currently about 0.5 to 0.9 W/m2."},

    {"snow_cover",
     "Snow cover extent",
     "million km2",
     {"cryosphere", "sea_ice"},
     {"global", "arctic", "midlatitudes"},
     "Fractional snow cover from visible and near-infrared reflectance, aggregated to monthly extent.",
     "Monthly anomalies by calendar month; Mann-Kendall and Sen's slope.",
     "Seasonal snow extent; an albedo-feedback indicator and a freshwater storage term."},

    {"ozone",
     "Total column ozone",
     "DU",
     {"ozone", "atmosphere"},
     {"global", "arctic", "antarctic", "midlatitudes"},
     "OMI total column ozone from backscattered ultraviolet radiance, with cloud and surface "
     "albedo corrections.",
     "Monthly anomalies; STL decomposition and Mann-Kendall/Sen's slope.",
     "Stratospheric ozone protects the biosphere from ultraviolet radiation; its recovery is a "
     "direct test of the Montreal Protocol."},

    {"aerosol_optical_depth",
     "Aerosol optical depth",
     "unitless",
     {"aerosol", "atmosphere"},
     {"global", "tropics", "midlatitudes", "sahel", "amazon"},
     "OMI UV aerosol index and optical depth from near-UV radiances.",
     "Monthly anomalies; STL decomposition and Mann-Kendall/Sen's slope.",
     "Aerosol loading affects radiation, clouds and air quality; trends mix emissions, transport "
     "and humidity."},

    {"soil_moisture",
     "Surface soil moisture",
     "m3/m3",
     {"hydrology", "drought"},
     {"global", "tropics", "midlatitudes", "sahel", "amazon"},
     "SMAP L-band brightness-temperature retrieval of surface soil moisture.",
     "Monthly anomalies; STL decomposition and Mann-Kendall/Sen's slope.",
     "Root-zone water availability; a drought and land-atmosphere coupling indicator."},

    {"ocean_color",
     "Ocean chlorophyll-a",
     "mg/m3",
     {"ocean", "biosphere"},
     {"global", "tropics", "north_atlantic", "midlatitudes", "arctic"},
     "MODIS ocean-colour chlorophyll-a from blue/green reflectance ratios.",
     "Monthly anomalies; STL decomposition and Mann-Kendall/Sen's slope.",
     "Phytoplankton biomass at the base of the marine food web and a component of the ocean "
     "carbon cycle."},

    {"sea_level",
     "Global mean sea level",
     "mm",
     {"ocean", "warming"},
     {"global", "tropics", "midlatitudes"},
     "Satellite altimetry and tide-gauge reconstruction, corrected for glacial isostatic adjustment.",
     "STL decomposition and Mann-Kendall/Sen's slope, with acceleration checked through the "
     "harmonic polynomial term.",
     "Thermal expansion plus land-ice loss; the integrative indicator of ocean and cryosphere change."},
};

const std::vector<std::string> kProcesses = {"temperature", "warming",  "sea_ice",  "cryosphere",
                                             "vegetation",  "drought",  "precipitation",
                                             "carbon",      "greenhouse", "ocean",
                                             "energy_balance", "climate_change", "ozone",
                                             "aerosol", "hydrology", "biosphere"};

const std::vector<std::string> kRegions = {"global", "arctic", "antarctic", "tropics",
                                           "north_atlantic", "midlatitudes", "amazon", "sahel"};

}  // namespace

AnalysisMethod catalog_analysis_method(const std::string& id) {
    const std::string target = lowered(id);
    // Well-mixed greenhouse gases carry a large seasonal cycle on top of a smooth,
    // accelerating trend; the NOAA/CMDL polynomial plus harmonic fit is their method,
    // so the trend is not a by-product of an STL decomposition.
    if (target == "co2" || target == "ch4") return AnalysisMethod::HarmonicFit;
    return AnalysisMethod::SeasonalStl;
}

const std::vector<CatalogEntry>& catalog_entries() { return kEntries; }

const CatalogEntry* find_catalog_entry(const std::string& id) {
    const std::string target = lowered(id);
    for (const auto& entry : kEntries)
        if (lowered(entry.id) == target) return &entry;
    return nullptr;
}

std::vector<std::string> catalog_processes() { return kProcesses; }

std::vector<std::string> catalog_regions() { return kRegions; }

std::string canonical_process(const std::string& name) {
    const std::string target = lowered(name);
    for (const auto& process : kProcesses)
        if (process == target) return process;
    return std::string();
}

std::string canonical_region(const std::string& name) {
    const std::string target = lowered(name);
    for (const auto& region : kRegions)
        if (region == target) return region;
    return std::string();
}

}  // namespace climate
