// Self-checking test program for the climate analysis tool. Build and run with
// `make test`. Any failure prints a diagnostic and the process exits non-zero.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "aerosol.hpp"
#include "analysis.hpp"
#include "assim.hpp"
#include "atmos.hpp"
#include "catalog.hpp"
#include "common.hpp"
#include "cryo.hpp"
#include "csv.hpp"
#include "earth.hpp"
#include "fft.hpp"
#include "land.hpp"
#include "landuse.hpp"
#include "detect.hpp"
#include "climateidx.hpp"
#include "attrib.hpp"
#include "ml.hpp"
#include "copula.hpp"
#include "hmm.hpp"
#include "regime.hpp"
#include "causal.hpp"
#include "magnitude.hpp"
#include "coherence.hpp"
#include "change.hpp"
#include "uncertainty.hpp"
#include "metrology.hpp"
#include "geostat.hpp"
#include "intercal.hpp"
#include "homogen.hpp"
#include "validation.hpp"
#include "regional.hpp"
#include "lagrangian.hpp"
#include "power.hpp"
#include "quantile.hpp"
#include "ensemble.hpp"
#include "flow.hpp"
#include "ocean2.hpp"
#include "kernel.hpp"
#include "spatial.hpp"
#include "harmonic.hpp"
#include "fetch.hpp"
#include "intent.hpp"
#include "json.hpp"
#include "linalg.hpp"
#include "netcdf.hpp"
#include "ocean.hpp"
#include "ozone.hpp"
#include "pde.hpp"
#include "physics.hpp"
#include "loess.hpp"
#include "report.hpp"
#include "retrieval.hpp"
#include "series.hpp"
#include "stats.hpp"
#include "stl.hpp"
#include "trend.hpp"

using namespace climate;

namespace {

int failures = 0;

void check(bool condition, const char* expression, const char* file, int line) {
    if (condition) return;
    std::cerr << "FAIL " << file << ":" << line << ": " << expression << "\n";
    ++failures;
}

void check_near(double actual, double expected, double tolerance, const char* expression,
                const char* file, int line) {
    if (std::isfinite(actual) && std::fabs(actual - expected) <= tolerance) return;
    std::cerr << "FAIL " << file << ":" << line << ": " << expression << " = " << actual
              << ", expected " << expected << " +/- " << tolerance << "\n";
    ++failures;
}

#define CHECK(expr) check((expr), #expr, __FILE__, __LINE__)
#define CHECK_NEAR(actual, expected, tolerance) \
    check_near((actual), (expected), (tolerance), #actual, __FILE__, __LINE__)

void test_linalg() {
    std::vector<std::vector<double>> a = {{2.0, 1.0}, {1.0, 3.0}};
    std::vector<double> x;
    CHECK(solve_linear_system(a, {3.0, 5.0}, x));
    CHECK_NEAR(x[0], 0.8, 1e-12);
    CHECK_NEAR(x[1], 1.4, 1e-12);

    // Least squares fit of y = 2 + 3t.
    std::vector<std::vector<double>> design;
    std::vector<double> response;
    for (int i = 0; i < 10; ++i) {
        design.push_back({1.0, static_cast<double>(i)});
        response.push_back(2.0 + 3.0 * i);
    }
    std::vector<double> coefficients;
    CHECK(least_squares(design, response, coefficients));
    CHECK_NEAR(coefficients[0], 2.0, 1e-9);
    CHECK_NEAR(coefficients[1], 3.0, 1e-9);
}

void test_stats() {
    std::vector<double> values = {1.0, 2.0, 3.0, 4.0};
    CHECK_NEAR(mean_value(values), 2.5, 1e-12);
    CHECK_NEAR(median_value(values), 2.5, 1e-12);
    CHECK_NEAR(quantile_value(values, 0.5), 2.5, 1e-12);
    CHECK_NEAR(stdev_value(values), 1.2909944487358056, 1e-12);
    CHECK_NEAR(normal_cdf(0.0), 0.5, 1e-12);
    CHECK_NEAR(normal_cdf(1.959963984540054), 0.975, 1e-10);
    CHECK_NEAR(normal_quantile(0.975), 1.959963984540054, 1e-9);
    CHECK_NEAR(normal_quantile(0.5), 0.0, 1e-12);

    std::vector<double> alternating = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
    CHECK_NEAR(autocorrelation(alternating, 1), -1.0, 1e-12);
}

void test_csv_and_json() {
    const std::string csv = "a,b,c\n1,\"two, with comma\",NA\n2,\"quote \"\"x\"\"\",4\n";
    const std::string path = "test_tmp.csv";
    write_text_file(path, csv);
    const Table table = read_csv(path);
    CHECK(table.header.size() == 3);
    CHECK(table.rows.size() == 2);
    CHECK(table.column("b") == 1);
    CHECK(table.cell(0, 1) == "two, with comma");
    CHECK(table.cell(1, 1) == "quote \"x\"");
    double number = 0.0;
    CHECK(!parse_number("NA", number));
    CHECK(!parse_number("", number));
    CHECK(parse_number("4.5", number));
    CHECK_NEAR(number, 4.5, 1e-12);
    std::remove(path.c_str());

    const std::string text = "{\"process\": \"warming\", \"region\": \"arctic\", "
                             "\"start_year\": 1979, \"end_year\": 2024, \"notes\": [\"a\", \"b\"]}";
    const json::Value root = json::parse(text);
    CHECK(root.is_object());
    CHECK(root.find("process") != nullptr);
    CHECK(root.find("process")->as_string("") == "warming");
    const json::Value* notes = root.find("notes");
    CHECK(notes != nullptr && notes->is_array() && notes->items.size() == 2);
    const json::Value reparsed = json::parse(json::dump(root, 2));
    CHECK(reparsed.find("region")->as_string("") == "arctic");
    CHECK(json::escape("line\nbreak\"quoted\"") == "\"line\\nbreak\\\"quoted\\\"\"");
}

void test_fft() {
    std::vector<Complex> data(12);
    for (int i = 0; i < 12; ++i) data[static_cast<std::size_t>(i)] = Complex(std::sin(0.7 * i), 0.0);
    const std::vector<Complex> original = data;
    fft(data, false);
    fft(data, true);
    for (std::size_t i = 0; i < data.size(); ++i)
        CHECK_NEAR(data[i].real(), original[i].real(), 1e-9);

    // A slow sine survives the low-pass filter; a fast one is removed.
    std::vector<double> slow(120, 0.0);
    std::vector<double> fast(120, 0.0);
    for (int i = 0; i < 120; ++i) {
        slow[static_cast<std::size_t>(i)] = std::sin(2.0 * PI * i / 40.0);
        fast[static_cast<std::size_t>(i)] = std::sin(2.0 * PI * i / 3.0);
    }
    const std::vector<double> slow_filtered = lowpass_filter(slow, 10.0);
    const std::vector<double> fast_filtered = lowpass_filter(fast, 10.0);
    double slow_error = 0.0;
    double fast_amplitude = 0.0;
    for (std::size_t i = 20; i < 100; ++i) {
        slow_error = std::max(slow_error, std::fabs(slow_filtered[i] - slow[i]));
        fast_amplitude = std::max(fast_amplitude, std::fabs(fast_filtered[i]));
    }
    CHECK(slow_error < 0.05);
    CHECK(fast_amplitude < 0.05);
}

void test_loess_and_stl() {
    std::vector<double> linear(80, 0.0);
    for (std::size_t i = 0; i < linear.size(); ++i) linear[i] = 3.0 + 0.5 * static_cast<double>(i);
    const std::vector<double> smoothed = loess_smooth(linear, 21, 1);
    for (std::size_t i = 0; i < linear.size(); ++i)
        CHECK_NEAR(smoothed[i], linear[i], 1e-8);

    const int period = 12;
    const int n = 240;
    std::vector<double> series(n, 0.0);
    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / 12.0;
        series[static_cast<std::size_t>(i)] =
            10.0 + 0.5 * t + 3.0 * std::sin(2.0 * PI * i / period);
    }
    STLParams params;
    params.period = period;
    params.seasonal_window = 13;
    params.trend_window = 25;
    params.low_pass_window = 13;
    const STLResult result = stl_decompose(series, period, params);
    double seasonal_error = 0.0;
    for (int i = 24; i < n - 24; ++i)
        seasonal_error = std::max(seasonal_error,
                                  std::fabs(result.seasonal[static_cast<std::size_t>(i)] -
                                            3.0 * std::sin(2.0 * PI * i / period)));
    CHECK(seasonal_error < 0.25);
    CHECK_NEAR(result.trend[120], 10.0 + 0.5 * 10.0, 0.15);
}

void test_harmonic() {
    const int n = 480;
    std::vector<double> time(n);
    std::vector<double> values(n);
    for (int i = 0; i < n; ++i) {
        const double t = 1980.0 + i / 12.0;
        time[static_cast<std::size_t>(i)] = t;
        values[static_cast<std::size_t>(i)] =
            340.0 + 1.5 * (t - 2000.0) + 0.01 * (t - 2000.0) * (t - 2000.0) +
            3.0 * std::cos(2.0 * PI * (t - 1980.0)) - 2.0 * std::sin(2.0 * PI * (t - 1980.0));
    }
    const HarmonicModel model = fit_harmonic(time, values, 2, 4, 1.0, 0, 3.0);
    double max_error = 0.0;
    for (int i = 0; i < n; ++i)
        max_error = std::max(max_error, std::fabs(harmonic_evaluate(model, time[static_cast<std::size_t>(i)]) -
                                                  values[static_cast<std::size_t>(i)]));
    CHECK(max_error < 1e-6);
    CHECK_NEAR(model.r_squared, 1.0, 1e-9);
    CHECK_NEAR(harmonic_seasonal_amplitude(model), 2.0 * std::sqrt(13.0), 1e-3);

    const std::vector<HarmonicComponent> components = harmonic_components(model);
    CHECK(!components.empty());
    CHECK_NEAR(components[0].amplitude, std::sqrt(13.0), 1e-6);

    const double t0 = 2005.3;
    const double step = 1e-4;
    const double numeric_derivative =
        (harmonic_evaluate(model, t0 + step) - harmonic_evaluate(model, t0 - step)) / (2.0 * step);
    CHECK_NEAR(harmonic_derivative(model, t0), numeric_derivative, 1e-5);
    CHECK_NEAR(harmonic_trend_derivative(model, t0), 1.5 + 0.02 * (t0 - 2000.0), 1e-6);
}

void test_trend() {
    std::vector<double> increasing(50);
    for (std::size_t i = 0; i < increasing.size(); ++i) increasing[i] = 0.3 * i;
    const MKResult mk = mann_kendall(increasing);
    CHECK(mk.s > 0.0);
    CHECK(mk.z > 0.0);
    CHECK(mk.significant);

    std::vector<double> constant(50, 4.0);
    const MKResult flat = mann_kendall(constant);
    CHECK_NEAR(flat.s, 0.0, 1e-12);
    CHECK_NEAR(flat.p_value, 1.0, 1e-12);

    std::vector<double> time(60);
    std::vector<double> linear(60);
    for (std::size_t i = 0; i < 60; ++i) {
        time[i] = static_cast<double>(i) / 12.0;
        linear[i] = 3.0 + 2.0 * time[i];
    }
    const SenResult sen = sens_slope(time, linear);
    CHECK(sen.computed);
    CHECK_NEAR(sen.slope, 2.0, 1e-12);
    CHECK(sen.lower <= sen.slope && sen.slope <= sen.upper);

    // Prewhitening must not destroy a strong linear trend.
    const std::vector<double> prewhitened = trend_free_prewhitening(time, linear);
    const MKResult mk_prewhitened = mann_kendall(prewhitened);
    CHECK(mk_prewhitened.significant);
    CHECK(mk_prewhitened.z > 0.0);
}

void test_retrieval() {
    const double factor = nlsst_features(300.0, 299.1, 300.5, 0.4, 0.0, 0.2)[2];
    CHECK(std::isfinite(factor));

    // Calibration must recover the generating coefficients.
    const double truth[7] = {1.2, 0.995, 0.00035, 0.02, 0.01, -0.05, 0.04};
    std::vector<NLSSTSample> samples;
    for (int i = 0; i < 600; ++i) {
        NLSSTSample sample;
        // The predictors vary independently, so the calibration problem is
        // well conditioned and the generating coefficients are identifiable.
        sample.bt11_k = 298.0 + 0.01 * i;
        sample.bt12_k = sample.bt11_k - 0.9 - 0.0005 * i;
        sample.first_guess_k = 300.0 + 0.4 * std::sin(0.011 * i) + 0.002 * (i % 17);
        sample.sat_zenith_rad = deg_to_rad(5.0 + 0.05 * (i % 60));
        sample.mirror = (i % 2 == 0) ? 0.0 : 1.0;
        sample.theta_star_rad = deg_to_rad(120.0 * std::sin(0.017 * i));
        const std::vector<double> features =
            nlsst_features(sample.bt11_k, sample.bt12_k, sample.first_guess_k,
                           sample.sat_zenith_rad, sample.mirror, sample.theta_star_rad);
        double buoy = 0.0;
        for (std::size_t j = 0; j < 7; ++j) buoy += truth[j] * features[j];
        sample.buoy_sst_k = buoy;
        samples.push_back(sample);
    }
    const NLSSTModel model = calibrate_nlsst(samples);
    CHECK_NEAR(model.coefficients[1], truth[1], 1e-4);
    CHECK_NEAR(model.coefficients[3], truth[3], 1e-4);
    CHECK_NEAR(model.coefficients[4], truth[4], 1e-4);
    CHECK(model.rmse < 1e-4);

    CHECK_NEAR(ndvi(0.1, 0.5), 0.6666666666666666, 1e-12);
    CHECK_NEAR(ndvi(0.5, 0.1), -0.6666666666666666, 1e-12);
    CHECK(!std::isfinite(ndvi(0.0, 0.0)));  // zero denominator yields a non-finite value

    SeaIceTiePoints tie;
    SeaIceConcentrations open = nasa_team_mixture_invert(tie, tie.tb[0][0], tie.tb[0][1], tie.tb[0][2]);
    CHECK_NEAR(open.total(), 0.0, 1e-9);
    SeaIceConcentrations first = nasa_team_mixture_invert(tie, tie.tb[1][0], tie.tb[1][1], tie.tb[1][2]);
    CHECK_NEAR(first.total(), 1.0, 1e-9);
    CHECK_NEAR(first.first_year, 1.0, 1e-9);

    const PolarGradientModel empirical = calibrate_polar_gradient(tie, 10);
    const SeaIceConcentrations empirical_first =
        nasa_team_polar_gradient(empirical, tie.tb[1][0], tie.tb[1][1], tie.tb[1][2]);
    CHECK(empirical_first.total() > 0.9);
    const SeaIceConcentrations empirical_open =
        nasa_team_polar_gradient(empirical, tie.tb[0][0], tie.tb[0][1], tie.tb[0][2]);
    CHECK(empirical_open.total() < 0.1);

    CHECK_NEAR(gpi(0.5), 1.5, 1e-12);
    CHECK_NEAR(gpi(2.0), 3.0, 1e-12);  // fraction is clamped
    std::vector<PrecipEstimate> estimates = {{"a", 2.0, 1.0}, {"b", 4.0, 1.0}};
    const MergedPrecip merged = gpcp_merge(estimates);
    CHECK(merged.valid);
    CHECK_NEAR(merged.value, 3.0, 1e-12);
    estimates[1].error_sd = 2.0;
    CHECK_NEAR(gpcp_merge(estimates).value, 2.4, 1e-12);
    (void)gpi_with_water_vapour(1.0, 3.0, 0.1, 3.0);

    CHECK_NEAR(carbon_growth_rate_pgc_per_year(400.0, 410.0, 2000.0, 2010.0), 2.124, 1e-12);
    CHECK_NEAR(net_surface_radiation(0.2, 200.0, 0.9, 300.0, 0.98, 288.0),
               (1.0 - 0.2) * 200.0 + 0.9 * 300.0 - STEFAN_BOLTZMANN * 0.98 * std::pow(288.0, 4.0),
               1e-6);
    CHECK_NEAR(toa_energy_imbalance(2.0, -1.0, 1.0), 1.0, 1e-12);

    // Total spherical area equals 4 pi R^2.
    double total_area = 0.0;
    for (int lat = -90; lat < 90; ++lat)
        total_area += spherical_cell_area_km2(lat, lat + 1, 0.0, 360.0);
    CHECK_NEAR(total_area, 4.0 * PI * EARTH_RADIUS_KM * EARTH_RADIUS_KM, 1e6);
}

void test_series_and_analysis() {
    std::vector<double> values = {1.0, std::nan(""), 3.0, std::nan(""), 5.0};
    std::size_t filled = 0;
    const std::vector<double> filled_values = interpolate_missing(values, &filled);
    CHECK(filled == 2);
    CHECK_NEAR(filled_values[1], 2.0, 1e-12);
    CHECK_NEAR(filled_values[3], 4.0, 1e-12);

    const Climatology climatology = phase_climatology({1.0, 2.0, 3.0, 2.0, 3.0, 4.0}, 2);
    CHECK_NEAR(climatology.mean[0], 7.0 / 3.0, 1e-12);
    CHECK_NEAR(climatology.mean[1], 8.0 / 3.0, 1e-12);
    CHECK(climatology.count[0] == 3 && climatology.count[1] == 3);

    // Synthetic monthly series: +0.02 per year trend and a 0.5 seasonal cycle.
    const int n = 480;
    std::vector<double> time(n);
    std::vector<double> series(n);
    for (int i = 0; i < n; ++i) {
        const double t = 1980.0 + i / 12.0;
        time[static_cast<std::size_t>(i)] = t;
        series[static_cast<std::size_t>(i)] =
            5.0 + 0.02 * (t - 1980.0) + 0.5 * std::sin(2.0 * PI * i / 12.0);
    }
    SeriesOptions options;
    options.period = 12;
    options.seasonal_window = 13;
    options.trend_window = 25;
    options.low_pass_window = 13;
    SeriesAnalysis analysis = analyze_series(time, series, options);
    CHECK(analysis.decomposed);
    CHECK(analysis.mk.significant);
    CHECK(analysis.mk.z > 0.0);
    CHECK_NEAR(analysis.sen.slope, 0.02, 1e-6);
    CHECK_NEAR(analysis.sen_per_decade, 0.2, 1e-6);
    CHECK_NEAR(analysis.seasonal_amplitude, 1.0, 0.08);

    const std::string json_report = series_report_json(analysis, "test", "K", "memory");
    CHECK(json::parse(json_report).find("trend") != nullptr);
}

void test_catalog_and_intent() {
    CHECK(find_catalog_entry("sst") != nullptr);
    CHECK(find_catalog_entry("not_a_variable") == nullptr);
    CHECK(canonical_process("sea ice") == "");
    CHECK(canonical_process("sea_ice") == "sea_ice");
    CHECK(canonical_region("Arctic") == "arctic");

    const Intent sea_ice = parse_intent_text("How has arctic sea ice changed since 1979?");
    CHECK(sea_ice.valid);
    CHECK(sea_ice.process == "sea_ice");
    CHECK(sea_ice.region == "arctic");
    CHECK(sea_ice.start_year == 1979);
    CHECK(sea_ice.end_year == 2025);
    Selection selection = select_variables(sea_ice);
    CHECK(!selection.entries.empty());
    bool found_extent = false;
    for (const auto& entry : selection.entries)
        if (entry.id == "sea_ice_extent") found_extent = true;
    CHECK(found_extent);

    const Intent warming = parse_intent_text("is the arctic warming?", 2025);
    CHECK(warming.valid);
    CHECK(warming.process == "warming");
    CHECK(warming.region == "arctic");

    const Intent changing = parse_intent_text("is the arctic changing?", 2025);
    CHECK(changing.valid);
    CHECK(changing.process == "climate_change");
    CHECK(changing.region == "arctic");
    const Selection broad = select_variables(changing);
    bool broad_has_ice = false;
    bool broad_has_sst = false;
    for (const auto& entry : broad.entries) {
        if (entry.id == "sea_ice_extent") broad_has_ice = true;
        if (entry.id == "sst") broad_has_sst = true;
    }
    CHECK(broad_has_ice);
    CHECK(broad_has_sst);

    const Intent precipitation = parse_intent_text("rainfall in the sahel over the last 30 years");
    CHECK(precipitation.valid);
    CHECK(precipitation.process == "precipitation");
    CHECK(precipitation.region == "sahel");
    CHECK(precipitation.start_year == 1995);
    CHECK(precipitation.end_year == 2025);

    CHECK(!parse_intent_text("what is the capital of France?").valid);

    const Intent from_json = parse_intent_json(
        "{\"process\": \"carbon\", \"region\": \"global\", \"start_year\": 1980, \"end_year\": 2024}");
    CHECK(from_json.valid);
    CHECK(from_json.process == "carbon");
    CHECK(select_variables(from_json).entries.size() == 2);

    CHECK(!parse_intent_json("{\"process\": \"volcanoes\"}").valid);
    CHECK(!parse_intent_json("{\"variable\": \"magic\"}").valid);

    const Intent explicit_variable = parse_intent_json("{\"variable\": \"ndvi\", \"region\": \"amazon\"}");
    CHECK(explicit_variable.valid);
    CHECK(select_variables(explicit_variable).entries.size() == 1);

    const std::string ir = intent_to_json(precipitation);
    CHECK(json::parse(ir).find("process")->as_string("") == "precipitation");
}

void test_co2_analysis() {
    const int n = 540;
    std::vector<double> time(n);
    std::vector<double> values(n);
    for (int i = 0; i < n; ++i) {
        const double t = 1980.0 + i / 12.0;
        time[static_cast<std::size_t>(i)] = t;
        const double dt = t - 1980.0;
        values[static_cast<std::size_t>(i)] =
            338.0 + 1.55 * dt + 0.012 * dt * dt + 3.0 * std::sin(2.0 * PI * t) +
            1.0 * std::sin(4.0 * PI * t);
    }
    const CO2Analysis analysis = analyze_co2(time, values);
    CHECK(analysis.model.r_squared > 0.999);
    CHECK(analysis.seasonal_amplitude > 6.5 && analysis.seasonal_amplitude < 7.5);
    CHECK(analysis.growth_rate_mean > 1.5);
    CHECK(analysis.growth_rate_end > analysis.growth_rate_start);
    CHECK(analysis.mk.significant);
    CHECK(analysis.growth_rate_end_pgc_per_year > 0.0);
}

void test_amoc() {
    std::vector<VelocitySample> samples;
    for (int depth = 0; depth <= 1000; depth += 100) {
        for (int lon = 0; lon <= 10; lon += 1) {
            VelocitySample sample;
            sample.latitude_deg = 0.0;
            sample.longitude_deg = static_cast<double>(lon);
            sample.depth_m = static_cast<double>(depth);
            sample.velocity_m_s = 1.0;
            samples.push_back(sample);
        }
    }
    const double transport = amoc_streamfunction_sv(samples, 0.0);
    CHECK(transport > 0.0);
    const double rough = 1.0 * (10.0 * PI / 180.0) * EARTH_RADIUS_KM * 1000.0 * 1000.0 / 1.0e6;
    CHECK_NEAR(transport, rough, 0.1 * rough);
}

void test_netcdf() {
    // Build a minimal CDF-1 file: time(3), lat(2), and sst(time, lat).
    struct Writer {
        std::vector<unsigned char> bytes;
        void u32(std::uint32_t value) {
            for (int shift = 24; shift >= 0; shift -= 8)
                bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xFF));
        }
        void name(const std::string& text) {
            u32(static_cast<std::uint32_t>(text.size()));
            for (char c : text) bytes.push_back(static_cast<unsigned char>(c));
            while (bytes.size() % 4 != 0) bytes.push_back(0);
        }
        void f64(double value) {
            std::uint64_t bits = 0;
            std::memcpy(&bits, &value, sizeof(bits));
            for (int shift = 56; shift >= 0; shift -= 8)
                bytes.push_back(static_cast<unsigned char>((bits >> shift) & 0xFF));
        }
        std::size_t placeholder() {
            const std::size_t position = bytes.size();
            u32(0);
            return position;
        }
        void patch(std::size_t position, std::uint32_t value) {
            bytes[position] = static_cast<unsigned char>((value >> 24) & 0xFF);
            bytes[position + 1] = static_cast<unsigned char>((value >> 16) & 0xFF);
            bytes[position + 2] = static_cast<unsigned char>((value >> 8) & 0xFF);
            bytes[position + 3] = static_cast<unsigned char>(value & 0xFF);
        }
    };

    Writer writer;
    writer.bytes = {static_cast<unsigned char>('C'), static_cast<unsigned char>('D'),
                    static_cast<unsigned char>('F'), 1};
    writer.u32(0);  // numrecs
    writer.u32(10);
    writer.u32(2);  // two dimensions
    writer.name("time");
    writer.u32(3);
    writer.name("lat");
    writer.u32(2);
    writer.u32(0);
    writer.u32(0);  // no global attributes
    writer.u32(11);
    writer.u32(3);  // three variables

    writer.name("time");
    writer.u32(1);
    writer.u32(0);
    writer.u32(0);
    writer.u32(0);
    writer.u32(6);
    writer.u32(24);
    const std::size_t time_begin = writer.placeholder();

    writer.name("lat");
    writer.u32(1);
    writer.u32(1);
    writer.u32(0);
    writer.u32(0);
    writer.u32(6);
    writer.u32(16);
    const std::size_t lat_begin = writer.placeholder();

    writer.name("sst");
    writer.u32(2);
    writer.u32(0);
    writer.u32(1);
    writer.u32(12);
    writer.u32(1);  // one attribute
    writer.name("units");
    writer.u32(2);
    writer.u32(1);
    writer.bytes.push_back(static_cast<unsigned char>('K'));
    while (writer.bytes.size() % 4 != 0) writer.bytes.push_back(0);
    writer.u32(6);
    writer.u32(48);
    const std::size_t sst_begin = writer.placeholder();

    while (writer.bytes.size() % 4 != 0) writer.bytes.push_back(0);
    const std::uint32_t header_size = static_cast<std::uint32_t>(writer.bytes.size());
    writer.patch(time_begin, header_size);
    writer.patch(lat_begin, header_size + 24);
    writer.patch(sst_begin, header_size + 40);

    for (int i = 0; i < 3; ++i) writer.f64(i);
    writer.f64(10.0);
    writer.f64(20.0);
    for (int i = 1; i <= 6; ++i) writer.f64(i);

    const std::string path = "test_climate.nc";
    {
        std::ofstream output(path, std::ios::binary);
        output.write(reinterpret_cast<const char*>(writer.bytes.data()),
                     static_cast<std::streamsize>(writer.bytes.size()));
    }

    CHECK(!is_hdf5_file(path));
    const NetCDFFile file = NetCDFFile::open(path);
    CHECK(file.version() == 1);
    CHECK(file.dimensions().size() == 2);
    CHECK(file.dimensions()[0].name == "time");
    CHECK(file.dimensions()[0].length == 3);
    CHECK(file.variables().size() == 3);

    const NcVariable* sst = file.find_variable("sst");
    CHECK(sst != nullptr);
    CHECK(sst->dimension_ids.size() == 2);
    CHECK(sst->shape.size() == 2 && sst->shape[0] == 3 && sst->shape[1] == 2);
    CHECK(file.attribute_text(*sst, "units") == "K");

    const std::vector<double> values = file.read_variable("sst");
    CHECK(values.size() == 6);
    for (std::size_t i = 0; i < values.size(); ++i) CHECK_NEAR(values[i], static_cast<double>(i + 1), 1e-12);

    CHECK_NEAR(decimal_year_from_time(0.0, "days since 2000-01-01"), 2000.0, 1e-9);
    CHECK_NEAR(decimal_year_from_time(366.0, "days since 2000-01-01"), 2001.0, 1e-6);
    CHECK_NEAR(decimal_year_from_time(12.0, "months since 2000-01-01"), 2001.0, 1e-9);

    const NetCDFExportResult exported =
        export_netcdf_to_csv(path, "sst", "timeseries", "test_climate.csv");
    CHECK(exported.ok);
    if (exported.ok) {
        const Table table = read_csv("test_climate.csv");
        CHECK(table.rows.size() == 3);
        CHECK(table.column("value") == 1);
        double first = 0.0;
        CHECK(parse_number(table.cell(0, 1), first));
        const double weight_low = std::cos(deg_to_rad(10.0));
        const double weight_high = std::cos(deg_to_rad(20.0));
        CHECK_NEAR(first, (1.0 * weight_low + 2.0 * weight_high) / (weight_low + weight_high), 1e-9);
    }

    const NetCDFExportResult bad_grid =
        export_netcdf_to_csv(path, "sst", "grid", "test_climate_grid.csv");
    CHECK(!bad_grid.ok);
    CHECK(!bad_grid.error.empty());

    std::remove(path.c_str());
    std::remove("test_climate.csv");
}

void test_physics() {
    CHECK_NEAR(equilibrium_temperature_change(3.7, -1.0), 3.7, 1e-12);
    CHECK_NEAR(climate_sensitivity_2x(-2.0), 1.85, 1e-12);
    CHECK_NEAR(co2_radiative_forcing(560.0, 280.0), 5.35 * std::log(2.0), 1e-12);
    CHECK_NEAR(co2_radiative_forcing(400.0, 400.0), 0.0, 1e-12);

    FeedbackSet feedback;
    feedback.planck = -3.2;
    feedback.water_vapour = -1.8;
    feedback.lapse_rate = 0.0;
    feedback.albedo = 0.35;
    feedback.cloud = 0.6;
    CHECK_NEAR(feedback.total(), -4.05, 1e-12);

    const FeedbackSet defaults;
    CHECK(defaults.total() < -1.0 && defaults.total() > -1.5);
    CHECK(climate_sensitivity_2x(defaults.total()) > 2.5 &&
          climate_sensitivity_2x(defaults.total()) < 3.5);

    // Mie theory in the Rayleigh limit: Qsca = (8/3) x^4 |K|^2.
    const double x = 0.01;
    const double m = 1.33;
    const double k = (m * m - 1.0) / (m * m + 2.0);
    const double expected = (8.0 / 3.0) * std::pow(x, 4) * k * k;
    const MieResult rayleigh = mie_sphere(x, std::complex<double>(m, 0.0));
    CHECK(std::fabs(rayleigh.qsca - expected) / expected < 1e-3);
    CHECK(std::fabs(rayleigh.qext - rayleigh.qsca) / rayleigh.qsca < 1e-3);
    CHECK_NEAR(rayleigh.single_scattering_albedo, 1.0, 1e-6);

    const MieResult absorbing = mie_sphere(5.0, std::complex<double>(1.33, 0.01));
    CHECK(absorbing.qext > 0.0);
    CHECK(absorbing.qsca > 0.0);
    CHECK(absorbing.qsca <= absorbing.qext + 1e-12);
    CHECK(absorbing.single_scattering_albedo > 0.0 && absorbing.single_scattering_albedo <= 1.0);
    CHECK(absorbing.g > 0.0 && absorbing.g < 1.0);
    // Cross-checked against an independent Mie implementation.
    CHECK_NEAR(absorbing.qext, 3.484147367, 1e-6);
    CHECK_NEAR(absorbing.qsca, 3.288906472, 1e-6);
    CHECK_NEAR(absorbing.g, 0.640422993, 1e-5);
    CHECK_NEAR(absorbing.single_scattering_albedo, 0.943963078, 1e-6);
    const MieResult small_absorbing = mie_sphere(1.0, std::complex<double>(1.33, 0.01));
    CHECK_NEAR(small_absorbing.qext, 0.121821799, 1e-7);
    CHECK_NEAR(small_absorbing.qsca, 0.093317470, 1e-7);

    const WindStress stress = bulk_wind_stress(1.225, 1.3e-3, 10.0, 0.0);
    CHECK_NEAR(stress.tau_x, 0.15925, 1e-9);
    CHECK_NEAR(stress.tau_y, 0.0, 1e-12);
    CHECK_NEAR(stress.magnitude(), 0.15925, 1e-9);

    const double f = coriolis_parameter(30.0);
    CHECK_NEAR(f, 7.2921159e-5, 1e-15);
    const EkmanTransport ekman = ekman_transport(0.1, 0.0, 1027.0, f);
    CHECK_NEAR(ekman.m_y, -0.1 / (1027.0 * f), 1e-12);
    CHECK_NEAR(ekman.m_x, 0.0, 1e-15);

    const GeostrophicVelocity velocity = geostrophic_velocity(1025.0, 1e-4, 0.0, 30.0);
    CHECK_NEAR(velocity.v, 1e-4 / (1025.0 * f), 1e-12);
    CHECK_NEAR(velocity.u, 0.0, 1e-15);

    const ThermalWind shear = thermal_wind(1025.0, 0.0, 1e-6, 30.0);
    CHECK(shear.du_dz > 0.0);
    CHECK_NEAR(shear.du_dz, 9.80665 * 1e-6 / (f * 1025.0), 1e-12);
    CHECK_NEAR(shear.dv_dz, 0.0, 1e-15);

    CHECK_NEAR(savi(0.1, 0.5, 0.5), (0.5 - 0.1) / (0.5 + 0.1 + 0.5) * 1.5, 1e-12);
    CHECK_NEAR(evi(0.1, 0.5, 0.05), 2.5 * (0.5 - 0.1) / (0.5 + 6.0 * 0.1 - 7.5 * 0.05 + 1.0),
               1e-12);
    CHECK(!std::isfinite(savi(0.0, 0.0, 0.0)));
}

void test_fetch() {
    CHECK(find_remote_dataset("sst") != nullptr);
    CHECK(find_remote_dataset("sea_ice_extent") != nullptr);
    CHECK(find_remote_dataset("not_a_dataset") == nullptr);
    CHECK(shell_quote("plain") == "'plain'");
    CHECK(shell_quote("it's") == "'it'\\''s'");
    CHECK(build_command({"python3", "tools/a b.py", "--x", "y"}) ==
          "'python3' 'tools/a b.py' '--x' 'y'");
    CHECK(!resolve_fetch_wrapper().empty());
}

void test_earth_land_spatial() {
    // Planck feedback and pattern-corrected ECS.
    CHECK_NEAR(planck_feedback(255.0), -3.7609, 1e-3);
    CHECK_NEAR(pattern_corrected_ecs(-1.3, 0.5), 4.625, 1e-12);

    // Two-layer energy balance model.
    TwoLayerModel model;
    CHECK_NEAR(two_layer_equilibrium(model, 3.7), 2.96, 1e-12);
    const TwoLayerState state = integrate_two_layer(model, 3.7, 20000.0, 2.0);
    CHECK_NEAR(state.surface, 2.96, 1e-3);
    CHECK_NEAR(state.deep, 2.96, 1e-3);
    CHECK_NEAR(ocean_heat_uptake(model, state), 0.0, 1e-5);
    CHECK_NEAR(transient_climate_response(model), 1.8974358974, 1e-9);
    CHECK_NEAR(equilibration_ratio(model), 1.5600000000, 1e-6);

    // Carbon-cycle feedbacks.
    CarbonFeedback carbon;
    carbon.beta_land = 0.3;
    carbon.gamma_land = -0.02;
    carbon.beta_ocean = 0.6;
    carbon.gamma_ocean = -0.05;
    CHECK_NEAR(land_carbon_change(carbon, 100.0, 2.0), 29.96, 1e-12);
    CHECK_NEAR(ocean_carbon_change(carbon, 100.0, 2.0), 59.9, 1e-12);
    CHECK_NEAR(airborne_fraction(2.5, 10.0), 0.531, 1e-12);

    // Permafrost decay with Q10 temperature dependence.
    CHECK_NEAR(permafrost_decay_rate(0.01, 2.0, 10.0), 0.02, 1e-15);
    CHECK_NEAR(permafrost_remaining_fraction(0.01, 2.0, 10.0, 100.0), 0.1353352832,
               1e-10);

    // Planck function inversion and Schwarzschild transfer.
    const double wavelength = 10e-6;
    const double radiance = planck_radiance(wavelength, 300.0);
    CHECK_NEAR(radiance, 9.9240e6, 5e2);
    CHECK_NEAR(brightness_temperature(wavelength, radiance), 300.0, 1e-9);
    CHECK_NEAR(optical_depth(0.5, 2.0, 1.0), 1.0, 1e-15);
    CHECK_NEAR(transmittance(1.0), std::exp(-1.0), 1e-15);
    CHECK_NEAR(weighting_function(1.0, 0.01, 0.0), 0.01 / 9.80665, 1e-15);
    CHECK_NEAR(layer_radiance(10.0, 20.0, std::log(2.0)), 15.0, 1e-12);

    // Two-stream: pure absorption, then a conservative scattering layer.
    const TwoStreamResult absorbed = two_stream_layer(1.0, 0.0, 0.0, 60.0, 0.0);
    CHECK_NEAR(absorbed.reflectance, 0.0, 1e-15);
    CHECK_NEAR(absorbed.direct_transmittance, std::exp(-2.0), 1e-12);
    CHECK_NEAR(absorbed.total_transmittance, std::exp(-2.0), 1e-12);
    CHECK_NEAR(absorbed.absorptance, 1.0 - std::exp(-2.0), 1e-12);

    const TwoStreamResult scattering = two_stream_layer(5.0, 1.0, 0.85, 60.0, 0.0);
    CHECK(scattering.reflectance > 0.0);
    CHECK_NEAR(scattering.reflectance + scattering.total_transmittance + scattering.absorptance,
               1.0, 1e-6);
    const TwoStreamResult reflected = two_stream_layer(5.0, 1.0, 0.85, 60.0, 0.3);
    CHECK_NEAR(reflected.reflectance + reflected.total_transmittance + reflected.absorptance, 1.0,
               1e-6);

    // Waves.
    CHECK_NEAR(shallow_water_speed(9.80665, 100.0), 31.3156, 1e-4);
    const EquatorialWave kelvin = kelvin_wave(2e-7, 9.80665, 100.0);
    CHECK_NEAR(kelvin.frequency, 6.2631e-6, 1e-9);
    const std::vector<double> frequencies = equatorial_wave_frequencies(2e-7, 2.3e-11, 30.0, 0);
    CHECK(!frequencies.empty());
    for (const double omega : frequencies) {
        const double residual = omega * omega * omega / (30.0 * 30.0) -
                                omega * (4e-14 + 2.3e-11 / 30.0) - 2.3e-11 * 2e-7;
        CHECK_NEAR(residual, 0.0, 1e-16);
    }

    // Turbulence.
    CHECK_NEAR(refractive_index_spectrum(1e-14, 1.0), 3.3e-16, 1e-20);
    CHECK_NEAR(dissipation_rate_from_spectrum(1.0, 1.0, 10.0), 1.540408, 1e-5);
    CHECK_NEAR(kolmogorov_microscale(1.5e-5, 1e-3), 1.3554e-3, 1e-6);
    CHECK(kolmogorov_spectrum(1.0, 1e-3, 100.0) > 0.0);
    CHECK_NEAR(stable_boundary_layer_ratio(1.0), 1.6129032258, 1e-9);

    // Land surface energy balance and Penman-Monteith.
    CHECK_NEAR(land_net_radiation(0.2, 1000.0, 1.0, 1.0, 288.0, 288.0), 800.0, 1e-9);
    CHECK_NEAR(sensible_heat_flux(1.2, 1005.0, 290.0, 288.0, 20.0), 120.6, 1e-9);
    CHECK_NEAR(aerodynamic_resistance(0.4, 0.3, 10.0, 0.0, 0.01, 0.0), 57.5646, 1e-3);
    const double latent = penman_monteith_latent_heat(1.0, 100.0, 1.225, 1005.0, 0.0, 50.0,
                                                      0.066, 50.0);
    const double sensible = penman_monteith_sensible_heat(1.0, 100.0, 1.225, 1005.0, 0.0, 50.0,
                                                          0.066, 50.0);
    CHECK_NEAR(latent, 88.3392, 1e-3);
    CHECK_NEAR(latent + sensible, 100.0, 1e-9);
    CHECK_NEAR(surface_resistance_soil(100.0, 0.4, 0.1, 0.3), 337.5, 1e-9);

    // Photosynthesis and stomatal conductance.
    CHECK_NEAR(farquhar_rubisco(100.0, 30.0, 4.0, 40.0, 21.0, 21000.0), 37.1216, 1e-3);
    CHECK_NEAR(farquhar_electron_transport(150.0, 30.0, 4.0), 25.6579, 1e-3);
    CHECK_NEAR(farquhar_assimilation(37.1216, 25.6579, 1.0), 24.6579, 1e-3);
    CHECK_NEAR(ball_berry_conductance(8.0, 20.0, 101325.0, 40.0, 4.0, 1.5, 1.5), 225166.6667,
               1e-3);
    CHECK_NEAR(medlyn_conductance(4.0, 20.0, 40.0, 1.0), 2.5, 1e-12);
    CHECK_NEAR(leuning_conductance(0.01, 10.0, 20.0, 30.0, 4.0, 0.8, 1.0), 6.163846, 1e-6);

    // Carbon balance, vegetation structure and soil water.
    CHECK_NEAR(npp_lpj(100.0, 20.0), 60.0, 1e-12);
    CHECK_NEAR(maintenance_respiration(1.0, 2.0, 30.0, 1.0, 4.0, 6.0, 0.0, 0.0, 0.0), 20.0,
               1e-12);
    CHECK_NEAR(lai_from_leaf_carbon(2.0, 2.5), 5.0, 1e-12);
    CHECK_NEAR(xylem_conductance(-1.0, -1.0), 0.5, 1e-12);
    CHECK_NEAR(soil_water_step(10.0, 5.0, 2.0, 1.0, 1.0, 1.0), 11.0, 1e-12);
    CHECK_NEAR(soil_thermal_inertia(100.0, 0.3, 0.1, 500.0), 200.0, 1e-12);
    CHECK_NEAR(vcmax_temperature(100.0, 25.0), 100.0, 1e-9);
    CHECK(vcmax_temperature(100.0, 35.0) > 100.0);

    // Spatial statistics.
    const std::vector<double> values = {1.0, 2.0, 3.0};
    const Weights identity = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    CHECK_NEAR(getis_ord_gi(values, identity, 1), 0.0, 1e-12);
    CHECK_NEAR(local_moran_i(values, identity, 1), 0.0, 1e-12);
    CHECK_NEAR(geodetector_q({0.0, 0.0, 10.0, 10.0}, {0, 0, 1, 1}), 1.0, 1e-12);
    CHECK_NEAR(geodetector_q({1.0, 2.0, 3.0}, {0, 0, 0}), 0.0, 1e-12);
    CHECK_NEAR(walker_threshold(0.05, 10), 0.0051162, 1e-7);
    const std::vector<double> adjusted = benjamini_hochberg({0.01, 0.04, 0.03, 0.005});
    CHECK_NEAR(adjusted[0], 0.02, 1e-12);
    CHECK_NEAR(adjusted[1], 0.04, 1e-12);
    CHECK_NEAR(adjusted[2], 0.04, 1e-12);
    CHECK_NEAR(adjusted[3], 0.02, 1e-12);
    CHECK(benjamini_hochberg_rejections({0.01, 0.04, 0.03, 0.005}, 0.05).size() == 4);
    CHECK_NEAR(poisson_scan_log_likelihood(10, 5.0, 100, 50.0), 69.314718, 1e-5);
}

void test_atmos_ocean() {
    // CAPE/CIN from a parcel profile: +2 K up to 1 km, neutral at 2 km, cold above.
    const std::vector<ParcelLevel> profile = {{0.0, 302.0, 300.0},
                                              {1000.0, 302.0, 300.0},
                                              {2000.0, 300.0, 300.0},
                                              {3000.0, 298.0, 300.0}};
    const CapeResult cape = cape_cin(profile);
    CHECK_NEAR(cape.cape_j_per_kg, 98.0665, 1e-3);
    CHECK(cape.found_lfc);
    CHECK(cape.found_lnb);
    CHECK_NEAR(cape.lfc_m, 0.0, 1e-12);
    CHECK_NEAR(cape.lnb_m, 2000.0, 1e-12);
    CHECK_NEAR(cape.cin_j_per_kg, 0.0, 1e-12);

    // Photolysis integral: constant phi sigma F over 0.3-0.4 um.
    CHECK_NEAR(photolysis_rate({0.3e-6, 0.4e-6}, {1.0, 1.0}, {1e-20, 1e-20}, {1e15, 1e15}),
               1e-12, 1e-24);
    CHECK_NEAR(falloff_rate_constant(300.0, 0.0), 0.0, 1e-30);
    CHECK_NEAR(falloff_rate_constant(300.0, 1e19), 6.7693178899e-13, 1e-23);
    CHECK_NEAR(aqueous_sulfate_oxidation_rate(1.0, 2.0, 4.0, 2.0, 3.0, 5.0, 0.5), 18.0, 1e-12);

    // Aerosol coagulation and ion balance.
    CHECK_NEAR(free_molecular_coagulation(1e-7, 1e-7, 298.0, 1000.0), 8.887925e-15, 1e-18);
    CHECK_NEAR(coagulation_sink({1.0, 2.0}, {3.0, 4.0}), 11.0, 1e-12);
    CHECK_NEAR(cloud_radiative_susceptibility(0.5, 2.0), 1.0, 1e-12);
    CHECK_NEAR(aci_radiative_forcing(0.5, 100.0, 0.01, 2.0, 2.0), 4.0, 1e-12);
    CHECK_NEAR(ion_equilibrium_number(1.0, 1.0, 0.0, 0.0), 1.0, 1e-12);
    CHECK_NEAR(ion_equilibrium_number(1.0, 1.0, 2.0, 1.0), 0.4142135624, 1e-9);
    CHECK_NEAR(recombination_limited_ions(1.0, 1.0), 1.0, 1e-12);

    // Ocean stratification, mixed layer, eddies and mixing.
    CHECK_NEAR(buoyancy_frequency_squared(9.80665, 1025.0, -1e-3), 9.56746e-6, 1e-10);
    CHECK_NEAR(mixed_layer_depth({0.0, 100.0}, {1e-4, 1e-4}), 50.0, 1e-9);
    CHECK_NEAR(ekman_velocity(0.1, 1e-4, 50.0), 20.0, 1e-12);
    CHECK_NEAR(barotropic_compensation(20.0, 50.0, 1000.0), -1.0, 1e-12);
    CHECK_NEAR(eddy_diffusivity_pv(1.0, 1e-11, 2e-11), 1e11, 1.0);
    CHECK_NEAR(eddy_diffusivity_gm(1e-7, 1e-3), 1e-4, 1e-15);
    CHECK_NEAR(internal_wave_frequency_squared(3.64605795e-5, 1e-4, 1e-3, 1e-2), 1.00133e-6,
               1e-12);
    CHECK_NEAR(gradient_richardson(9.80665, 1025.0, -1e-3, 1e-2, 0.0), 0.0956746341, 1e-9);
    CHECK_NEAR(osborn_diffusivity(0.2, 1e-8, 1e-4), 2e-5, 1e-15);
    CHECK_NEAR(dissipation_from_shear(1e-6, 1e-4), 7.5e-10, 1e-20);
    CHECK_NEAR(nori_viscosity(-0.1, 1e-3, 1e-1, 1e-5, 0.25, 0.05),
               0.0964387, 1e-6);
    CHECK_NEAR(nori_viscosity(0.1, 1e-3, 1e-1, 1e-5, 0.25, 0.05), 0.000604, 1e-9);
    CHECK_NEAR(nori_viscosity(0.3, 1e-3, 1e-1, 1e-5, 0.25, 0.05), 1e-5, 1e-12);

    // Carbon, heat fluxes and biogeochemistry.
    CHECK_NEAR(martin_export_flux(10.0, 1000.0, -0.858), 1.38675583, 1e-6);
    CHECK_NEAR(air_sea_co2_flux(5e-5, 0.03, 400.0, 300.0), 1.5e-4, 1e-15);
    CHECK_NEAR(dissolved_inorganic_carbon(2000.0, 100.0, 10.0), 2110.0, 1e-12);
    CHECK_NEAR(latent_heat_flux_bulk(1.2, 2.5e6, 1.2e-3, 10.0, 0.01, 0.0), 360.0, 1e-9);
    CHECK_NEAR(sensible_heat_flux_bulk(1.2, 1005.0, 1.2e-3, 10.0, 2.0, 0.0), 28.944, 1e-9);
    CHECK_NEAR(net_air_sea_heat_flux(100.0, -50.0, -20.0, -10.0), 20.0, 1e-12);
    CHECK_NEAR(iron_limitation(10.0, 1.0, 1.0, 1.0, 100.0, 100.0, 0.1, 0.5), 1.0 / 6.0, 1e-12);
    CHECK(free_iron_ligand_speciation(1e-9, 1e-10, 1e10) > 0.0);
    CHECK(free_iron_ligand_speciation(1e-9, 1e-10, 1e10) < 1e-10);
    CHECK_NEAR(mcdr_potential_ratio(2000.0, 2300.0, 0.5, 1.0), -0.4347826087, 1e-9);

    // Waves and sea-ice rheology.
    CHECK_NEAR(pierson_moskowitz_spectrum(1.0, 10.0), 0.3929092, 1e-6);
    CHECK_NEAR(doppler_shifted_frequency(1.0, 2.0, 3.0), 7.0, 1e-12);
    CHECK_NEAR(hibler_ice_pressure(2e4, 2.0, 0.9, 20.0), 5413.4113, 1e-3);
}

void test_closures() {
    // Updraft momentum budget and drag length.
    UpdraftParameters updraft;
    CHECK_NEAR(updraft_drag_length(updraft), 316.2278, 1e-3);
    CHECK_NEAR(updraft_momentum_tendency(1.0, 0.05, updraft), 0.0055387142, 1e-9);
    CHECK_NEAR(updraft_drag_sensitivity(1.0, updraft), -1.234568e-5, 1e-11);
    const UpdraftColumn column = integrate_updraft({0.0, 100.0, 200.0}, {0.01, 0.01, 0.01}, updraft,
                                                   0.0);
    CHECK(column.w_squared[1] > 0.0);
    CHECK(column.w_squared[2] > column.w_squared[1]);

    // Aerosol moment closure.
    CHECK_NEAR(mass_based_coagulation(1.0, 1.0, 1.0), 5.656854249, 1e-9);
    CHECK_NEAR(lognormal_moment_decay(1.0, 1.0, 1.0), 0.3928332311, 1e-9);

    // IDEMIX energetics.
    CHECK_NEAR(idemix_mu0(2.0, 1e-3, 1e-4, 5.24e-3, 7.29e-5), 0.4025694455, 1e-9);
    CHECK_NEAR(idemix_dissipation(0.4025694455, 1e-4, 1e-3, 0.01, 1.0), 0.0120498006,
               1e-9);
    CHECK_NEAR(idemix_dissipation_derivative(0.4025694455, 1e-4, 1e-3, 0.01, 1.0),
               0.0240996012, 1e-9);
    CHECK_NEAR(idemix_vertical_term(0.01, 15.0 * 86400.0, 1.0), -129.6, 1e-9);

    // KPP bulk Richardson number and its Vt sensitivity.
    CHECK_NEAR(kpp_bulk_richardson(0.0, -0.01, 100.0, 0.0, 0.0, 0.01), 10000.0, 1e-6);
    CHECK_NEAR(kpp_bulk_richardson_sensitivity(0.0, -0.01, 100.0, 0.0, 0.0, 0.01), -2e6, 1e-1);

    // Stochastic martingale pressure acceleration.
    const HorizontalVector acceleration = stochastic_pressure_acceleration(1.0, 0.0, 1000.0);
    CHECK_NEAR(acceleration.x, -0.001, 1e-12);
    CHECK_NEAR(acceleration.y, 0.0, 1e-12);

    // BBM rheology.
    CHECK_NEAR(bbm_effective_viscosity(1e5, 0.2, 2.0, 20.0, 0.9), 1172.2009, 1e-3);
    CHECK_NEAR(bbm_viscosity_damage_derivative(1e5, 0.2, 2.0, 20.0, 0.9), -2930.5022,
               1e-3);
    CHECK_NEAR(mohr_coulomb_yield(-1e5, 1e5, 30.0), 42264.9731, 1e-4);
    CHECK_NEAR(bbm_elastic_stress(0.5, 1e9, 1e-3), 5e5, 1e-6);
    CHECK_NEAR(bbm_friction_stress(1e4, 1e6, 2.0, 3.0), 3.0, 1e-12);
    CHECK_NEAR(bbm_friction_stress(-1e4, 1e6, 2.0, 3.0), 2.0, 1e-12);
    CHECK_NEAR(bbm_friction_stress(-1e7, 1e6, 2.0, 3.0), 3.0, 1e-12);

    // TEOS-10 generic Gibbs evaluation and sound speed from a known quadratic.
    std::vector<double> water(56, 0.0);
    water[1 * 7 + 0] = 1.0;
    CHECK_NEAR(teos10_gibbs_water(40.0, 0.0, water), 1.0, 1e-12);
    CHECK_NEAR(teos10_gibbs_water(80.0, 0.0, water), 2.0, 1e-12);
    CHECK_NEAR(teos10_gibbs_saline(40.0, 0.0, 0.0, {{2, 0, 0, 1.0}}), 1.0, 1e-12);

    const auto gibbs = [](double, double t, double p) {
        return 0.5 * t * t + 2.0 * t * p - 1.0 * p * p;
    };
    CHECK_NEAR(gibbs_sound_speed(gibbs, 35.0, 100.0, 1.0), 80.83316151, 1e-3);

    // Pure-water table against the official GSW reference values.
    const GibbsDerivatives pure = teos10_pure_water_derivatives(10.0, 1e7);
    CHECK_NEAR(pure.g, 9320.754144273566, 1e-6);
    CHECK_NEAR(pure.g_t, -150.0810634484812, 1e-7);
    CHECK_NEAR(pure.g_p, 0.0009955890810800538, 1e-16);
    CHECK_NEAR(sound_speed_from_derivatives(pure), 1463.531085477897, 1e-6);
    const GibbsDerivatives cold = teos10_pure_water_derivatives(0.0, 0.0);
    CHECK_NEAR(cold.g, 101.342743139674, 1e-10);
    CHECK_NEAR(cold.g_t, 0.1476445869773505, 1e-12);
    CHECK_NEAR(cold.g_p, 0.00100015695367145, 1e-16);
}

}  // namespace

void test_pde() {
    BlatterPattynParameters ice;
    ice.nx = 61;
    ice.nz = 31;
    ice.length = 50000.0;
    ice.thickness = 1000.0;
    ice.surface_slope = 0.01;
    ice.rate_factor = 1e-16;
    const BlatterPattynResult flow = solve_blatter_pattyn(ice);
    const double driving = 917.0 * 9.81 * 0.01;
    const double analytic_surface =
        2.0 * 1e-16 / 4.0 * std::pow(driving, 3.0) * std::pow(1000.0, 4.0);
    const std::size_t mid = (ice.nx / 2) * ice.nz + (ice.nz - 1);
    CHECK_NEAR(flow.velocity[mid], analytic_surface, analytic_surface * 0.12);
    CHECK(flow.iterations <= ice.max_iterations);

    StefanParameters conduction;
    conduction.cells = 101;
    conduction.length = 10.0;
    conduction.initial_temperature = 260.0;
    conduction.surface_temperature = 280.0;
    conduction.basal_temperature = 260.0;
    conduction.freezing_temperature = 200.0;
    conduction.time_step = 3600.0;
    conduction.steps = 28;
    const StefanResult linear = solve_stefan_1d(conduction);
    const double alpha = conduction.thermal_conductivity /
                         (conduction.density * conduction.specific_heat);
    const double time = conduction.time_step * static_cast<double>(conduction.steps);
    const double similarity = std::erfc(0.5 / (2.0 * std::sqrt(alpha * time)));
    CHECK_NEAR(linear.temperature[5], 260.0 + 20.0 * similarity, 1.0);
    CHECK_NEAR(linear.liquid_fraction[5], 1.0, 1e-9);

    StefanParameters freezing;
    freezing.initial_temperature = 268.0;
    freezing.surface_temperature = 278.0;
    freezing.basal_temperature = 268.0;
    freezing.steps = 48;
    const StefanResult thaw = solve_stefan_1d(freezing);
    CHECK(thaw.liquid_fraction[1] > 0.0);
    CHECK(thaw.liquid_fraction[1] < 1.0);
    CHECK(thaw.temperature[0] <= freezing.surface_temperature + 1e-9);

    PhaseFieldParameters fracture;
    fracture.cells = 101;
    fracture.length = 1.0;
    fracture.initial_phase.assign(fracture.cells, 1.0);
    fracture.tensile_energy = 0.0;
    fracture.fracture_toughness = 0.0;
    fracture.steps = 100;
    const std::vector<double> intact = solve_phase_field_1d(fracture);
    CHECK_NEAR(intact[50], 1.0, 1e-12);

    PhaseFieldParameters driven = fracture;
    driven.initial_phase.assign(driven.cells, 0.5);
    driven.tensile_energy = 1e5;
    driven.fracture_toughness = 0.0;
    driven.mobility = 1.0;
    driven.time_step = 1e-4;
    driven.steps = 1;
    const std::vector<double> healed = solve_phase_field_1d(driven);
    CHECK_NEAR(healed[50], 1.0, 1e-12);

    std::vector<double> density(20, 350.0);
    for (std::size_t i = 0; i < density.size(); ++i) density[i] += 5.0 * static_cast<double>(i);
    const std::vector<double> firn = evolve_firn_column(density, 1.0, 1e7, 10, 1e-3, 250.0);
    for (std::size_t i = 0; i < firn.size(); ++i) {
        CHECK(firn[i] >= density[i]);
        CHECK(firn[i] <= 917.0);
    }
}

void test_cryo() {
    CHECK_NEAR(glen_rate_factor(1e-16, 60000.0, 270.0), 2.469044184578022e-28, 1e-40);
    CHECK_NEAR(glen_strain_rate(1e-16, 1e5, 3.0), 0.1, 1e-12);
    CHECK_NEAR(effective_strain(0.1, 0.0, 0.0, 0.0, 0.0, 0.0), 0.0707106781, 1e-9);
    CHECK_NEAR(effective_stress(1e5, 0.0, 0.0, 0.0, 0.0, 0.0), 70710.6781, 1e-3);
    CHECK_NEAR(ice_viscosity(1e8, 0.1, 3.0), 5e9, 1e-3);
    CHECK_NEAR(sia_velocity(1e-16, 0.01, 1000.0), 1e-6, 1e-18);
    CHECK_NEAR(weertman_sliding(1e-3, 1e5, 1e4, 3.0, 1.0), 1e8, 1e-6);
    CHECK_NEAR(coulomb_sliding(0.5, 1e4), 5000.0, 1e-9);

    CHECK_NEAR(bulk_ice_stress(1.3, 1.5e-3, 10.0, 0.0, 0.0, 0.0, 0.0).x, 0.195, 1e-12);
    const HorizontalStress air{0.1, 0.0};
    const HorizontalStress ocean{0.0, 0.2};
    double du = 0.0, dv = 0.0;
    sea_ice_momentum_tendency(du, dv, 0.1, 0.0, 100.0, 1e-4, air, ocean, 0.0, 0.0,
                              HorizontalStress{});
    CHECK_NEAR(du, 0.001, 1e-12);
    CHECK_NEAR(dv, 0.00199, 1e-12);

    CHECK_NEAR(surface_melt_energy(0.5, 300.0, 200.0, 1.0, 273.15), 34.34217769919542, 1e-9);
    CHECK_NEAR(stefan_front_velocity(2.0, 5.0, 1.0, 1.0, 1000.0, 3.34e5), 2.6946107784e-8,
               1e-18);
    CHECK_NEAR(stefan_number(2000.0, 1500.0, 275.0, 273.0, 1.0, 1.0, 1000.0, 0.3, 3.34e5),
               0.0598802395, 1e-10);
    CHECK_NEAR(active_layer_depth(1.0, 1.5, 1000.0, 3.34e5), 0.0947736035, 1e-9);

    CHECK_NEAR(basal_melt_three_equation(1027.0, 917.0, 2.5e-3, 2.5e-5, 3.34e5, 2009.0,
                                         273.15, 271.25, 3985.0, 0.05, 272.65),
               1.1827226367e-9, 1e-16);
    CHECK_NEAR(basal_melt_simplified(3985.0, 3.34e5, 1027.0, 917.0, 1e-5, 1.0, -0.5, -1.9),
               1.8707295333e-7, 1e-15);
    CHECK(plume_melt(1.0, -2.0, 1.0, 7.86e-4, 34.0, 9.81, 0.05, 0.1, 1000.0, 0.5, 1.0, 0.0025,
                     0.1, 0.1, 1.0) > 0.0);

    CHECK_NEAR(subglacial_flux(1e-3, 0.1, 1.0, 2.0, -0.01), 1e-8, 1e-18);
    CHECK_NEAR(hydraulic_potential(1000.0, 9.81, 100.0, 5e5), 1481000.0, 1e-6);
    CHECK_NEAR(subglacial_melt(0.05, 1e5, 1e-4, 1e-4, 0.01, 4218.0, 1000.0, 0.0, 917.0,
                               3.34e5),
               3.2813329720e-8, 1e-16);

    CHECK_NEAR(calving_flux(1000.0, 1e4), 1e7, 1e-6);
    CHECK_NEAR(phase_field_notch(500.0, 1e5), 33000.00309374985, 1e-9);
    CHECK_NEAR(firn_densification_rate(1e5, 1e-3, 917.0, 400.0, 60000.0, 250.0),
               1.3842340922940609e-7, 1e-16);
    CHECK_NEAR(overburden_pressure(400.0, 10.0), 39226.6, 1e-9);
    CHECK_NEAR(herron_langway_gradient(1.0, 400.0, 917.0), 225.51799345692476, 1e-9);
    CHECK_NEAR(vertical_ice_velocity(0.1, 0.0, 0.0, 500.0, 1000.0, 1.0), -0.025, 1e-12);

    CHECK_NEAR(tipping_equilibrium(1.0, 2.0, 1.0, 1.0, 1.0), -1.0, 1e-12);
    CHECK_NEAR(tipping_critical_thickness(2.0, -1.0), 1.0, 1e-12);
    CHECK_NEAR(scalar_kalman_update(0.0, 10.0, 1.0, 1.0), 5.0, 1e-12);

    CHECK_NEAR(damaged_viscosity(1e9, 0.2), 8e8, 1e-3);
    CHECK_NEAR(damaged_glen_strain(1e-16, 0.5, 1e5, 3.0), 0.8, 1e-12);
    CHECK_NEAR(biot_effective_stress(1e5, 0.5, 4e4), 30000.0, 1e-9);
    CHECK_NEAR(maxwell_elasto_brittle_stress(0.2, 1e9, 1e-3), 8e5, 1e-6);
    const HorizontalStress floe = floe_acceleration({3.0, 4.0}, {1.0, 1.0}, {1.0, 1.0}, 10.0);
    CHECK_NEAR(floe.x, 0.5, 1e-12);
    CHECK_NEAR(floe.y, 0.6, 1e-12);
    CHECK_NEAR(roethlisberger_channel_rate(1e-4, 2.0, 1e5, 3.0, 1e-20, 1000.0, 917.0),
               8.905125409e-5, 1e-15);
    CHECK_NEAR(channel_energy_melt(1.0, -100.0, 1.0, 0.0, 0.0, 1000.0, 4218.0, 3.34e5),
               2.994011976e-7, 1e-15);
    CHECK_NEAR(phase_field_energy_density(0.5, 0.0, 100.0, 10.0, 1000.0, 1.0), 160.0, 1e-9);
    CHECK_NEAR(phase_field_evolution(0.5, 0.0, 100.0, 1000.0, 1.0, 1.0), -400.0, 1e-9);
    CHECK_NEAR(plume_flux_derivative(1027.0, 1e-6), 1.027e-3, 1e-15);
    CHECK_NEAR(plume_momentum_derivative(2.0, 9.81, 10.0, 0.1, 1.0), 18.62, 1e-9);
    CHECK_NEAR(probabilistic_age_increment(1.2, 0.1, 0.5, 10.0), 240.0, 1e-9);

    CHECK_NEAR(hayhurst_stress(0.21, 0.63, 1e5, 8e4, -3e4), 66600.0, 1e-9);
    CHECK_NEAR(hayhurst_damage_rate(66600.0, 50000.0, 2.0), 1.774224, 1e-9);
    CHECK_NEAR(hayhurst_damage_rate_derivative(66600.0, 50000.0, 2.0), 5.328e-5, 1e-12);
    CHECK_NEAR(hayhurst_damage_rate(40000.0, 50000.0, 2.0), 0.0, 1e-12);

    const ContactForce contact =
        hertz_mindlin_contact(1e9, 0.3, 0.1, 0.5, 0.01, 0.02, 1.0, 1.0, 1e-4, 0.0, 1e-5, 0.0,
                              0.5);
    CHECK_NEAR(contact.normal_stiffness, 7853981.633974481, 1e-3);
    CHECK_NEAR(contact.normal, 785.3981633974482, 1e-9);
    CHECK_NEAR(contact.shear_stiffness, 3020762.166913262, 1e-3);
    CHECK_NEAR(contact.tangential, 30.207621669132624, 1e-9);
    const ContactForce damped =
        hertz_mindlin_contact(1e9, 0.3, 0.1, 0.5, 0.01, 0.02, 1.0, 1.0, 1e-4, 0.01, 1e-5,
                              0.0, 0.5);
    CHECK_NEAR(damped.normal, 789.3614708222966, 1e-3);

    const std::vector<double> unit_dissolution = {1e-6, 1e-6, 1e-6};
    const std::vector<double> unit_age = {2.0, 2.0, 2.0};
    const std::vector<double> unit_mineral = {1.0, 1.0, 1.0};
    CHECK_NEAR(silicate_weathering_flux(0.5, unit_dissolution, unit_age, unit_mineral, 1.0,
                                       -0.5, 1.0),
               1.0606601717798212e-6, 1e-15);

    CHECK_NEAR(firn_strain_softening(0.01, 0.0, -0.02, 0.0, 0.0, 0.0, 0.02),
               1.118033988749895, 1e-9);

    const std::vector<std::vector<double>> covariance = {{4.0, 1.0}, {1.0, 9.0}};
    const std::vector<double> residual = {2.0, 3.0};
    CHECK_NEAR(icechrono_cost(covariance, residual), 1.714285714285714, 1e-9);
    const std::vector<std::vector<double>> jacobian = {{1.0, 0.0}, {0.0, 1.0}};
    const std::vector<double> gradient = icechrono_gradient(covariance, residual, jacobian);
    CHECK_NEAR(gradient[0], 0.857142857142857, 1e-9);
    CHECK_NEAR(gradient[1], 0.571428571428571, 1e-9);

    CHECK_NEAR(leaf_turnover(0.02, 1.0), 7.805604701973806e-6, 1e-15);
    CHECK_NEAR(root_uptake(0.3, 10.0, 0.2, 5.0), 4.0, 1e-12);
    CHECK_NEAR(transpiration_moisture(0.3, 0.4, 2.0, 0.8, 5.0), 2.25, 1e-12);
    CHECK_NEAR(optimal_stomatal_objective(10.0, 2.0, 10.0), 8.0, 1e-12);
    CHECK_NEAR(hydraulic_profit(10.0, 20.0, 0.1, 0.5, 2.0), 6.0, 1e-12);
    CHECK_NEAR(growth_temperature_acclimation(20.0, 0.5, 10.0), 25.0, 1e-12);
    CHECK_NEAR(mesophyll_limitation(10.0, 12.0, 2.0), 20.0, 1e-9);
    CHECK_NEAR(nutrient_colimitation(0.6, 0.4), 0.4, 1e-12);
    CHECK_NEAR(nutrient_demand(100.0, 20.0, 2.0), 3.0, 1e-12);
    CHECK_NEAR(nitrogen_fixation(1e-3, 0.5, 10.0, 3600.0), 18.0, 1e-9);
    CHECK_NEAR(canopy_air_exchange(1.2, 0.01, 0.005, 100.0), -6e-5, 1e-12);
    CHECK_NEAR(multilayer_scalar_tendency(10.0, 4.0, 2.0, 2.0), -1.0, 1e-12);
    CHECK_NEAR(soil_production(1e-3, 1.0, 1.0), 0.00036787944117144236, 1e-15);
    CHECK_NEAR(soil_production_reference(1e-6, 0.5, 50000.0, 280.0), 2.723667150940544e-7,
               1e-15);
}

void test_new_gaps() {
    CHECK_NEAR(green_ampt_infiltration(10.0, 100.0, 0.3, 0.1, 50.0), 14.0, 1e-9);
    CHECK_NEAR(green_ampt_infiltration_sensitivity(10.0, 100.0, 0.3, 0.1, 50.0), -0.08,
               1e-12);
    CHECK_NEAR(green_ampt_ponding_time(10.0, 100.0, 0.3, 0.1, 20.0), 1.0, 1e-12);

    const double geometry = perirhizal_geometry(2.0);
    CHECK_NEAR(geometry, 17.51562645708963, 1e-9);
    CHECK_NEAR(perirhizal_effective_conductivity(2.0, 1.0, 0.0, -1.0), 1.0, 1e-12);
    CHECK_NEAR(perirhizal_conductance(1.0, geometry, 0.01), 1.100539268012316, 1e-9);
    CHECK_NEAR(root_system_conductance(2.0, 1.0, 0.5), 1.0, 1e-12);
    CHECK_NEAR(root_system_conductance_sensitivity(2.0, 1.0, 0.5), 0.5, 1e-12);

    CHECK_NEAR(canopy_layer_flux(400.0, 380.0, 100.0), 0.2, 1e-12);
    CHECK_NEAR(canopy_layer_assimilation(2.0, 1.0, 400.0, 200.0, 100.0, 10000.0, 0.01),
               -2.0099502487562186, 1e-9);
    CHECK_NEAR(canopy_layer_assimilation_sensitivity(2.0, 1.0, 400.0, 200.0, 100.0, 10000.0,
                                                     0.01),
               1.0099750006187963e-6, 1e-15);

    CHECK_NEAR(mineral_dissolution_constant(1e-6, 2.0, 0.5, 50000.0, 280.0),
               3.4433720030310447e-7, 1e-15);
    CHECK_NEAR(soil_pco2_maximum(2.0), 1.5266552331761296, 1e-12);
    CHECK_NEAR(soil_pco2(400.0, 1.5266552331761296, 288.15), 400.92313648907077, 1e-9);
    CHECK_NEAR(carbonate_weathering_flux(1e-6, 0.4, 2.0, 1e-3), 8e-10, 1e-18);
    CHECK_NEAR(stream_power_erosion(1e-4, 4.0, 0.1), 2e-5, 1e-12);
    CHECK_NEAR(petrogenic_carbon_flux(0.5, 0.02, 2e-5), 1e-7, 1e-15);
    CHECK_NEAR(sulfide_weathering_flux(0.5, 0.01, 2e-5), 1e-7, 1e-15);
    CHECK_NEAR(microbial_uptake(1.0, 10.0, 5.0, 5.0, 2.0, 2.0), 0.9523809523809523, 1e-12);
    CHECK_NEAR(maoc_saturation(100.0, 10.0, 1.0, 1.0), 66.89085029457019, 1e-9);
    CHECK_NEAR(soil_microbial_stabilization(0.5, 0.2), 0.4806, 1e-12);
    CHECK_NEAR(arc_volcanic_flux(100.0, 0.4), 60.0, 1e-12);
    CHECK_NEAR(mantle_degassing_flux(0.1, 0.5, 100.0), 5.0, 1e-12);
    CHECK_NEAR(organic_phosphorus(1000.0, 50.0), 20.0, 1e-12);
    CHECK_NEAR(plant_phosphorus_uptake(20.0, 10.0, 5.0), 13.333333333333334, 1e-9);
    CHECK_NEAR(epic_stable_phosphorus_fraction(0.5), 0.0955, 1e-12);
    CHECK_NEAR(epic_labile_phosphorus_loss(0.5), 0.3775, 1e-12);
    CHECK_NEAR(weathering_front_thickness(1e-9, 1e10), 3.1622776601683795, 1e-9);

    const std::vector<double> enkf = {2.0, 1.0};
    const std::vector<double> enviar = {4.0, 3.0};
    const std::vector<double> hybrid = hybrid_analysis(enkf, enviar, 0.25);
    CHECK_NEAR(hybrid[0], 3.5, 1e-12);
    CHECK_NEAR(hybrid[1], 2.5, 1e-12);
    const std::vector<double> hybrid_sensitivity = hybrid_analysis_sensitivity(enkf, enviar);
    CHECK_NEAR(hybrid_sensitivity[0], -2.0, 1e-12);

    const std::vector<std::vector<double>> members = {{0.0}, {2.0}};
    const std::vector<std::vector<double>> h = {{1.0}};
    const std::vector<std::vector<double>> r = {{1.0}};
    const EnsembleAnalysis eakf =
        ensemble_adjustment_filter(members, {1.0}, h, r);
    CHECK_NEAR(eakf.mean[0], 1.0, 1e-9);
    CHECK_NEAR(eakf.covariance[0][0], 2.0 / 3.0, 1e-9);
    CHECK_NEAR(eakf.members[0][0], 0.42264973081037416, 1e-9);
    CHECK_NEAR(eakf.members[1][0], 1.5773502691896257, 1e-9);

    FeedForwardNetwork linear;
    linear.weights = {{{1.0, 2.0}, {3.0, 4.0}}};
    linear.biases = {{0.0, 0.0}};
    linear.activation = {0};
    const std::vector<double> output = network_forward(linear, {1.0, 1.0});
    CHECK_NEAR(output[0], 3.0, 1e-12);
    CHECK_NEAR(output[1], 7.0, 1e-12);
    const std::vector<std::vector<double>> jac = network_jacobian(linear, {1.0, 1.0});
    CHECK_NEAR(jac[0][0], 1.0, 1e-12);
    CHECK_NEAR(jac[1][1], 4.0, 1e-12);

    FeedForwardNetwork relu;
    relu.weights = {{{1.0}}};
    relu.biases = {{-1.0}};
    relu.activation = {1};
    CHECK_NEAR(network_forward(relu, {0.5})[0], 0.0, 1e-12);
    CHECK_NEAR(network_jacobian(relu, {0.5})[0][0], 0.0, 1e-12);
    CHECK_NEAR(network_jacobian(relu, {2.0})[0][0], 1.0, 1e-12);

    WeatheringParameters weathering;
    weathering.cells = 11;
    weathering.thickness = 10.0;
    weathering.production_rate = 0.01;
    weathering.time_step = 10.0;
    weathering.steps = 10;
    const WeatheringResult inert = solve_weathering_column(weathering);
    CHECK_NEAR(inert.mineral_fraction[0], 1.0, 1e-12);
    CHECK_NEAR(inert.mineral_fraction[10], 1.0, 1e-12);
    CHECK(inert.particle_age[5] > 0.0);
    weathering.dissolution.assign(weathering.cells, 1e-3);
    const WeatheringResult weathered = solve_weathering_column(weathering);
    CHECK(weathered.mineral_fraction[5] < 1.0);
}

void test_aerosol_ozone() {
    const MieEfficiencies small = mie_efficiencies(1.0, 1.5, 0.01);
    CHECK_NEAR(small.q_ext, 0.242479335500, 1e-9);
    CHECK_NEAR(small.q_sca, 0.213638571599, 1e-9);
    CHECK_NEAR(small.q_abs, 0.028840763901, 1e-9);
    const MieEfficiencies medium = mie_efficiencies(5.0, 1.33, 0.001);
    CHECK_NEAR(medium.q_ext, 3.579847888555, 1e-8);
    CHECK_NEAR(medium.q_sca, 3.558767512512, 1e-8);
    const MieEfficiencies absorbing = mie_efficiencies(10.0, 1.5, 0.1);
    CHECK_NEAR(absorbing.q_ext, 2.459790528444, 1e-8);
    CHECK_NEAR(absorbing.q_sca, 1.235144209371, 1e-8);

    const MieEfficiencies cross = mie_cross_sections(1e-6, 1e-6, 1.5, 0.01);
    CHECK_NEAR(cross.q_ext / (3.14159265358979323846 * 1e-12),
               mie_efficiencies(2.0 * 3.14159265358979323846, 1.5, 0.01).q_ext, 1e-9);
    CHECK_NEAR(aerosol_optical_depth({1.0, 2.0}, {100.0, 50.0}), 200.0, 1e-12);
    CHECK_NEAR(single_scattering_albedo(8.0, 10.0), 0.8, 1e-12);
    CHECK_NEAR(angstrom_optical_depth(1.0, 0.5, 0.55, 1.5), 1.1536897329871667, 1e-12);
    CHECK_NEAR(angstrom_exponent(1.0, std::pow(0.5 / 0.55, -2.0), 0.55, 0.5), 2.0, 1e-12);
    CHECK_NEAR(direct_aerosol_forcing(1361.0, 0.8, 0.5, 0.1, 0.95, 0.7, 0.1), -5.2281454,
               1e-7);

    double d_up = 0.0, d_down = 0.0;
    two_stream_derivatives(10.0, 5.0, 0.2, 0.9, 100.0, 0.5, 1.0, 0.5, 0.7, 0.3, d_up, d_down);
    CHECK_NEAR(d_up, -34.730085818, 1e-3);
    CHECK_NEAR(d_down, -18.098846306, 1e-3);
    CHECK_NEAR(polarized_reflectance(1.0, 1.0, 1000.0, 0.5),
               3.14159265358979323846 / 500.0 * std::sqrt(2.0), 1e-12);

    CHECK_NEAR(kohler_activated_droplets(100.0, 1.0, 0.1), 9.090909090909092, 1e-12);
    CHECK_NEAR(aci_index(-0.1, 0.2), 0.5, 1e-12);
    CHECK_NEAR(aci_optical_depth_index(0.3, 0.2), 1.5, 1e-12);
    CHECK_NEAR(aci_forcing_decomposition(1.0, -2.0, -0.5, -0.2), -1.7, 1e-12);
    CHECK_NEAR(aci_residual_twomey(2.0, 0.5, 0.3), 1.2, 1e-12);
    CHECK_NEAR(ml_aci_forcing(300.0, 0.5, 0.6, 0.1, 0.2, 0.1), -0.18, 1e-12);

    CHECK_NEAR(chapman_ozone_production(1e-5, 1e12), 1e7, 1e-3);
    CHECK_NEAR(chapman_odd_oxygen_loss(2e-15, 1e6, 1e12), -2e3, 1e-9);
    CHECK_NEAR(nox_odd_oxygen_loss(1e-15, 1e9, 1e6), -2.0, 1e-12);
    CHECK_NEAR(chapman_ozone_ratio(4e-5, 4e-6), 3.1622776601683795, 1e-12);
    CHECK_NEAR(methane_fractional_release(10.0, 6.0), 0.4, 1e-12);
    CHECK_NEAR(age_of_air(0.4), 3.6536, 1e-9);
    CHECK_NEAR(stratospheric_water_vapour(3.0, 1.7, 0.7), 5.0, 1e-12);
    CHECK_NEAR(stratospheric_water_vapour_change(0.4, 0.2), 0.16, 1e-12);
    CHECK_NEAR(water_vapour_mass_perturbation({1.0, 2.0}, {100.0, 50.0}, 18.0, 29.0),
               200.0 * 18.0 / 29.0, 1e-9);
}

void test_resolved_gaps() {
    const double flux = schoof_grounding_line_flux(1e-16, 917.0, 1027.0, 9.81, 3.0, 3.0, 1000.0,
                                                   1.0);
    CHECK_NEAR(flux, 16471727584.013472, 1e-3);
    CHECK_NEAR(grounding_line_migration(flux, 1000.0), 16471727.584013471, 1e-6);
    CHECK_NEAR(schoof_flux_thickness_sensitivity(flux, 1000.0, 3.0, 3.0), 78240706.02406399,
               1e-5);

    CHECK_NEAR(kpp_shape_function(0.5, 1.0, 2.0, 3.0), 1.375, 1e-12);
    CHECK_NEAR(kpp_velocity_scale(2.0, 0.5, 1.0, 2.0, 3.0), 2.75, 1e-12);

    CHECK_NEAR(soil_pco2_maximum_sensitivity(2.0), 0.21066209327045182, 1e-12);
    CHECK_NEAR(stream_power_erosion_runoff_sensitivity(1e-4, 4.0, 0.1), 2.5e-6, 1e-15);
    CHECK_NEAR(stream_power_erosion_slope_sensitivity(1e-4, 4.0), 2e-4, 1e-15);
    CHECK_NEAR(soil_production_reference_sensitivity(1e-6, 50000.0, 280.0), 5.447334301881088e-7,
               1e-18);

    CHECK_NEAR(van_genuchten_effective_saturation(0.0, 1.4, 2.0), 1.0, 1e-12);
    CHECK_NEAR(van_genuchten_effective_saturation(1.0, 1.4, 2.0), 0.5812381937190964, 1e-12);
    CHECK_NEAR(van_genuchten_conductivity(1.0, 1e-5, 1.4, 2.0), 2.645128133968696e-7, 1e-18);
    CHECK_NEAR(van_genuchten_matric_flux_potential(1.0, 1e-5, 1.4, 2.0), 3.3944315543502036e-6,
               1e-18);
    CHECK_NEAR(perirhizal_conductance_van_genuchten(1.0, 17.51562645708963, 2.0, 1.0, 1e-5, 1.4,
                                                    2.0),
               0.00028683924701039445, 1e-15);
}

void test_advanced_ozone_landuse() {
    CHECK_NEAR(oxygen_damping_rate(0.1, 1e6, 2e6, 3e9, 4e8, 1e7, 1e-12, 1e-12, 1e-12, 1e-12,
                                   1e-12),
               0.10682300000000002, 1e-15);
    CHECK_NEAR(ozone_damping_rate(0.1, 1e6, 2e6, 3e6, 1e-12, 2e-12, 3e-12), 0.100014, 1e-15);
    CHECK_NEAR(odd_oxygen_chemical_tendency(0.1, 2.0, 0.2, 3.0), -0.8, 1e-12);
    CHECK_NEAR(effective_ox_damping(0.2, 0.2, 0.1, 0.1, 1.0), 0.11, 1e-12);
    CHECK_NEAR(water_vapour_feedback(1.0, 0.5), 2.0, 1e-12);
    CHECK_NEAR(tropical_adjustment(5.0, 2.0), 3.0, 1e-12);
    const std::vector<double> delta_t =
        temperature_adjustment({{2.0, 0.0}, {0.0, 4.0}}, {2.0, 4.0});
    CHECK_NEAR(delta_t[0], -1.0, 1e-12);
    CHECK_NEAR(delta_t[1], -1.0, 1e-12);
    CHECK_NEAR(methane_water_trend(0.5, 0.1), 0.1, 1e-12);
    CHECK_NEAR(cold_trap_entry_change(2.0), 0.9, 1e-12);
    CHECK_NEAR(ozone_kernel_forcing({1.0, 2.0}, {3.0, 4.0}), 11.0, 1e-12);
    CHECK_NEAR(stratospheric_adjusted_forcing(1.0, 2.0), 3.0, 1e-12);

    CHECK_NEAR(mixed_albedo({{0.5, 1.0, 1.0, 0.2}, {0.5, 1.0, 1.0, 0.6}}), 0.4, 1e-12);
    CHECK_NEAR(annual_albedo_change({0.3, 0.5}, {0.2, 0.4}), 0.1, 1e-12);
    CHECK_NEAR(albedo_radiative_forcing(10.0, 0.1), 1.0, 1e-12);
    CHECK_NEAR(clear_sky_albedo_forcing(0.1, 1000.0), -70.0, 1e-12);
    CHECK_NEAR(global_conversion_forcing({{0.1}}, {{0.5}}, {{1e6}}, {{10.0}}, 5.1e8),
               0.000980392156862745, 1e-15);
    CHECK_NEAR(canopy_resistance(100.0, 4.0), 50.0, 1e-12);
    CHECK_NEAR(canopy_penman_monteith(1.0, 100.0, 1.2, 1005.0, 1.0, 50.0, 0.066, 50.0),
               109.64664310954063, 1e-9);
    CHECK_NEAR(surface_temperature_decomposition(200.0, 0.0, 0.2, 0.0, 0.95, 290.0, 0.0, 2.5e6,
                                                 0.0, 1.2, 0.01, 300.0, 0.0, 50.0, 0.0),
               0.0, 1e-12);
    CHECK_NEAR(bookkeeping_carbon_flux(1.0, 100.0, 200.0, 0.5, 0.2, 0.1), 100.0, 1e-12);
    CHECK_NEAR(stock_difference(120.0, 100.0, 5.0), 4.0, 1e-12);
    CHECK_NEAR(immediate_emission_factor(10.0), 0.5, 1e-12);
    CHECK_NEAR(gain_loss_stock_change(5.0, 2.0, 0.5, 10.0, 8.0, 2.0), 8.0, 1e-12);

    double vegetation = 0.0, litter = 0.0, soil = 0.0, product = 0.0;
    land_use_pool_tendencies(vegetation, litter, soil, product, 10.0, 2.0, 3.0, 1.0, 1.0, 2.0,
                             0.5, 0.5);
    CHECK_NEAR(vegetation, 4.0, 1e-12);
    CHECK_NEAR(litter, 0.0, 1e-12);
    CHECK_NEAR(soil, 1.5, 1e-12);
    CHECK_NEAR(product, 0.5, 1e-12);
    CHECK_NEAR(net_land_use_flux(10.0, 5.0, 20.0, 2.0, 8.0, 4.0, 18.0), 7.0, 1e-12);

    const std::vector<double> ridge = ridge_regression({{1.0, 0.0}, {0.0, 1.0}}, {2.0, 3.0}, 1.0);
    CHECK_NEAR(ridge[0], 1.0, 1e-12);
    CHECK_NEAR(ridge[1], 1.5, 1e-12);
    CHECK_NEAR(linear_prediction(ridge, {2.0, 4.0}, 0.5), 8.5, 1e-12);
}

void test_change_detection() {
    CHECK_NEAR(gev_cdf(0.0, 0.0, 1.0, 0.0), 0.36787944117144233, 1e-12);
    CHECK_NEAR(gev_return_level(0.5, 0.0, 1.0, 0.0), 0.36651292058166435, 1e-12);
    CHECK_NEAR(gp_cdf(1.0, 1.0, 0.0), 0.6321205588285577, 1e-12);
    CHECK_NEAR(ou_stationary_variance(2.0, 2.0), 1.0, 1e-12);
    CHECK_NEAR(ou_lag1_autocorrelation(2.0, 0.5), 0.36787944117144233, 1e-12);
    CHECK_NEAR(ar1_estimate({1.0, 2.0, 3.0, 4.0}), 1.0, 1e-12);
    CHECK_NEAR(rate_mismatch_criterion(2.0, 0.2, 1.0, 0.1), 4.0, 1e-12);
    CHECK_NEAR(regime_shift_index({1.0, 2.0, 3.0}, 0.0, 1.0), 2.0, 1e-12);
    CHECK_NEAR(cusum_statistic({1.0, 2.0, 3.0, 4.0}), 2.0, 1e-12);
    CHECK_NEAR(optimal_fingerprint({{1.0}, {2.0}}, {2.0, 4.0}, {{1.0, 0.0}, {0.0, 1.0}})[0], 2.0,
               1e-12);
    CHECK_NEAR(marchenko_pastur_bounds_lower(1.0, 4.0), 1.0, 1e-12);
    CHECK_NEAR(marchenko_pastur_bounds_upper(1.0, 4.0), 9.0, 1e-12);
    CHECK_NEAR(marchenko_pastur_density(5.0, 1.0, 4.0), 0.12732395447351627, 1e-12);

    CHECK_NEAR(area_mean_anomaly({1.0, 2.0}, {1.0, 1.0}), 1.5, 1e-12);
    CHECK_NEAR(relative_nino34(1.5, 0.0, 0.5, 0.0), 1.0, 1e-12);
    CHECK_NEAR(amo_relative(1.0, 0.4), 0.6, 1e-12);
    CHECK_NEAR(pdo_projection({1.0, 2.0}, {3.0, 4.0}), 11.0, 1e-12);
    CHECK_NEAR(amoc_index({{1.0, 2.0}, {3.0, 1.0}}), 3.0, 1e-12);
    CHECK_NEAR(ocean_heat_content(1030.0, 3980.0, {1.0, 2.0}, {0.0, 0.0}, {100.0, 50.0}),
               819880000.0, 1e-6);
    CHECK_NEAR(layer_heat_content_anomaly(1030.0, 3980.0, 1.0, 100.0), 409940000.0, 1e-6);
    CHECK_NEAR(earth_energy_imbalance(340.0, 240.0, 100.0), 0.0, 1e-12);
    CHECK_NEAR(energy_imbalance_decomposition(1.0, 1.0, 0.5), 0.5, 1e-12);
    CHECK_NEAR(sea_level_budget(1.0, 2.0, 3.0, 4.0), 10.0, 1e-12);
    CHECK_NEAR(barystatic_sea_level(1.0, 2.0, 3.0, 4.0, 0.5), 10.5, 1e-12);
    CHECK_NEAR(quadratic_sea_level(0.0, 1.0, 0.5, 2.0), 4.0, 1e-12);
    CHECK_NEAR(quadratic_sea_level_acceleration(0.5), 1.0, 1e-12);
    CHECK_NEAR(broken_linear_sea_level(0.0, 1.0, 2.0, 3.0, 2.0), 5.0, 1e-12);
    CHECK_NEAR(groundwater_storage_anomaly(10.0, 1.0, 2.0, 3.0), 4.0, 1e-12);
    CHECK_NEAR(recharge_convolution(0.0, {1.0, 2.0}, 1.0), 2.36787944117144233, 1e-12);
    CHECK_NEAR(storage_depletion({1.0, 2.0}, {10.0, 10.0}, 0.1), 3.0, 1e-12);
    CHECK_NEAR(heston_variance_step(0.04, 1.0, 0.04, 0.1, 0.1, 0.0), 0.04, 1e-12);
    const std::vector<double> stationary =
        markov_stationary_distribution({{0.9, 0.1}, {0.2, 0.8}});
    CHECK_NEAR(stationary[0], 2.0 / 3.0, 1e-9);
    CHECK_NEAR(regime_expected_value({1.0, 2.0}, {0.25, 0.75}), 1.75, 1e-12);

    DynamicLinearState state;
    state.update(1.0, 1.0, 0.0, 0.0);
    CHECK_NEAR(state.level, 2.0 / 3.0, 1e-12);
    CHECK_NEAR(state.rate, 1.0 / 3.0, 1e-12);
}

void test_attribution_and_ml() {
    CHECK_NEAR(prp_feedback(2.0, 1.0, 3.0, 1.0), 0.5, 1e-12);
    CHECK_NEAR(prp_feedback_two_sided(2.0, 1.0, 3.0, 2.0, 1.0, 1.0), 1.0, 1e-12);
    CHECK_NEAR(fSST_erf({2.0, 4.0}, {1.0, 1.0}), 2.0, 1e-12);
    CHECK_NEAR(erf_decomposition(1.0, {0.5, 0.5}), 2.0, 1e-12);
    CHECK_NEAR(cloud_adjustment_factor(0.5, 0.3, 0.1, 0.2, 0.5, 0.2, 0.3, 1.0), 0.195,
               1e-12);
    CHECK_NEAR(surface_temperature_eight_term(200.0, 0.0, 0.2, 0.0, 0.95, 290.0, 0.0, 2.5e6,
                                              0.0, 1.2, 0.01, 300.0, 0.0, 50.0, 0.0, 0.0, 0.0),
               0.0, 1e-12);
    CHECK_NEAR(surface_temperature_eight_term(200.0, 0.0, 0.2, 0.0, 0.95, 290.0, 0.0, 2.5e6,
                                              0.0, 1.2, 0.01, 300.0, 0.0, 50.0, 0.0, -5.0, 0.0),
               5.0 / (1.2 * 0.01 / 50.0 + 4.0 * 0.95 * 5.670374419e-8 * 300.0 * 300.0 * 300.0),
               1e-12);
    CHECK_NEAR(stratospheric_swv_feedback_target().median, 0.25, 1e-15);
    CHECK_NEAR(full_column_swv_feedback_target().spread, 0.007, 1e-15);

    const Matrix predictors = {{0.0}, {1.0}, {2.0}, {3.0}};
    const std::vector<double> response = {0.0, 1.0, 2.0, 3.0};
    RegressionTree tree;
    tree.fit(predictors, response, 8, 1, 1, 7);
    CHECK_NEAR(tree.predict({1.0}), 1.0, 1e-12);
    CHECK_NEAR(tree.importance[0], 1.0, 1e-12);

    RandomForest forest;
    forest.fit(predictors, response, 25, 8, 2, 1, 11);
    CHECK_NEAR(forest.predict({1.5}), 1.5, 0.6);
    CHECK_NEAR(forest.deweathered_prediction({2.0}), 2.0, 0.6);

    CHECK_NEAR(causal_forest_split_criterion(2.0, 2.0, 1.0, 3.0), 1.0, 1e-12);
    CHECK_NEAR(cate_estimate({3.0, 5.0}, {1.0, 2.0}), 2.5, 1e-12);
}

void test_copula_hmm() {
    CHECK_NEAR(copula_parameter(CopulaFamily::Gumbel, 0.0, 1.0, 0.0, 0.0), 2.0, 1e-12);
    CHECK_NEAR(copula_parameter(CopulaFamily::Clayton, 0.0, 1.0, 0.0, 0.0), 0.0, 1e-12);
    CHECK_NEAR(copula_parameter(CopulaFamily::HuelserReiss, 0.0, 1.0, 0.0, 0.0), 1.0, 1e-12);
    CHECK_NEAR(copula_parameter_derivative(CopulaFamily::Gumbel, 0.0, 1.0, 0.0, 0.0), 1.0,
               1e-12);
    CHECK_NEAR(copula_parameter_derivative(CopulaFamily::Gumbel, 0.0, 0.0, 0.5, 1.0),
               std::exp(0.5), 1e-12);

    CHECK_NEAR(honest_emse_criterion(4.0, 2.0, 2.0, 2.0), 0.0, 1e-12);

    GaussianHmm single;
    single.initial = {1.0};
    single.transition = {{1.0}};
    single.means = {0.0};
    single.variances = {1.0};
    const HmmPosterior single_posterior = hmm_forward_backward(single, {0.0, 0.0});
    CHECK_NEAR(single_posterior.log_likelihood, 2.0 * std::log(1.0 / std::sqrt(2.0 * M_PI)),
               1e-12);
    CHECK_NEAR(single_posterior.gamma[1][0], 1.0, 1e-12);

    GaussianHmm two;
    two.initial = {0.5, 0.5};
    two.transition = {{0.9, 0.1}, {0.1, 0.9}};
    two.means = {0.0, 5.0};
    two.variances = {1.0, 1.0};
    const HmmPosterior posterior = hmm_forward_backward(two, {0.0, 0.0, 5.0, 5.0});
    CHECK(posterior.gamma[0][0] > 0.99);
    CHECK(posterior.gamma[2][1] > 0.99);
    CHECK(posterior.xi.size() == 3);
    CHECK_NEAR(posterior.gamma[1][0] + posterior.gamma[1][1], 1.0, 1e-12);

    GaussianHmm initial;
    initial.initial = {0.5, 0.5};
    initial.transition = {{0.9, 0.1}, {0.1, 0.9}};
    initial.means = {1.0, 4.0};
    initial.variances = {1.0, 1.0};
    const GaussianHmm trained = baum_welch(initial, {0.0, 0.0, 0.0, 5.0, 5.0, 5.0}, 50, 1e-9);
    CHECK(trained.means[0] < 2.0);
    CHECK(trained.means[1] > 3.0);
    CHECK_NEAR(trained.transition[0][0] + trained.transition[0][1], 1.0, 1e-9);

    const Matrix hovmoller = hovmoller_average({{{1.0, 2.0}, {3.0, 4.0}}});
    CHECK_NEAR(hovmoller[0][0], 2.0, 1e-12);
    CHECK_NEAR(hovmoller[0][1], 3.0, 1e-12);
    const std::vector<double> derivative = phase_space_derivative({0.0, 1.0, 4.0}, 1.0);
    CHECK_NEAR(derivative[0], 1.0, 1e-12);
    CHECK_NEAR(derivative[1], 2.0, 1e-12);
    CHECK_NEAR(derivative[2], 3.0, 1e-12);
}

void test_regimes_and_causal() {
    CHECK_NEAR(joe_copula(0.5, 0.5, 2.0), 0.3385621722338523, 1e-12);
    CHECK(std::isfinite(symmetric_stable_sample(1.5, 0.7, 1.0)));
    CHECK_NEAR(climate_feedback_manifold(2.0, 1.0, 2.0), 2.5, 1e-12);
    CHECK_NEAR(temperature_tendency(1.0, 2.0, 1.0, 2.0), -4.0, 1e-12);
    CHECK_NEAR(logistic_transition(0.0, 1.0, 0.0), 0.5, 1e-12);
    CHECK_NEAR(exponential_transition(0.0, 1.0, 0.0), 0.0, 1e-12);
    CHECK_NEAR(multiple_lstr(1.0, 2.0, 3.0, 0.5, 0.2), 2.6, 1e-12);

    CHECK_NEAR(stochastic_volatility_observation(0.0, 2.0), 2.0, 1e-12);
    CHECK_NEAR(stochastic_volatility_step(0.0, 0.0, 0.9, 0.1, 1.0), 0.1, 1e-12);
    CHECK_NEAR(stochastic_volatility_mean_reversion(0.9, 1.0), 0.10536051565782628, 1e-12);
    CHECK_NEAR(stochastic_volatility_intercept(0.01, 0.9), 0.1, 1e-12);

    CHECK_NEAR(gjr_garch_variance(0.01, 0.05, 0.1, 0.8, 0.04, -0.2), 0.048, 1e-12);
    CHECK_NEAR(beta_weight(1, 2, 1.0, 1.0), 1.0, 1e-12);
    CHECK_NEAR(long_run_log_variance(1.0, 2.0, {0.5, 0.5}, {1.0, 3.0}), 5.0, 1e-12);
    CHECK_NEAR(garch_midas_variance(0.5, 1.0, 2.0, 0.1, 0.8), 0.55, 1e-12);

    CHECK_NEAR(asymmetric_laplace_density(0.0, 0.0, 1.0, 0.5), 0.25, 1e-12);
    CHECK_NEAR(average_seasonal_trends({1.0, 2.0, 3.0}), 2.0, 1e-12);
    CHECK_NEAR(seasonal_trend_standard_error({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}),
               0.5773502691896258, 1e-12);

    CHECK_NEAR(beast_seasonal(0.0, 4.0, {1.0}, {0.0}), 1.0, 1e-12);
    CHECK_NEAR(beast_seasonal(1.0, 4.0, {1.0}, {0.0}), 0.0, 1e-12);
    CHECK_NEAR(beast_steps(5.0, {2.0, 4.0}, {1.0, 2.0}), 3.0, 1e-12);

    CHECK_NEAR(gev_log_likelihood(0.0, 0.0, 1.0, 0.0), -1.0, 1e-12);
    CHECK(significant_eigenvalue_count({10.0, 1.0, 0.5}, 1.0, 4.0) == 1);
    const Matrix reconstruction =
        rmt_reconstruction({2.0, 3.0}, {{1.0, 0.0}, {0.0, 1.0}});
    CHECK_NEAR(reconstruction[0][0], 2.0, 1e-12);
    CHECK_NEAR(reconstruction[1][1], 3.0, 1e-12);
    CHECK_NEAR(total_least_squares({{1.0}, {2.0}}, {2.0, 4.0})[0], 2.0, 1e-12);

    const Matrix design = sem_interaction_design({1.0, 2.0}, {3.0, 4.0});
    CHECK_NEAR(design[1][2], 8.0, 1e-12);
    const std::vector<double> demeaned = two_way_demean({1.0, 3.0, 5.0, 7.0}, {0, 0, 1, 1}, {0, 1, 0, 1});
    for (double value : demeaned) CHECK_NEAR(value, 0.0, 1e-12);
}

void test_spatial_advanced() {
    const std::vector<double> values = {1.0, 2.0, 3.0};
    const Weights weights = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    CHECK(lisa_classification(values, weights, 0) == LisaClass::LowLow);
    CHECK(lisa_classification(values, weights, 2) == LisaClass::HighHigh);
    CHECK_NEAR(gaussian_kernel(0.0, 2.0), 1.0, 1e-12);
    CHECK_NEAR(gaussian_kernel(2.0, 2.0), std::exp(-0.5), 1e-12);

    const std::vector<double> gwr =
        gwr_coefficients({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}}, {2.0, 3.0, 5.0},
                         {0.0, 0.0, 0.0}, 1.0);
    CHECK_NEAR(gwr[0], 2.0, 1e-12);
    CHECK_NEAR(gwr[1], 3.0, 1e-12);

    CHECK_NEAR(matern32_covariance(0.0, 2.0, 1.0), 2.0, 1e-12);
    CHECK_NEAR(matern32_covariance(1.0, 2.0, 1.0), 2.0 * (1.0 + std::sqrt(3.0)) * std::exp(-std::sqrt(3.0)),
               1e-12);
    CHECK_NEAR(probit_category_probability(0.0, -1e300, 0.0, 1.0), 0.5, 1e-12);
    CHECK_NEAR(hotelling_functional_scan(1, 1, {1.0}, {0.0}, {{1.0}}), 0.5, 1e-12);
    CHECK_NEAR(crisp_criterion({1.0, 3.0}, {1, 0}, 0.0, 3.0), 0.0, 1e-12);
    CHECK_NEAR(crisp_criterion({1.0, 3.0}, {1, 0}, 2.0, 3.0), 6.0, 1e-12);
    CHECK_NEAR(log_gaussian_cox_intensity(1.0, {1.0}, {2.0}, 0.5), 3.5, 1e-12);
    CHECK_NEAR(deviance_information_criterion(10.0, 8.0), 12.0, 1e-12);

    CHECK_NEAR(quadratic_trend_level(1.0, 2.0, 0.5, 2.0), 7.0, 1e-12);
    CHECK_NEAR(quadratic_trend_rate(2.0, 0.5, 2.0), 4.0, 1e-12);
    CHECK_NEAR(quadratic_trend_acceleration(0.5), 1.0, 1e-12);
    CHECK_NEAR(linearithmic_trend_level(1.0, 2.0, std::exp(1.0)), 1.0 + 2.0 * std::exp(1.0), 1e-12);
    CHECK_NEAR(linearithmic_trend_rate(2.0, std::exp(1.0)), 4.0, 1e-12);
    CHECK_NEAR(linearithmic_trend_acceleration(2.0, std::exp(1.0)), 2.0 / std::exp(1.0), 1e-12);
    CHECK_NEAR(step_magnitude(5.0, 3.0), 2.0, 1e-12);
    CHECK_NEAR(integrated_change({0.0, 1.0, 3.0}, 1.0), 2.5, 1e-12);
    CHECK_NEAR(return_level_shift(5.0, 2.0), 3.0, 1e-12);
    CHECK_NEAR(newey_west_variance({1.0, 1.0, 1.0}, 1), 0.0, 1e-12);

    std::vector<double> series(16);
    for (std::size_t i = 0; i < series.size(); ++i)
        series[i] = std::sin(0.5 * static_cast<double>(i));
    const ComplexField field = morlet_cwt(series, 1.0, {1.0, 2.0, 4.0});
    CHECK(field.size() == 3);
    CHECK(field[0].size() == series.size());
    CHECK_NEAR(wavelet_coherence(field, field, 8, 1, 3, 1), 1.0, 1e-12);
}

void test_magnitude_hierarchy() {
    CHECK_NEAR(piecewise_linear_value(0.0, 1.0, 2.0, 3.0, 2.0), 5.0, 1e-12);
    CHECK_NEAR(piecewise_rate_before(1.0), 1.0, 1e-12);
    CHECK_NEAR(piecewise_rate_after(1.0, 2.0), 3.0, 1e-12);
    CHECK_NEAR(effective_return_level(0.5, 0.0, 1.0, 0.0), 0.36651292058166435, 1e-12);
    CHECK_NEAR(return_period_derivative(1.0, 0.5, 0.1), 0.4, 1e-12);

    CHECK_NEAR(ar2_step(1.0, 2.0, 0.5, 0.3, 0.0), 1.1, 1e-12);
    CHECK_NEAR(ar2_stability(0.5, 0.3), 0.8, 1e-12);
    CHECK_NEAR(oht_amplitude_term(2.0, 3.0), 6.0, 1e-12);
    CHECK_NEAR(oht_circulation_term(4.0, 5.0), 20.0, 1e-12);
    CHECK_NEAR(oht_covariance_term(2.0, 5.0), 10.0, 1e-12);

    CHECK_NEAR(fisher_information_gaussian_mean(4.0), 0.25, 1e-12);
    CHECK_NEAR(information_speed(0.25, 2.0), 1.0, 1e-12);
    CHECK_NEAR(clipped_standardised_value(5.0, 0.0, 1.0, 2.0), 2.0, 1e-12);
    CHECK_NEAR(clipped_standardised_value(-5.0, 0.0, 1.0, 2.0), -2.0, 1e-12);
    CHECK_NEAR(local_contrast(0.0, 2.0, 1.0, 1.0, 0.0, 4.0), 2.8284271247461903, 1e-12);
    CHECK_NEAR(max_cusum_statistic({-1.0, 3.0, -2.0}), 3.0, 1e-12);
    CHECK(higher_criticism({0.001, 0.5}) > 20.0);

    CHECK_NEAR(gaussian_covariance_log_likelihood_ratio(10, 10, {{1.0}}, {{1.0}}, {{1.0}}), 0.0,
               1e-12);
    CHECK_NEAR(optimal_transport_cost_1d({0.0, 2.0}, {1.0, 3.0}), 1.0, 1e-12);

    const std::vector<double> slope =
        functional_slope({{0.0, 0.0}, {1.0, 2.0}, {2.0, 4.0}}, {0.0, 1.0, 2.0});
    CHECK_NEAR(slope[0], 1.0, 1e-12);
    CHECK_NEAR(slope[1], 2.0, 1e-12);
    CHECK_NEAR(functional_slope_norm(slope, 1.0), 5.0, 1e-12);

    CHECK_NEAR(combined_standard_uncertainty({2.0, 3.0}, {1.0, 1.0}), std::sqrt(13.0), 1e-12);
    CHECK_NEAR(correlated_standard_uncertainty({1.0, 1.0}, {1.0, 1.0}, {{1.0, 0.0}, {0.0, 1.0}}),
               std::sqrt(2.0), 1e-12);
    CHECK_NEAR(correlated_mean_variance({1.0, 1.0}, {{1.0, 0.0}, {0.0, 1.0}}), 0.5, 1e-12);
    CHECK_NEAR(two_component_error(3.0, 4.0), 5.0, 1e-12);
    CHECK_NEAR(reduced_chi_square({1.0, 3.0}, {0.0, 0.0}, {1.0, 1.0}), 10.0, 1e-12);
    CHECK_NEAR(root_mean_square_error({1.0, 3.0}, {0.0, 0.0}), std::sqrt(5.0), 1e-12);
    CHECK_NEAR(mean_bias({1.0, 3.0}, {1.0, 1.0}), 1.0, 1e-12);
    CHECK_NEAR(normalised_mean_bias({1.0, 3.0}, {1.0, 1.0}), 1.0, 1e-12);
}

void test_metrology() {
    const Matrix covariance = estimator_covariance({{1.0, 0.0}}, {{4.0, 1.0}, {1.0, 9.0}});
    CHECK_NEAR(covariance[0][0], 4.0, 1e-12);

    const Matrix structured =
        structured_measurement_covariance({{1.0}}, {{4.0}}, {{1.0}}, {{2.0}});
    CHECK_NEAR(structured[0][0], 7.0, 1e-12);

    const Matrix gls = gls_estimator({{1.0}, {2.0}}, {{1.0, 0.0}, {0.0, 1.0}});
    CHECK_NEAR(gls[0][0], 0.2, 1e-12);
    CHECK_NEAR(gls[0][1], 0.4, 1e-12);
    const Matrix ols = ols_estimator({{1.0}, {2.0}});
    CHECK_NEAR(ols[0][0], 0.2, 1e-12);

    const Matrix monte_carlo = monte_carlo_covariance({{1.0, 2.0}, {3.0, 4.0}});
    CHECK_NEAR(monte_carlo[0][0], 2.0, 1e-12);
    CHECK_NEAR(monte_carlo[0][1], 2.0, 1e-12);

    CHECK_NEAR(ndvi(0.2, 0.5), 0.42857142857142855, 1e-12);
    const std::vector<double> partials = ndvi_partials(0.5, 0.2);
    CHECK_NEAR(partials[0], 0.8163265306122449, 1e-12);
    CHECK_NEAR(partials[1], -2.0408163265306123, 1e-12);
    CHECK_NEAR(ndvi_uncertainty(0.5, 0.2, 0.01, 0.01, 0.0), 0.0219802, 1e-6);

    const std::vector<double> prediction = simple_kriging({2.0}, {{1.0}}, {{0.5}}, 0.0);
    CHECK_NEAR(prediction[0], 1.0, 1e-12);
    const Matrix kriging = kriging_error_covariance({{1.0}}, {{0.5}});
    CHECK_NEAR(kriging[0][0], 0.25, 1e-12);

    CHECK_NEAR(station_homogenisation_uncertainty(0.5, 4), 0.25, 1e-12);
    CHECK_NEAR(sampling_uncertainty(1.0, 30), 0.18257418583505536, 1e-12);
    CHECK_NEAR(bootstrap_mean_standard_error({5.0, 5.0, 5.0}, 20, 1), 0.0, 1e-12);
    CHECK_NEAR(simultaneous_nadir_bias({1.0, 3.0}, {0.0, 2.0}), 1.0, 1e-12);
    CHECK_NEAR(diurnal_harmonic(0.0, 10.0, 2.0, 0.0, 0.0, 0.0), 12.0, 1e-12);
}

void test_geostat_intercal() {
    CHECK_NEAR(spde_matern_variance(1.0, 1.0, 1, 1), 0.5, 1e-12);
    CHECK_NEAR(spde_matern_variance(2.0, 1.0, 2, 2), 1.0 / (16.0 * M_PI), 1e-12);
    CHECK_NEAR(circle_overlap_fraction(0.0, 1.0), 1.0, 1e-12);
    CHECK_NEAR(circle_overlap_fraction(2.0, 1.0), 0.0, 1e-12);
    const Matrix precision =
        spde_precision({{1.0, 0.0}, {0.0, 1.0}}, {{2.0, 0.0}, {0.0, 2.0}}, 1.0, 1.0, 1);
    CHECK_NEAR(precision[0][0], 3.0, 1e-12);
    CHECK_NEAR(precision[1][1], 3.0, 1e-12);
    CHECK_NEAR(gmrf_conditional_mean({{2.0, 0.0}, {0.0, 2.0}}, {2.0, 4.0})[1], 2.0, 1e-12);
    CHECK_NEAR(gmrf_conditional_covariance({{2.0, 0.0}, {0.0, 2.0}})[0][0], 0.5, 1e-12);

    const std::vector<double> gp_mean = gp_posterior_mean({0.0, 1.0}, {0.0, 1.0}, {0.5}, 1.0, 1.0, 0.0);
    CHECK_NEAR(gp_mean[0], 0.549318, 1e-5);
    const std::vector<double> gp_variance = gp_posterior_variance({0.0, 1.0}, {0.5}, 1.0, 1.0, 0.01);
    CHECK(gp_variance[0] > 0.0 && gp_variance[0] < 1.0);

    CHECK_NEAR(mad_envelope_statistic({0.0, 0.0}, {1.0, -1.0}, {1.0, 1.0}), 1.0, 1e-12);
    CHECK_NEAR(mad_envelope_p_value({0.0, 0.0}, {{0.0, 0.0}, {1.0, -1.0}}), 1.0, 1e-12);
    CHECK_NEAR(mad_envelope_p_value({10.0, 10.0}, {{0.0, 0.0}, {1.0, -1.0}}), 0.0, 1e-12);
    CHECK_NEAR(mad_global_statistic({1.0}, {0.0}, {2.0}), 0.5, 1e-12);

    CHECK_NEAR(gtwr_distance(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), std::sqrt(3.0), 1e-12);
    const std::vector<double> gtwr =
        gtwr_coefficients({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}}, {2.0, 3.0, 5.0},
                          {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0.0, 0.0, 1.0);
    CHECK_NEAR(gtwr[0], 2.0, 1e-12);
    CHECK_NEAR(gtwr[1], 3.0, 1e-12);
    CHECK_NEAR(gtwr_aicc(10, 1.0, 1.0), 34.093056378, 1e-6);

    CHECK_NEAR(spherical_variogram(0.0, 0.5, 2.0, 1.0), 0.5, 1e-12);
    CHECK_NEAR(spherical_variogram(1.0, 0.5, 2.0, 1.0), 2.5, 1e-12);
    CHECK_NEAR(exponential_variogram(0.0, 0.5, 2.0, 1.0), 0.5, 1e-12);
    CHECK_NEAR(gaussian_variogram(1.0, 0.5, 2.0, 1.0), 0.5 + 2.0 * (1.0 - std::exp(-1.0)),
               1e-12);
    CHECK_NEAR(matern_variogram(1.0, 0.5, 2.0, 1.0, 0.5), 0.5 + 2.0 * (1.0 - std::exp(-1.0)),
               1e-12);

    const Matrix weighted = weighted_covariance({{1.0, 0.0}, {0.0, 1.0}}, {1.0, 1.0});
    CHECK_NEAR(weighted[0][0], 0.5, 1e-12);
    const Matrix refined = weighted_pca_residual({{1.0, 0.0}, {-1.0, 0.0}}, {1.0, 1.0}, 1);
    CHECK_NEAR(refined[0][0], 0.0, 1e-12);
    CHECK_NEAR(refined[1][0], 0.0, 1e-12);

    CHECK_NEAR(gsics_leo_radiance(10.0, 1.0, 2.0), 4.5, 1e-12);
    CHECK_NEAR(gsics_uncertainty(10.0, 1.0, 2.0, 0.1, 0.2, 0.0), 0.45276925690687085, 1e-12);
    CHECK_NEAR(msu_delta_radiance(1.0, 2.0, 3.0, 1.0), 5.0, 1e-12);
    CHECK_NEAR(msu_nonlinearity(1.0, 0.5, 3.0, 1.0), 2.0, 1e-12);
    CHECK_NEAR(msu_z_factor(2.0, 0.9, 0.1, 0.2), 2.24, 1e-12);
    CHECK_NEAR(msu_calibrated_radiance(10.0, 5.0, 2.0, 2.24), 9.48, 1e-12);
    CHECK_NEAR(diurnal_drift_correction(10.0, 0.5, 2.0), 11.0, 1e-12);

    CHECK_NEAR(double_difference(10.0, 8.0, 9.0, 8.0), 1.0, 1e-12);
    CHECK_NEAR(predicted_brightness(8.0, 9.0, 8.0), 9.0, 1e-12);
    CHECK_NEAR(diurnal_anomaly(1.0, 2.0, 0.0, 6.0), 3.0, 1e-12);
    CHECK_NEAR(seasonal_amplitude(1.0, 2.0, 0.0, 3), 3.0, 1e-12);

    CHECK_NEAR(cwv_retrieval({200.0}, {1.0, 2.0}),
               1.0 + 2.0 * std::log(88.0 / 288.0), 1e-12);
    CHECK_NEAR(cwv_channel_sensitivity(200.0, 2.0), -2.0 / 88.0, 1e-12);
    CHECK_NEAR(sst_first_step({200.0}, {1.0, 2.0}, {0.5}),
               1.0 + 2.0 * std::log(88.0 / 288.0) +
                   0.5 * std::log(88.0 / 288.0) * std::log(88.0 / 288.0),
               1e-12);
    CHECK_NEAR(sst_channel_sensitivity(200.0, 2.0, 0.5),
               -(2.0 + std::log(88.0 / 288.0)) / 88.0, 1e-12);

    CHECK_NEAR(ssa_recalibration(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0), 3.0, 1e-12);
    CHECK_NEAR(ssa_derivative(0.0, 1.0, 0.0, 2.0, 3.0, 0.0), 7.0, 1e-12);

    CHECK_NEAR(clarreo_weight(1.0, 2.0, 3.0), 1.0 / 6.0, 1e-12);
    CHECK_NEAR(clarreo_bias({1.0, 3.0}, {1.0, 1.0}), 2.0, 1e-12);
    CHECK_NEAR(clarreo_uncertainty({1.0, 1.0}), 1.0 / std::sqrt(2.0), 1e-12);
    CHECK_NEAR(clarreo_temporal_variance(0.5, 4.0, 2.0), 2.0, 1e-12);
    CHECK_NEAR(clarreo_inflate(1.0, 3.0), 2.0, 1e-12);

    CHECK_NEAR(obcs_radiance(2.0, 10.0, 1.0, 0.5), 17.5, 1e-12);
    CHECK_NEAR(angular_standard_deviation({1.0, 3.0}), 1.0, 1e-12);
    CHECK_NEAR(bias_regression_slope({1.0, 2.0}, {1.0, 2.0}), 1.0, 1e-12);
    CHECK_NEAR(sno_variance_decomposition(1.0, 2.0, 3.0), 6.0, 1e-12);
    CHECK_NEAR(sno_noise_fraction(3.0, 6.0), 50.0, 1e-12);
    CHECK_NEAR(sno_t_statistic(1.0, 2.0, 4), 1.0, 1e-12);
}

void test_homogenisation() {
    CHECK_NEAR(absolute_mean_bias({1.0, 4.0}, {0.0, 0.0}), 2.5, 1e-12);
    CHECK_NEAR(factor_of_exceedance({2.0, 0.0}, {1.0, 1.0}), 0.0, 1e-12);
    CHECK_NEAR(relative_rmse({2.0, 4.0}, {4.0, 8.0}, {0.0, 0.0}), 0.5, 1e-12);
    CHECK_NEAR(acmant_objective({0.0, 2.0}, {1.0, 1.0}, {}, 1.0), 2.0, 1e-12);

    std::vector<double> mean, variance;
    homogenisation_ensemble({{1.0, 2.0}, {3.0, 4.0}}, mean, variance);
    CHECK_NEAR(mean[0], 2.0, 1e-12);
    CHECK_NEAR(variance[1], 2.0, 1e-12);

    CHECK_NEAR(residual_cusum({1.0, 1.0, 1.0, 1.0}, 2), 0.0, 1e-12);
    CHECK_NEAR(student_t_log_likelihood({0.0}, 3.0, 1.0),
               std::lgamma(2.0) - std::lgamma(1.5) - 0.5 * std::log(3.0 * M_PI), 1e-12);
    CHECK_NEAR(singh_maddala_density(1.0, 1.0, 1.0, 1.0), 0.25, 1e-12);
    CHECK_NEAR(sparse_jump_objective({{0.0}, {10.0}}, {0, 1}, {{0.0}, {10.0}}, 1.0), 1.0, 1e-12);
    CHECK_NEAR(cautious_correction(10.0, 2.0), 8.0, 1e-12);
    CHECK_NEAR(surface_urban_heat_island(5.0, 2.0), 3.0, 1e-12);
    CHECK_NEAR(excursion_probability(0.5, 0.0, 0.2), 0.9937903346742238, 1e-12);
    CHECK_NEAR(variational_lower_bound(2.0, 3.0), -1.0, 1e-12);

    CHECK_NEAR(msu_linear_radiance(1.0, 2.0, 5.0, 1.0), 9.0, 1e-12);
    CHECK_NEAR(msu_nonlinear_radiance(0.0, 2.0, 5.0, 1.0, 1.0, 1.0, 3.0), 10.0, 1e-12);
    CHECK_NEAR(msu_drift_term(2.0, 1.0, 3.0), 8.0, 1e-12);
    CHECK_NEAR(msu_drift_budget(1.0, 2.0, 1.0, 0.0, 1.0, 1.0, 5.0, 5.0), -3.0, 1e-12);
    CHECK_NEAR(single_difference(3.0, 1.0), 2.0, 1e-12);
    CHECK_NEAR(double_difference_simple(2.0, 1.0), 1.0, 1e-12);
    CHECK_NEAR(clarero_overlap_fraction(1.0, 1.0), 2.0 / 3.0 - std::sqrt(3.0) / (2.0 * M_PI),
               1e-12);
    CHECK_NEAR(crisp_cusum_statistic(6.0, 8.0, 2, 2), -1.0, 1e-12);
    CHECK_NEAR(crisp_cusum_statistic(1.0, 1.0, 0, 2), 0.0, 1e-12);
    CHECK(crisp_argmax({0.5, -2.0, 1.0}) == 1);
    CHECK_NEAR(msu_tlt_synthesis(1.0, 2.0, 3.0), 0.602, 1e-12);
    CHECK_NEAR(msu_z_counts(2.0, 5.0, 1.0, 3.0), 32.0, 1e-12);
    CHECK_NEAR(aci_susceptibility(0.5, 2.0), 1.0, 1e-12);
    CHECK_NEAR(nd_relative_bias(0.1, -0.5, 0.5), -0.45, 1e-12);
    CHECK_NEAR(lwp_error_variance(2.0, 3.0, 0.01, 0.04, -0.005), 0.34, 1e-12);
    CHECK_NEAR(dml_treatment_effect({1.0, 2.0}, {1.0, 2.0}), 1.0, 1e-12);
    CHECK_NEAR(g_computation_effect({1.0, 2.0, 3.0}, {0.2, 0.3, 0.5}), 2.3, 1e-12);
    CHECK_NEAR(natural_effect_total(1.0, 2.0), 3.0, 1e-12);
    CHECK_NEAR(mean_squared_error({1.0, 3.0}, {0.0, 0.0}), 5.0, 1e-12);
    CHECK_NEAR(coefficient_of_determination({1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}), 1.0, 1e-12);

    CHECK_NEAR(sno_scene_variance(1.0, 2.0), 3.0, 1e-12);
    CHECK_NEAR(clarreo_bias_variance({1.0, 1.0}), 0.5, 1e-12);
    CHECK_NEAR(truths_total_uncertainty(1.0, 2.0, 3.0), std::sqrt(14.0), 1e-12);
    double gain = 0.0, innovation = 0.0, innovation_variance = 0.0;
    const double drift = kalman_random_walk(0.0, 1.0, 0.5, 1.0, 2.0, 2.0, gain, innovation,
                                           innovation_variance);
    CHECK_NEAR(drift, 0.8571428571428571, 1e-12);
    CHECK_NEAR(gain, 0.42857142857142855, 1e-12);
    CHECK_NEAR(innovation_variance, 3.5, 1e-12);
    CHECK_NEAR(adversarial_loss({1.0}, {0.0}), 0.0, 1e-12);
    CHECK_NEAR(cycle_consistency_loss({1.0, 3.0}, {2.0, 5.0}), 1.5, 1e-12);
    CHECK_NEAR(cyclegan_objective(1.0, 2.0, 3.0, 1.0), 6.0, 1e-12);
}

void test_variable_validation() {
    CHECK_NEAR(standardised_difference_variance({1.0, -1.0}, {0.0, 0.0}, 1.0, 1.0, 0.0), 0.5,
               1e-12);
    CHECK_NEAR(quadruple_collocation_covariance(3.0, 1.0, 1.0, 2.0), 1.0, 1e-12);
    CHECK_NEAR(linear_ice_fraction(30.0, 100.0), 30.0, 1e-12);
    CHECK_NEAR(sea_ice_quadratic_bias(0.0, 0.1, 0.01, 0.5), 0.0525, 1e-12);
    CHECK_NEAR(sea_ice_bias_derivative(0.1, 0.01, 0.5), 0.11, 1e-12);
    CHECK_NEAR(inverse_variance_merge({1.0, 2.0}, {1.0, 4.0}), 1.2, 1e-12);
    CHECK_NEAR(precipitation_bias_sum(1.0, 2.0, 3.0), 6.0, 1e-12);
    CHECK_NEAR(observation_error_covariance({{2.0}}, {{3.0}})[0][0], 5.0, 1e-12);
    CHECK_NEAR(flux_error_decomposition(1.0, 2.0, 3.0), std::sqrt(14.0), 1e-12);
    CHECK_NEAR(ocean_colour_variance({1.0, 2.0}, {{1.0, 0.0}, {0.0, 1.0}}), 5.0, 1e-12);
    CHECK_NEAR(measurement_uncertainty(3.0, 4.0), 5.0, 1e-12);
    CHECK_NEAR(partially_linear_effect({1.0, 2.0}, {1.0, 2.0}), 1.0, 1e-12);
}

void test_regional_sign_and_clustering() {
    CHECK_NEAR(sign_difference_statistic(7, 3), 0.4, 1e-12);
    CHECK_NEAR(sign_test_percentage(0, 10), -100.0, 1e-12);
    CHECK_NEAR(energy_budget_trend(1.0, 2.0, 0.5, 0.5, 0.0), 2.0, 1e-12);
    CHECK_NEAR(budget_contribution(1.0, 2.0), 50.0, 1e-12);
    CHECK_NEAR(contribution_sensitivity(1.0, 2.0, true), 25.0, 1e-12);
    CHECK_NEAR(sign_test_percentage(10, 10), 0.0, 1e-12);
    CHECK_NEAR(pooling_standard_error_ratio(1.0, 3.0), 0.5, 1e-12);
    CHECK_NEAR(sign_test_binomial_p_value(7, 10), 0.171875, 1e-12);
    CHECK_NEAR(effective_sample_size(100.0, 0.5), 100.0 / 3.0, 1e-12);
    CHECK_NEAR(sign_test_empirical_p_value(7, {5, 7, 8, 9}), 0.75, 1e-12);
    CHECK_NEAR(north_test_standard_error(4.0, 100.0), 4.0 * std::sqrt(0.02), 1e-12);
    CHECK_NEAR(normal_quantile(0.975), 1.959963984540054, 1e-9);
    CHECK_NEAR(sign_test_critical_value(100.0, 0.5, 0.05), 1.6448536269514722 / std::sqrt(100.0 / 3.0),
               1e-12);
    CHECK_NEAR(sign_test_critical_value_derivative(100.0, 0.5, 0.05),
               normal_quantile(0.95) / 10.0 / (0.3535533905932738 * 1.224744871391589),
               1e-12);

    CHECK_NEAR(zonal_mean_at_latitude({1.0, 2.0, 3.0, 4.0}, {0, 0, 1, 1}, 0), 1.5, 1e-12);
    const std::vector<double> asymmetric =
        zonally_asymmetric_trend({1.0, 2.0, 3.0, 4.0}, {0, 0, 1, 1});
    CHECK_NEAR(asymmetric[0], -0.5, 1e-12);
    CHECK_NEAR(asymmetric[3], 0.5, 1e-12);
    CHECK_NEAR(regional_mean({1.0, 2.0, 3.0, 4.0}, {1, 0, 1, 0}), 2.0, 1e-12);
    CHECK_NEAR(sign_agreement({1.0, -1.0, 1.0}, {1.0, 1.0, 1.0}, {1, 1, 1}), 2.0 / 3.0, 1e-12);

    const double support = 1.1;
    CHECK_NEAR(gev_location_score(1.0, 0.0, 1.0, 0.1),
               (1.1 / support - std::pow(support, -11.0)), 1e-12);
    CHECK_NEAR(gev_trend_score(1.0, 0.0, 1.0, 0.1, 2.0),
               2.0 * (1.1 / support - std::pow(support, -11.0)), 1e-12);
    CHECK_NEAR(shrinkage_log_prior(1.0, 0.0, 2.0), -0.25, 1e-12);
    CHECK_NEAR(shrinkage_derivative(1.0, 0.0, 2.0), -0.5, 1e-12);
    CHECK_NEAR(pooled_posterior_variance(0.4, 0.6), 0.24, 1e-12);
    CHECK_NEAR(pooled_cluster_mean({1.0, 3.0}, {4.0, 1.0}), 1.4, 1e-12);

    const EigenSystem eigen = symmetric_eigen({{2.0, 1.0}, {1.0, 2.0}});
    CHECK_NEAR(eigen.values[0], 3.0, 1e-12);
    CHECK_NEAR(eigen.values[1], 1.0, 1e-12);
    const Matrix correlation = temporal_correlation_matrix({{1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}});
    CHECK_NEAR(correlation[0][1], 1.0, 1e-12);
    const Matrix loadings = {{0.9, 0.1}, {0.8, 0.2}, {0.1, 0.9}, {0.2, 0.8}};
    Matrix rotation;
    const Matrix rotated = varimax_rotation(loadings, 100, 1e-10, rotation);
    CHECK(rotation.size() == 2);
    CHECK_NEAR(rotation[0][0] * rotation[0][0] + rotation[1][0] * rotation[1][0], 1.0, 1e-9);
    CHECK_NEAR(rotation[0][0] * rotation[0][1] + rotation[1][0] * rotation[1][1], 0.0, 1e-9);
    CHECK_NEAR(rotated[0][0] * rotated[0][0] + rotated[0][1] * rotated[0][1],
               loadings[0][0] * loadings[0][0] + loadings[0][1] * loadings[0][1], 1e-9);
    const Matrix identity = {{1.0, 0.0}, {0.0, 1.0}};
    const std::vector<double> rpc = reof_rpc_trends({{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}},
                                                    identity, identity);
    CHECK_NEAR(rpc[0], 1.0, 1e-12);
    CHECK_NEAR(reof_back_projection({1.0, 2.0}, {{1.0, 1.0}}, identity, 0), 3.0, 1e-12);

    CHECK_NEAR(weatherhead_trend_standard_error(1.0, 0.5, 100.0), std::sqrt(3.0) / 10.0, 1e-12);
    CHECK_NEAR(detection_time(1.0, 0.5, 0.1, 0.05), 10.485, 1e-3);
    CHECK_NEAR(detection_time_autocorrelation_derivative(1.0, 0.5, 0.1, 0.05), 9.31937, 1e-4);
    CHECK_NEAR(weatherhead_se_autocorrelation_derivative(1.0, 0.5, 100.0), 0.23094, 1e-4);
    CHECK_NEAR(weatherhead_se_noise_derivative(1.0, 0.5, 100.0), 0.173205, 1e-4);
    CHECK_NEAR(detection_time_noise_derivative(1.0, 0.5, 0.1, 0.05), 6.99, 1e-2);
    CHECK_NEAR(trend_pattern_correlation({1.0, 2.0, 3.0}, {-1.0, -2.0, -3.0}), -1.0, 1e-12);
    CHECK_NEAR(effective_spatial_dof(100.0, 0.5), 100.0 / 3.0, 1e-12);
    CHECK_NEAR(correlation_t_statistic(0.5, 100.0 / 3.0), 3.23179, 1e-4);
    CHECK_NEAR(coherence_height({100.0, 70.0, 50.0}, {0.8, 0.2, 0.1}, 0.3), 70.0, 1e-12);
    CHECK_NEAR(segmented_trend_difference(0.2, -0.75), 0.95, 1e-12);
    CHECK_NEAR(bayesian_posterior(0.2, 0.1, 0.3), 0.06 / 0.13, 1e-12);
    CHECK_NEAR(bayesian_log_odds(0.2, 0.1, 0.3), std::log(3.0 / 7.0) + std::log(2.0), 1e-12);
    CHECK_NEAR(posterior_prior_derivative(0.2, 0.1, 0.3), 1.18343, 1e-5);
    CHECK_NEAR(probability_weighted_centroid({0.0, 10.0}, {1.0, 1.0}), 5.0, 1e-12);
    CHECK_NEAR(centroid_probability_sensitivity(10.0, 5.0, 2.0), 2.5, 1e-12);
    CHECK_NEAR(kolmogorov_smirnov_d_plus({1.0, 2.0}, {0.0, 1.0}), 0.5, 1e-12);
    CHECK_NEAR(optical_trend_precision(0.5, {1.0, 1.0}, 1.0, 0.0), 2.5, 1e-12);
    CHECK_NEAR(quantile_influence(2.0, 1.0, 2.0, 0.9, true, 2.0), 0.225, 1e-12);
    CHECK_NEAR(fused_trend({0.5, 0.5}, {1.0, 3.0}), 2.0, 1e-12);
    CHECK_NEAR(fused_weight_sensitivity({1.0, 3.0}, 2.0, 0), -1.0, 1e-12);
    CHECK_NEAR(vertical_trend_difference(1.0, -1.0), 2.0, 1e-12);
    CHECK_NEAR(trend_signal_to_noise(2.0, 1.0), 2.0, 1e-12);
    CHECK_NEAR(potts_energy({0, 0, 0, 0}, 2, 2, {0.0}, {1.0}), 4.0, 1e-12);
    CHECK_NEAR(parallel_tempering_swap_log_ratio(10.0, 8.0, 1.0, 2.0), -1.0, 1e-12);
    std::vector<int> exchange_labels(9, 0);
    std::vector<double> exchange_alpha{0.0}, exchange_beta{1.0};
    std::mt19937 exchange_generator(7);
    const bool accepted = potts_exchange_update(
        exchange_labels, 3, 3, exchange_alpha, exchange_beta,
        [](const std::vector<int>& values) {
            for (int value : values)
                if (value != 0) return -1e6;
            return 0.0;
        },
        0.0, exchange_generator);
    CHECK(accepted);
    for (int label : exchange_labels) CHECK(label == 0);
    CHECK_NEAR(elevation_trend_slope({1.0, 2.0, 3.0}, {2.0, 3.0, 4.0}), 1.0, 1e-12);
    CHECK_NEAR(elevation_leverage(3.0, {1.0, 2.0, 3.0}), 0.5, 1e-12);
    CHECK_NEAR(latitude_corrected_edw(1.0, 0.4, 0.5), 0.8, 1e-12);
    CHECK_NEAR(lapse_rate_trend(2.0, 1.0, 10.0), 0.1, 1e-12);
}

void test_lagrangian_objects() {
    const std::pair<double, double> eigenvalues = hessian_eigenvalues(-2.0, 0.0, -3.0);
    CHECK_NEAR(eigenvalues.first, -2.0, 1e-12);
    CHECK_NEAR(eigenvalues.second, -3.0, 1e-12);
    CHECK_NEAR(shape_index(-1.0, -2.0), 0.7951672353008665, 1e-12);
    CHECK_NEAR(shape_index(-2.0, -2.0), 1.0, 1e-12);
    CHECK_NEAR(shape_index(2.0, 2.0), -1.0, 1e-12);
    CHECK_NEAR(shape_index_k1_derivative(-1.0, -2.0), 0.25464790894703254, 1e-12);
    CHECK_NEAR(shape_index_k2_derivative(-1.0, -2.0), -0.12732395447351627, 1e-12);
    const Matrix smoothed = gaussian_smooth({{5.0, 5.0, 5.0}, {5.0, 5.0, 5.0}, {5.0, 5.0, 5.0}}, 1.0);
    CHECK_NEAR(smoothed[1][1], 5.0, 1e-12);
    CHECK_NEAR(linear_trend({0.0, 2.0, 4.0}), 2.0, 1e-12);
    CHECK_NEAR(object_frequency({1, 0, 1}), 2.0, 1e-12);
    CHECK_NEAR(mean_object_intensity({5.0, 7.0}, {1, 0}), 5.0, 1e-12);
    CHECK_NEAR(nearest_link_distance(0.0, 0.0, 3.0, 4.0, 6.0), 5.0, 1e-12);
    CHECK_NEAR(nearest_link_distance(0.0, 0.0, 3.0, 4.0, 4.0), -1.0, 1e-12);
    CHECK_NEAR(object_population_trend({1.0, 2.0, 3.0}), 1.0, 1e-12);
    CHECK_NEAR(along_track_trend({5.0, 3.0, 1.0}), -2.0, 1e-12);
}

void test_detectability_power() {
    CHECK_NEAR(variance_inflation_factor(0.5), 3.0, 1e-12);
    CHECK_NEAR(stationary_noise_variance(1.0, 0.5), 4.0 / 3.0, 1e-12);
    CHECK_NEAR(minimum_detectable_trend_power(1.0, 0.5, 100.0, 0.05, 0.8), 0.00485247, 1e-7);
    CHECK_NEAR(detection_time_power(1.0, 0.5, 0.1, 0.05, 0.8), 13.304, 1e-2);
    CHECK_NEAR(trend_power(0.2, 0.1, 0.05), 0.51601, 1e-4);
    CHECK_NEAR(signal_to_noise_ratio(0.1, 100.0, 1.0), 10.0, 1e-12);
    CHECK_NEAR(minimum_detectable_trend_from_snr(1.0, 0.2, 50.0), 0.004, 1e-12);
    CHECK_NEAR(mdt_autocorrelation_derivative(1.0, 0.5, 100.0, 0.05, 0.8), 0.00646996, 1e-7);
    CHECK_NEAR(mdt_noise_derivative(1.0, 0.5, 100.0, 0.05, 0.8), 0.00485247, 1e-7);
    CHECK_NEAR(detection_time_trend_derivative(1.0, 0.5, 0.1, 0.05, 0.8), -88.69, 1e-1);
    CHECK(variable_is_feasible(0.1, 0.2));
    CHECK(!variable_is_feasible(0.3, 0.2));
    CHECK_NEAR(effective_field_significance(0.05, 0.5, 10.0), 0.15718, 1e-4);
}

void test_quantile_clustering() {
    CHECK_NEAR(check_loss(2.0, 0.9), 1.8, 1e-12);
    CHECK_NEAR(check_loss(-2.0, 0.9), 0.2, 1e-12);
    CHECK_NEAR(sample_quantile({1.0, 2.0, 3.0, 4.0}, 0.5), 2.5, 1e-12);
    const QuantileFit fit = quantile_trend({0.0, 1.0, 2.0, 3.0}, {1.0, 3.0, 5.0, 7.0}, 0.5);
    CHECK_NEAR(fit.slope, 2.0, 1e-6);
    CHECK_NEAR(fit.intercept, 1.0, 1e-6);
    CHECK_NEAR(ald_density(0.0, 0.5, 1.0), 0.25, 1e-12);
    CHECK_NEAR(ald_subgradient(-1.0, 0.9), -0.1, 1e-12);
    CHECK_NEAR(ald_subgradient(1.0, 0.9), 0.9, 1e-12);
    const std::vector<double> conditional =
        potts_conditional_probabilities({2.0, 0.0}, {0.0, 0.0}, {1.0, 1.0});
    CHECK_NEAR(conditional[0], std::exp(2.0) / (std::exp(2.0) + 1.0), 1e-12);
    CHECK_NEAR(potts_log_odds_neighbour_derivative(1.5), 1.5, 1e-12);
    std::vector<int> labels(9, 0);
    std::mt19937 generator(7);
    potts_gibbs_sweep(labels, 3, 3, {0.0, 0.0}, {20.0, 20.0}, generator);
    for (int label : labels) CHECK(label == 0);
    potts_swendsen_wang_step(labels, 3, 3, {0.0, 0.0}, {20.0, 20.0}, generator);
    for (int label : labels) CHECK(label == labels[0]);
}

void test_null_ensembles_and_flow() {
    const Matrix correlation = spatial_correlation_matrix(2, 2, 0.5);
    CHECK_NEAR(correlation[0][1], 0.5, 1e-12);
    CHECK_NEAR(correlation[0][3], std::pow(0.5, std::sqrt(2.0)), 1e-12);
    Matrix factor;
    CHECK(cholesky_factor({{4.0, 2.0}, {2.0, 3.0}}, factor));
    CHECK_NEAR(factor[0][0], 2.0, 1e-12);
    CHECK_NEAR(factor[1][1], std::sqrt(2.0), 1e-12);
    Matrix field = simulate_zero_trend_field(2, 2, 4, 0.5, 0.5, 11);
    CHECK(field.size() == 4 && field[0].size() == 4);
    CHECK(!std::isnan(field[3][3]));
    Matrix zeros(5, std::vector<double>(4, 0.0));
    inject_uniform_trend(zeros, 100.0);
    CHECK(positive_trend_count(zeros) == 4);
    const std::vector<std::size_t> nulls =
        sign_test_null_ensemble(3, 3, 8, 0.4, 0.3, 5, 3);
    CHECK(nulls.size() == 5);
    CHECK_NEAR(central_difference([](double x) { return x * x; }, 3.0, 0.1), 6.0, 1e-12);

    auto strain = [](double x, double y, double) {
        return std::make_pair(0.1 * x, -0.1 * y);
    };
    CHECK_NEAR(ftle(0.0, 0.0, 0.0, 1.0, strain, 1e-4, 20), 0.1, 1e-6);
    const Matrix ftle_values = ftle_field(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 3, 3, strain, 1e-4, 10);
    CHECK(ftle_values.size() == 3 && ftle_values[0].size() == 3);
    auto uniform = [](double, double, double) { return std::make_pair(1.0, 0.0); };
    CHECK_NEAR(lagrangian_descriptor(0.0, 0.0, 0.0, 1.0, uniform, 10), 2.0, 1e-12);
    CHECK_NEAR(occurrence_percent(3, 10), 30.0, 1e-12);
    CHECK_NEAR(regime_shift_difference(2.0, 1.0), 1.0, 1e-12);
    CHECK_NEAR(year_to_year_variability({1.0, 2.0, 3.0}), 1.0, 1e-12);
    CHECK_NEAR(wilks_correlation(2.0, -1.0, 1.0), std::exp(-2.0), 1e-12);
    const Matrix occurrence = simulate_occurrence_field(2, 2, 3, 0.5, 1.0, 5);
    CHECK_NEAR(occurrence[2][1], 1.0, 1e-12);
    const Matrix no_occurrence = simulate_occurrence_field(2, 2, 3, 0.5, 0.0, 5);
    CHECK_NEAR(no_occurrence[1][0], 0.0, 1e-12);
    const Matrix gamma_field = simulate_gamma_field(2, 2, 2.0, 1.0, 9);
    CHECK(gamma_field[0][0] > 0.0);
    Matrix trended = {{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}};
    detrend_field(trended);
    CHECK_NEAR(linear_trend({trended[0][0], trended[1][0], trended[2][0]}), 0.0, 1e-12);
    CHECK_NEAR(linear_trend({trended[0][1], trended[1][1], trended[2][1]}), 0.0, 1e-12);
    const std::vector<double> smoothed = cyclic_moving_average({0.0, 10.0, 0.0}, 3);
    CHECK_NEAR(smoothed[0], 10.0 / 3.0, 1e-12);
    CHECK(can_stitch(3.0, 5.0, 4.0, 6.0));
    CHECK(!can_stitch(5.0, 5.0, 4.0, 6.0));
    const auto tracks = stitch_tracks({{{{0.0, 0.0}}}, {{{1.0, 0.0}}}, {{{10.0, 10.0}}}}, 2.0, 1.0);
    CHECK(tracks.size() == 2);
    CHECK(tracks[0].size() == 2);
    CHECK_NEAR(jet_trend_correction(1.0, 0.5, 2.0), 2.0, 1e-12);
    const Matrix constant(3, std::vector<double>(3, 5.0));
    CHECK_NEAR(sobel_gradient(constant, 0)[1][1], 0.0, 1e-12);
    CHECK_NEAR(total_sobel_signal(0.0, 3.0, 4.0, 1.0), 5.0, 1e-12);
    CHECK_NEAR(total_sobel_signal(3.0, 4.0, 0.0, 1.0), 7.0, 1e-12);
    CHECK_NEAR(calibration_factor(2.0, 4.0), 0.5, 1e-12);
    const std::vector<double> suppressed = non_maximum_suppression_1d({1.0, 3.0, 2.0});
    CHECK_NEAR(suppressed[0], 0.0, 1e-12);
    CHECK_NEAR(suppressed[1], 3.0, 1e-12);
    const std::vector<int> hysteresis = hysteresis_threshold_1d({0.9, 0.6, 0.2}, 0.8, 0.5);
    CHECK(hysteresis[0] == 1 && hysteresis[1] == 2 && hysteresis[2] == 0);
    CHECK_NEAR(cohens_d(2.0, 10.0, 2.0), 4.0, 1e-12);
    CHECK(connected_component_count({1, 1, 0, 0}, 2, 1, 2) == 1);
    CHECK(connected_component_count({1, 0, 0, 1}, 2, 1, 2) == 2);
    const Matrix ridge_field = {{5.0, 5.0, 5.0, 5.0, 5.0},
                                {5.0, 8.0, 6.0, 7.0, 5.0},
                                {5.0, 5.0, 5.0, 5.0, 5.0}};
    const auto ridge = constrained_ridge_path(ridge_field, 2, 1.0);
    CHECK(ridge.size() == 2);
    const Matrix closed = {{1100.0, 1100.0, 1100.0},
                           {1100.0, 900.0, 1100.0},
                           {1100.0, 1100.0, 1100.0}};
    CHECK(has_closed_contour(closed, 1, 1, 200.0, 1));
    const Matrix open = {{1000.0, 1000.0, 1000.0},
                         {1000.0, 900.0, 1000.0},
                         {1000.0, 1000.0, 1000.0}};
    CHECK(!has_closed_contour(open, 1, 1, 200.0, 1));
}

void test_ocean_profiles() {
    CHECK_NEAR(dynamic_class_mean(2.0, 1.0, 3.0, 0.5), 5.5, 1e-12);
    CHECK_NEAR(static_region_mixture(0.5, 2.0, 0.5, -2.0), 0.0, 1e-12);
    CHECK_NEAR(class_switch_level_shift(3.0, 1.0), 2.0, 1e-12);
    CHECK_NEAR(depth_warming_gradient(0.88, 0.64, 35.0), 0.00685714, 1e-8);
    const auto events = marine_heatwave_events({1.0, 2.0, 2.0, 0.1, 2.0, 2.0},
                                               {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, 2);
    CHECK(events.size() == 2);
    CHECK(events[0].first == 0 && events[0].second == 2);
    CHECK(events[1].first == 4 && events[1].second == 5);
    CHECK_NEAR(marine_heatwave_intensity({2.0, 4.0}, {1.0, 1.0}, 0, 1), 2.0, 1e-12);
    CHECK_NEAR(marine_heatwave_frequency(5, 2.0), 2.5, 1e-12);
    CHECK_NEAR(mean_event_duration({{0, 2}, {4, 5}}), 2.5, 1e-12);
    CHECK_NEAR(vertically_integrated_anomaly({1.0, 2.0, 3.0}, {10.0, 10.0, 10.0}), 60.0, 1e-12);
    CHECK_NEAR(marine_heatwave_cooccurrence({1, 0, 1}, {1, 0, 0}), 1.0 / 3.0, 1e-12);
    CHECK_NEAR(independence_delta({1, 0, 1}, {1, 0, 0}), 1.0 / 3.0, 1e-12);
    CHECK(is_subsurface_intensified(50.0, 30.0));
    CHECK(!is_subsurface_intensified(20.0, 30.0));
    CHECK(depth_of_maximum_anomaly({1.0, 5.0, 3.0}) == 1);

    Gaussian2D first;
    first.mean_x = 0.0;
    first.mean_y = 0.0;
    Gaussian2D second;
    second.mean_x = 2.0;
    second.mean_y = 0.0;
    CHECK_NEAR(gaussian_2d_density(first, 0.0, 0.0), 1.0 / (2.0 * M_PI), 1e-12);
    CHECK_NEAR(pcm_posterior(first, second, 1.0, 1.0, 1.0, 0.0), 0.5, 1e-12);
    const std::pair<double, double> gradient =
        pcm_posterior_gradient(first, second, 1.0, 1.0, 1.0, 0.0);
    CHECK_NEAR(gradient.first, -0.5, 1e-12);
    CHECK_NEAR(cluster_weighted_mean({2.0, 4.0}, {1.0, 3.0}), 3.5, 1e-12);
    CHECK_NEAR(kmeans_objective({{0.0}, {2.0}, {5.0}}, {0, 0, 1}, {{1.0}, {5.0}}), 2.0, 1e-12);
    CHECK_NEAR(cluster_heat_content({1.0, 1.0}, {1.0, 0.0}, {10.0, 20.0}, 1.0), 10.0, 1e-12);
    CHECK_NEAR(cluster_area({1.0, 2.0, 3.0}, {0, 1, 0}, 0), 4.0, 1e-12);
    const std::vector<double> coefficients =
        ridge_fingerprint({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}}, {1.0, 2.0, 3.0}, 0.0);
    CHECK_NEAR(coefficients[0], 1.0, 1e-12);
    CHECK_NEAR(coefficients[1], 2.0, 1e-12);
    CHECK_NEAR(attribution_ambiguity(2.0, 2.0), 0.5, 1e-12);
    CHECK_NEAR(tracer_anomaly_decomposition(0.1, 10.0, 0.5, 2.0, 1.0), 3.0, 1e-12);
    const std::vector<double> retrieved =
        retrieved_trend_profile({{1.0, 0.0}, {0.5, 0.5}}, {4.0, 0.0});
    CHECK_NEAR(retrieved[0], 4.0, 1e-12);
    CHECK_NEAR(retrieved[1], 2.0, 1e-12);
    const Matrix propagated = trend_covariance({{1.0, 0.0}, {0.5, 0.5}},
                                               {{1.0, 0.0}, {0.0, 2.0}},
                                               {{0.0, 0.0}, {0.0, 0.0}});
    CHECK_NEAR(propagated[0][0], 1.0, 1e-12);
    CHECK_NEAR(propagated[0][1], 0.5, 1e-12);
    CHECK_NEAR(propagated[1][1], 0.75, 1e-12);
    CHECK_NEAR(toe_relative_error(12.0, 10.0, 10.0), 0.2, 1e-12);
    CHECK_NEAR(kernel_sensitivity({{1.0, 0.0}, {0.5, 0.5}}, 1, 0), 0.5, 1e-12);
    const std::vector<double> apriori =
        apriori_trend_sensitivity({{1.0, 0.0}, {0.5, 0.5}}, {1.0, 2.0});
    CHECK_NEAR(apriori[0], 0.0, 1e-12);
    CHECK_NEAR(apriori[1], 0.5, 1e-12);
}

int main() {
    test_detectability_power();
    test_quantile_clustering();
    test_null_ensembles_and_flow();
    test_ocean_profiles();
    test_variable_validation();
    test_regional_sign_and_clustering();
    test_lagrangian_objects();
    test_homogenisation();
    test_geostat_intercal();
    test_metrology();
    test_magnitude_hierarchy();
    test_spatial_advanced();
    test_regimes_and_causal();
    test_copula_hmm();
    test_attribution_and_ml();
    test_change_detection();
    test_advanced_ozone_landuse();
    test_resolved_gaps();
    test_aerosol_ozone();
    test_new_gaps();
    test_cryo();
    test_pde();
    test_linalg();
    test_stats();
    test_csv_and_json();
    test_fft();
    test_loess_and_stl();
    test_harmonic();
    test_trend();
    test_retrieval();
    test_series_and_analysis();
    test_catalog_and_intent();
    test_co2_analysis();
    test_amoc();
    test_netcdf();
    test_physics();
    test_earth_land_spatial();
    test_atmos_ocean();
    test_closures();
    test_fetch();

    if (failures == 0) {
        std::cout << "all tests passed\n";
        return 0;
    }
    std::cerr << failures << " check(s) failed\n";
    return 1;
}
