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

#include "analysis.hpp"
#include "catalog.hpp"
#include "common.hpp"
#include "csv.hpp"
#include "fft.hpp"
#include "harmonic.hpp"
#include "fetch.hpp"
#include "intent.hpp"
#include "json.hpp"
#include "linalg.hpp"
#include "netcdf.hpp"
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

}  // namespace

int main() {
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
    test_fetch();

    if (failures == 0) {
        std::cout << "all tests passed\n";
        return 0;
    }
    std::cerr << failures << " check(s) failed\n";
    return 1;
}
