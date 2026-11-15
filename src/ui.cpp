#include "ui.hpp"

#include <arpa/inet.h>
#include <netinet/in.h>
#include <signal.h>
#include <sys/socket.h>
#include <unistd.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "analysis.hpp"
#include "catalog.hpp"
#include "common.hpp"
#include "csv.hpp"
#include "fetch.hpp"
#include "formula_lab.hpp"
#include "infotheory.hpp"
#include "intent.hpp"
#include "json.hpp"
#include "physics.hpp"
#include "retrieval.hpp"
#include "spatial.hpp"
#include "sstunc.hpp"
#include "trend.hpp"
#include "ui_page.hpp"

namespace climate {
namespace {

struct BadRequest : std::runtime_error {
    explicit BadRequest(const std::string& message) : std::runtime_error(message) {}
};

struct Request {
    std::string method;
    std::string path;
    std::string query;
    std::string body;
    std::map<std::string, std::string> headers;
};

struct Response {
    int status = 200;
    std::string content_type = "application/json";
    std::string body;
};

std::string status_text(int status) {
    switch (status) {
        case 200: return "OK";
        case 204: return "No Content";
        case 400: return "Bad Request";
        case 404: return "Not Found";
        case 405: return "Method Not Allowed";
        case 413: return "Payload Too Large";
        case 500: return "Internal Server Error";
        default: return "OK";
    }
}

std::string lower_copy(std::string text) {
    for (auto& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return text;
}

std::string url_decode(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '+') {
            out.push_back(' ');
        } else if (text[i] == '%' && i + 2 < text.size()) {
            const auto hex = [](char c) -> int {
                if (c >= '0' && c <= '9') return c - '0';
                if (c >= 'a' && c <= 'f') return c - 'a' + 10;
                if (c >= 'A' && c <= 'F') return c - 'A' + 10;
                return -1;
            };
            const int high = hex(text[i + 1]);
            const int low = hex(text[i + 2]);
            if (high >= 0 && low >= 0) {
                out.push_back(static_cast<char>(high * 16 + low));
                i += 2;
            } else {
                out.push_back(text[i]);
            }
        } else {
            out.push_back(text[i]);
        }
    }
    return out;
}

std::string query_parameter(const std::string& query, const std::string& key) {
    std::size_t begin = 0;
    while (begin <= query.size()) {
        const std::size_t end = query.find('&', begin);
        const std::string pair = query.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
        const std::size_t equals = pair.find('=');
        if (equals != std::string::npos && pair.substr(0, equals) == key)
            return url_decode(pair.substr(equals + 1));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return std::string();
}

std::string read_text_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw BadRequest("cannot open " + path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

json::Value number_array(const std::vector<double>& values) {
    json::Value array = json::Value::array();
    for (double value : values) array.push(json::Value::from_number(value));
    return array;
}

json::Value string_array(const std::vector<std::string>& values) {
    json::Value array = json::Value::array();
    for (const auto& value : values) array.push(json::Value::from_string(value));
    return array;
}

double value_number(const json::Value& body, const std::string& key, double fallback) {
    const json::Value* found = body.find(key);
    return found == nullptr ? fallback : found->as_number(fallback);
}

std::string value_string(const json::Value& body, const std::string& key,
                         const std::string& fallback = "") {
    const json::Value* found = body.find(key);
    return found == nullptr ? fallback : found->as_string(fallback);
}

bool value_bool(const json::Value& body, const std::string& key, bool fallback) {
    const json::Value* found = body.find(key);
    return found == nullptr ? fallback : found->as_bool(fallback);
}

std::vector<double> numeric_list(const json::Value* value) {
    std::vector<double> out;
    if (value == nullptr) return out;
    if (value->is_array()) {
        for (const auto& item : value->items) {
            if (item.is_number() && std::isfinite(item.number)) out.push_back(item.number);
        }
        return out;
    }
    if (!value->is_string()) return out;
    std::string text = value->text;
    for (char& c : text) {
        if (c == ',' || c == ';' || c == '\n' || c == '\t') c = ' ';
    }
    std::istringstream stream(text);
    double number = 0.0;
    while (stream >> number) {
        if (std::isfinite(number)) out.push_back(number);
    }
    return out;
}

json::Value mk_json(const MKResult& mk) {
    json::Value value = json::Value::object();
    value.set("n", json::Value::from_number(static_cast<double>(mk.n)));
    value.set("s", json::Value::from_number(mk.s));
    value.set("variance", json::Value::from_number(mk.variance));
    value.set("variance_corrected", json::Value::from_number(mk.variance_corrected));
    value.set("correction_factor", json::Value::from_number(mk.correction_factor));
    value.set("z", json::Value::from_number(mk.z));
    value.set("p_value", json::Value::from_number(mk.p_value));
    value.set("tau", json::Value::from_number(mk.tau));
    value.set("critical_value", json::Value::from_number(mk.critical_value));
    value.set("significant", json::Value::from_bool(mk.significant));
    value.set("direction", json::Value::from_string(mk.direction));
    return value;
}

json::Value sen_json(const SenResult& sen) {
    json::Value value = json::Value::object();
    value.set("slope", json::Value::from_number(sen.slope));
    value.set("intercept", json::Value::from_number(sen.intercept));
    value.set("lower", json::Value::from_number(sen.lower));
    value.set("upper", json::Value::from_number(sen.upper));
    value.set("alpha", json::Value::from_number(sen.alpha));
    value.set("pairs", json::Value::from_number(static_cast<double>(sen.pairs)));
    value.set("computed", json::Value::from_bool(sen.computed));
    return value;
}

json::Value ensemble_json(const TrendEnsemble& ensemble) {
    json::Value value = json::Value::object();
    json::Value members = json::Value::array();
    for (const auto& member : ensemble.members) {
        json::Value item = json::Value::object();
        item.set("name", json::Value::from_string(member.name));
        item.set("formula", json::Value::from_string(member.formula));
        item.set("slope", json::Value::from_number(member.slope));
        item.set("intercept", json::Value::from_number(member.intercept));
        item.set("rmse", json::Value::from_number(member.rmse));
        item.set("p_value", json::Value::from_number(member.p_value));
        item.set("score", json::Value::from_number(member.score));
        item.set("significant", json::Value::from_bool(member.significant));
        item.set("computed", json::Value::from_bool(member.computed));
        item.set("selected", json::Value::from_bool(member.selected));
        members.push(item);
    }
    value.set("members", members);
    value.set("selected", json::Value::from_number(static_cast<double>(ensemble.selected)));
    value.set("agreement", json::Value::from_number(ensemble.agreement));
    value.set("robust_slope", json::Value::from_number(ensemble.robust_slope));
    value.set("robust_p_value", json::Value::from_number(ensemble.robust_p_value));
    value.set("consensus", json::Value::from_string(ensemble.consensus));
    return value;
}

struct SeriesInput {
    std::vector<double> time;
    std::vector<double> values;
    std::string variable = "series";
    std::string units;
};

SeriesInput series_from_body(const json::Value& body) {
    std::string csv;
    if (body.has("file")) {
        const std::string name = value_string(body, "file");
        if (name.empty() || name.find('/') != std::string::npos ||
            name.find('\\') != std::string::npos || name.find("..") != std::string::npos)
            throw BadRequest("invalid data file name");
        csv = read_text_file("data/" + name);
    } else {
        csv = value_string(body, "csv");
    }
    if (csv.empty()) throw BadRequest("provide a CSV body or a server-side file");

    const Table table = parse_csv_text(csv);
    const std::string time_name = value_string(body, "time");
    const std::string value_name = value_string(body, "value");
    int time_column = time_name.empty() ? 0 : table.column(time_name);
    if (time_name.empty() && table.header.empty()) time_column = 0;
    if (time_column < 0) throw BadRequest("time column '" + time_name + "' was not found");

    int value_column = -1;
    if (!value_name.empty()) {
        value_column = table.column(value_name);
        if (value_column < 0) throw BadRequest("value column '" + value_name + "' was not found");
    } else {
        value_column = table.column("value");
        if (value_column < 0 && table.header.size() > 1) value_column = 1;
        if (value_column < 0) value_column = 0;
    }

    std::vector<std::pair<double, double>> rows;
    rows.reserve(table.rows.size());
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        double time = 0.0;
        if (!parse_number(table.cell(i, time_column), time)) continue;
        double value = std::numeric_limits<double>::quiet_NaN();
        parse_number(table.cell(i, value_column), value);
        rows.emplace_back(time, value);
    }
    if (rows.size() < 3) throw BadRequest("need at least three numeric observations");
    std::stable_sort(rows.begin(), rows.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });

    SeriesInput input;
    input.variable = value_name.empty() ? (table.header.size() > 1 ? table.header[1] : "series")
                                        : value_name;
    input.units = value_string(body, "units");
    input.time.reserve(rows.size());
    input.values.reserve(rows.size());
    for (const auto& row : rows) {
        input.time.push_back(row.first);
        input.values.push_back(row.second);
    }
    return input;
}

Response json_response(json::Value value, int status = 200) {
    Response response;
    response.status = status;
    response.body = json::dump(value);
    return response;
}

Response analyze_series_request(const json::Value& body) {
    const SeriesInput input = series_from_body(body);
    SeriesOptions options;
    options.period = std::max(2, static_cast<int>(value_number(body, "period", 12.0)));
    options.alpha = value_number(body, "alpha", 0.05);
    options.hamed_rao = value_bool(body, "hamed_rao", true);
    options.prewhiten = value_bool(body, "prewhiten", false);
    options.seasonal_window = std::max(7, options.period + 1);
    options.trend_window = std::max(21, 2 * options.period + 1);
    options.low_pass_window = std::max(7, options.period + 1);

    double spacing = 0.0;
    if (input.time.size() > 1) spacing = (input.time.back() - input.time.front()) /
                                         static_cast<double>(input.time.size() - 1);
    if (spacing > 0.0 && spacing < 0.02)
        options.time_units_per_year = 365.0;
    else if (spacing > 0.0 && spacing < 0.5)
        options.time_units_per_year = 12.0;

    SeriesAnalysis analysis = analyze_series(input.time, input.values, options);
    analysis.label = value_string(body, "label", input.variable);

    json::Value root = json::Value::object();
    root.set("variable", json::Value::from_string(input.variable));
    root.set("units", json::Value::from_string(input.units));
    root.set("label", json::Value::from_string(analysis.label));
    root.set("samples", json::Value::from_number(static_cast<double>(analysis.n)));
    root.set("missing_filled", json::Value::from_number(static_cast<double>(analysis.missing_filled)));
    root.set("start_time", json::Value::from_number(analysis.start_time));
    root.set("end_time", json::Value::from_number(analysis.end_time));
    root.set("spacing", json::Value::from_number(analysis.spacing));
    root.set("period", json::Value::from_number(static_cast<double>(analysis.period)));
    root.set("alpha", json::Value::from_number(analysis.alpha));
    root.set("decomposed", json::Value::from_bool(analysis.decomposed));

    json::Value summary = json::Value::object();
    summary.set("mean", json::Value::from_number(analysis.summary.mean));
    summary.set("sd", json::Value::from_number(analysis.summary.sd));
    summary.set("median", json::Value::from_number(analysis.summary.median));
    summary.set("minimum", json::Value::from_number(analysis.summary.minimum));
    summary.set("maximum", json::Value::from_number(analysis.summary.maximum));
    summary.set("q25", json::Value::from_number(analysis.summary.q25));
    summary.set("q75", json::Value::from_number(analysis.summary.q75));
    root.set("summary", summary);

    json::Value decomposition = json::Value::object();
    decomposition.set("seasonal_amplitude", json::Value::from_number(analysis.seasonal_amplitude));
    decomposition.set("seasonal_strength", json::Value::from_number(analysis.seasonal_strength));
    decomposition.set("trend_strength", json::Value::from_number(analysis.trend_strength));
    decomposition.set("latest_anomaly_sigma", json::Value::from_number(analysis.latest_anomaly));
    root.set("decomposition", decomposition);

    json::Value trend = json::Value::object();
    trend.set("mann_kendall", mk_json(analysis.mk));
    trend.set("sen", sen_json(analysis.sen));
    trend.set("sen_per_decade", json::Value::from_number(analysis.sen_per_decade));
    trend.set("ensemble", ensemble_json(analysis.ensemble));
    root.set("trend", trend);

    json::Value raw = json::Value::object();
    raw.set("mann_kendall", mk_json(analysis.mk_raw));
    raw.set("sen", sen_json(analysis.sen_raw));
    raw.set("sen_per_decade", json::Value::from_number(analysis.sen_raw_per_decade));
    root.set("raw_trend", raw);

    if (analysis.period >= 2) root.set("seasonal_mann_kendall", mk_json(analysis.mk_seasonal));

    json::Value series = json::Value::object();
    series.set("time", number_array(input.time));
    series.set("observed", number_array(input.values));
    root.set("series", series);

    json::Value components = json::Value::object();
    components.set("trend", number_array(analysis.decomposition.trend));
    components.set("seasonal", number_array(analysis.decomposition.seasonal));
    components.set("remainder", number_array(analysis.decomposition.remainder));
    root.set("components", components);
    root.set("warnings", string_array(analysis.warnings));
    return json_response(root);
}

Response analyze_co2_request(const json::Value& body) {
    const SeriesInput input = series_from_body(body);
    CO2Options options;
    options.polynomial_order = std::max(1, static_cast<int>(value_number(body, "poly", 2.0)));
    options.harmonic_count = std::max(1, static_cast<int>(value_number(body, "harmonics", 4.0)));
    options.alpha = value_number(body, "alpha", 0.05);
    options.reject_iterations =
        std::max(0, static_cast<int>(value_number(body, "reject_iterations", 4.0)));
    options.reject_threshold = value_number(body, "reject_threshold", 3.0);

    double spacing = 0.0;
    if (input.time.size() > 1) spacing = (input.time.back() - input.time.front()) /
                                         static_cast<double>(input.time.size() - 1);
    if (spacing > 0.0 && spacing < 0.5) options.time_units_per_year = 12.0;
    if (spacing > 0.0 && spacing < 0.02) options.time_units_per_year = 365.0;

    CO2Analysis analysis = analyze_co2(input.time, input.values, options);
    analysis.label = value_string(body, "label", input.variable);

    json::Value root = json::Value::object();
    root.set("variable", json::Value::from_string(input.variable));
    root.set("units", json::Value::from_string(input.units));
    root.set("label", json::Value::from_string(analysis.label));
    root.set("samples", json::Value::from_number(static_cast<double>(analysis.n)));
    root.set("missing_filled", json::Value::from_number(static_cast<double>(analysis.missing_filled)));
    root.set("start_time", json::Value::from_number(analysis.start_time));
    root.set("end_time", json::Value::from_number(analysis.end_time));

    json::Value model = json::Value::object();
    model.set("polynomial_order", json::Value::from_number(analysis.model.polynomial_order));
    model.set("harmonic_count", json::Value::from_number(analysis.model.harmonic_count));
    model.set("r_squared", json::Value::from_number(analysis.model.r_squared));
    model.set("residual_sd", json::Value::from_number(analysis.model.residual_sd));
    model.set("used_points", json::Value::from_number(static_cast<double>(analysis.model.used_points)));
    model.set("rejected_points",
              json::Value::from_number(static_cast<double>(analysis.model.rejected_points)));
    root.set("model", model);

    json::Value harmonics = json::Value::array();
    for (const auto& component : analysis.components) {
        json::Value item = json::Value::object();
        item.set("number", json::Value::from_number(component.number));
        item.set("amplitude", json::Value::from_number(component.amplitude));
        item.set("phase_degrees", json::Value::from_number(component.phase_degrees));
        harmonics.push(item);
    }
    root.set("harmonics", harmonics);
    root.set("seasonal_amplitude", json::Value::from_number(analysis.seasonal_amplitude));

    json::Value growth = json::Value::object();
    growth.set("linear_coefficient", json::Value::from_number(analysis.linear_coefficient));
    growth.set("quadratic_coefficient", json::Value::from_number(analysis.quadratic_coefficient));
    growth.set("start_per_year", json::Value::from_number(analysis.growth_rate_start));
    growth.set("end_per_year", json::Value::from_number(analysis.growth_rate_end));
    growth.set("mean_per_year", json::Value::from_number(analysis.growth_rate_mean));
    growth.set("last_five_years_per_year",
               json::Value::from_number(analysis.growth_rate_last_five_years));
    growth.set("end_pgc_per_year", json::Value::from_number(analysis.growth_rate_end_pgc_per_year));
    root.set("growth", growth);

    json::Value trend = json::Value::object();
    trend.set("mann_kendall", mk_json(analysis.mk));
    trend.set("sen", sen_json(analysis.sen));
    root.set("trend", trend);

    json::Value series = json::Value::object();
    series.set("time", number_array(input.time));
    series.set("observed", number_array(input.values));
    series.set("fitted", number_array(analysis.fitted));
    series.set("residual", number_array(analysis.residual));
    root.set("series", series);
    root.set("warnings", string_array(analysis.warnings));
    return json_response(root);
}

Response sst_request(const json::Value& body) {
    NLSSTModel model;
    std::string error;
    bool loaded = read_nlsst_model("data/nlsst_model.csv", model, error);
    if (!loaded) {
        std::vector<NLSSTSample> samples;
        std::string sample_error;
        if (read_nlsst_samples("data/sst_matchups.csv", samples, sample_error) &&
            !samples.empty()) {
            model = calibrate_nlsst(samples);
            loaded = true;
        } else {
            throw BadRequest("no NLSST model: " + sample_error);
        }
    }
    if (model.coefficients.size() != 7) throw BadRequest("NLSST model is incomplete");

    const double bt11 = value_number(body, "bt11", std::numeric_limits<double>::quiet_NaN());
    const double bt12 = value_number(body, "bt12", std::numeric_limits<double>::quiet_NaN());
    const double guess =
        value_number(body, "first_guess", std::numeric_limits<double>::quiet_NaN());
    const double zenith =
        value_number(body, "zenith", std::numeric_limits<double>::quiet_NaN());
    if (!std::isfinite(bt11) || !std::isfinite(bt12) || !std::isfinite(guess) ||
        !std::isfinite(zenith))
        throw BadRequest("bt11, bt12, first_guess and zenith are required");

    const double mirror = value_number(body, "mirror", 0.0);
    const double theta_star = deg_to_rad(value_number(body, "theta_star", 0.0));
    const double result =
        nlsst_retrieve(model, bt11, bt12, guess, deg_to_rad(zenith), mirror, theta_star);

    json::Value root = json::Value::object();
    json::Value model_json = json::Value::object();
    model_json.set("coefficients", number_array(model.coefficients));
    model_json.set("rmse", json::Value::from_number(model.rmse));
    model_json.set("bias", json::Value::from_number(model.bias));
    model_json.set("samples", json::Value::from_number(static_cast<double>(model.samples)));
    model_json.set("calibrated", json::Value::from_bool(model.calibrated));
    root.set("model", model_json);
    root.set("result_k", json::Value::from_number(result));
    root.set("result_c", json::Value::from_number(result - 273.15));
    return json_response(root);
}

Response vegetation_index_request(const json::Value& body) {
    const std::string name = lower_copy(value_string(body, "index", "ndvi"));
    const double red = value_number(body, "red", std::numeric_limits<double>::quiet_NaN());
    const double nir = value_number(body, "nir", std::numeric_limits<double>::quiet_NaN());
    const double blue = value_number(body, "blue", 0.0);
    const double soil = value_number(body, "soil", 0.5);
    if (!std::isfinite(red) || !std::isfinite(nir))
        throw BadRequest("red and nir reflectances are required");

    double value = 0.0;
    if (name == "ndvi")
        value = ndvi(red, nir);
    else if (name == "evi")
        value = evi(red, nir, blue);
    else if (name == "savi")
        value = savi(red, nir, soil);
    else
        throw BadRequest("index must be ndvi, evi or savi");

    json::Value root = json::Value::object();
    root.set("index", json::Value::from_string(name));
    root.set("value", json::Value::from_number(value));
    root.set("defined", json::Value::from_bool(std::isfinite(value)));
    root.set("red", json::Value::from_number(red));
    root.set("nir", json::Value::from_number(nir));
    root.set("blue", json::Value::from_number(blue));
    root.set("soil", json::Value::from_number(soil));
    return json_response(root);
}

Response seaice_request(const json::Value& body) {
    const double tb19v = value_number(body, "tb19v", std::numeric_limits<double>::quiet_NaN());
    const double tb19h = value_number(body, "tb19h", std::numeric_limits<double>::quiet_NaN());
    const double tb37v = value_number(body, "tb37v", std::numeric_limits<double>::quiet_NaN());
    if (!std::isfinite(tb19v) || !std::isfinite(tb19h) || !std::isfinite(tb37v))
        throw BadRequest("tb19v, tb19h and tb37v are required");

    SeaIceTiePoints tie_points;
    std::string error;
    read_sea_ice_tie_points("data/nasa_team_tiepoints.csv", tie_points, error);

    const SeaIceConcentrations mixture =
        nasa_team_mixture_invert(tie_points, tb19v, tb19h, tb37v);
    const PolarGradientModel model = calibrate_polar_gradient(tie_points);
    const SeaIceConcentrations empirical =
        nasa_team_polar_gradient(model, tb19v, tb19h, tb37v);

    const auto concentrations = [](const SeaIceConcentrations& c) {
        json::Value value = json::Value::object();
        value.set("open_water", json::Value::from_number(c.open_water));
        value.set("first_year", json::Value::from_number(c.first_year));
        value.set("multiyear", json::Value::from_number(c.multiyear));
        value.set("total", json::Value::from_number(c.total()));
        value.set("extent_indicator", json::Value::from_bool(c.extent_indicator() > 0.5));
        return value;
    };

    json::Value root = json::Value::object();
    root.set("pr", json::Value::from_number(polarization_ratio(tb19v, tb19h)));
    root.set("gr", json::Value::from_number(gradient_ratio(tb37v, tb19v)));
    root.set("mixture", concentrations(mixture));
    root.set("empirical", concentrations(empirical));
    json::Value model_json = json::Value::object();
    model_json.set("rmse_cf", json::Value::from_number(model.rmse_cf));
    model_json.set("rmse_cm", json::Value::from_number(model.rmse_cm));
    root.set("model", model_json);
    return json_response(root);
}

Response precip_request(const json::Value& body) {
    const double cold = value_number(body, "cold_fraction", std::numeric_limits<double>::quiet_NaN());
    if (!std::isfinite(cold)) throw BadRequest("cold_fraction is required");
    const double water = value_number(body, "water_vapour", 4.2);
    const double a1 = value_number(body, "a1", 0.1);
    const double a2 = value_number(body, "a2", 3.0);
    const double gpi_value = gpi(cold);
    const double adjusted = gpi_with_water_vapour(gpi_value, water, a1, a2);

    json::Value root = json::Value::object();
    root.set("cold_fraction", json::Value::from_number(cold));
    root.set("water_vapour", json::Value::from_number(water));
    root.set("a1", json::Value::from_number(a1));
    root.set("a2", json::Value::from_number(a2));
    root.set("gpi", json::Value::from_number(gpi_value));
    root.set("adjusted", json::Value::from_number(adjusted));
    root.set("factor", json::Value::from_number(gpi_value > 0.0 ? adjusted / gpi_value : 0.0));
    return json_response(root);
}

Response information_request(const json::Value& body) {
    const std::vector<double> x = numeric_list(body.find("x"));
    const std::vector<double> y = numeric_list(body.find("y"));
    if (x.size() < 3 || y.size() < 3) throw BadRequest("need at least three x and y samples");
    if (x.size() != y.size()) throw BadRequest("x and y must have the same length");
    const std::size_t bins = static_cast<std::size_t>(
        std::min(64.0, std::max(2.0, value_number(body, "bins", 8.0))));
    std::size_t neighbours = static_cast<std::size_t>(
        std::max(1.0, value_number(body, "neighbours", 3.0)));
    neighbours = std::min(neighbours, x.size() - 1);

    const double plugin = mutual_information_plugin(x, y, bins);
    const double correction = miller_madow_correction(bins, bins, x.size());
    const double ksg = ksg_mutual_information(x, y, neighbours);

    json::Value root = json::Value::object();
    root.set("samples", json::Value::from_number(static_cast<double>(x.size())));
    root.set("bins", json::Value::from_number(static_cast<double>(bins)));
    root.set("neighbours", json::Value::from_number(static_cast<double>(neighbours)));
    root.set("plugin", json::Value::from_number(plugin));
    root.set("miller_madow", json::Value::from_number(correction));
    root.set("adjusted", json::Value::from_number(plugin + correction));
    root.set("ksg", json::Value::from_number(ksg));
    return json_response(root);
}

Response uncertainty_request(const json::Value& body) {
    const std::vector<double> weights = numeric_list(body.find("weights"));
    const std::vector<double> sigmas = numeric_list(body.find("sigmas"));
    if (weights.empty() || weights.size() != sigmas.size())
        throw BadRequest("weights and sigmas must be non-empty and the same length");
    const double correlation =
        std::min(0.999, std::max(0.0, value_number(body, "correlation", 0.5)));

    double weight_sum = 0.0;
    double sigma_squared = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) {
        weight_sum += weights[i];
        sigma_squared += sigmas[i] * sigmas[i];
    }
    const double sigma_mean = std::sqrt(sigma_squared / static_cast<double>(sigmas.size()));
    const double uncorrelated = uncorrelated_variance(weights, sigmas);
    const double correlated = constant_covariance_variance(sigma_mean, correlation, weights);
    const double effective = effective_sample_size(weights, correlation);

    const std::vector<double> ensemble = numeric_list(body.find("ensemble"));
    double ensemble_sd = 0.0;
    if (ensemble.size() >= 2) ensemble_sd = ensemble_standard_deviation(ensemble);
    const double total = combine_uncertainty(uncorrelated, correlated, ensemble_sd);

    json::Value root = json::Value::object();
    root.set("weight_sum", json::Value::from_number(weight_sum));
    root.set("correlation", json::Value::from_number(correlation));
    root.set("uncorrelated_variance", json::Value::from_number(uncorrelated));
    root.set("uncorrelated_sd", json::Value::from_number(std::sqrt(std::max(0.0, uncorrelated))));
    root.set("correlated_variance", json::Value::from_number(correlated));
    root.set("correlated_sd", json::Value::from_number(std::sqrt(std::max(0.0, correlated))));
    root.set("ensemble_sd", json::Value::from_number(ensemble_sd));
    root.set("total_variance", json::Value::from_number(total));
    root.set("total_sd", json::Value::from_number(std::sqrt(std::max(0.0, total))));
    root.set("effective_sample_size", json::Value::from_number(effective));
    root.set("correlation_sensitivity",
             json::Value::from_number(correlation_sensitivity(sigma_mean, weights)));
    json::Value sensitivities = json::Value::array();
    for (std::size_t i = 0; i < weights.size(); ++i)
        sensitivities.push(json::Value::from_number(
            uncorrelated_weight_sensitivity(weights[i], sigmas[i])));
    root.set("weight_sensitivities", sensitivities);
    return json_response(root);
}

Response physics_request(const json::Value& body) {
    const std::string formula = value_string(body, "formula");
    json::Value root = json::Value::object();
    json::Value primary = json::Value::object();
    json::Value details = json::Value::array();
    const auto detail = [&details](const std::string& label, double value,
                                   const std::string& unit) {
        json::Value item = json::Value::object();
        item.set("label", json::Value::from_string(label));
        item.set("value", json::Value::from_number(value));
        item.set("unit", json::Value::from_string(unit));
        details.push(item);
    };

    if (formula == "toa") {
        const double forcing = value_number(body, "F", 3.7);
        const double feedback = value_number(body, "lambda", -1.2);
        const double delta = value_number(body, "delta_t", 1.2);
        const double value = toa_energy_imbalance(forcing, feedback, delta);
        root.set("title", json::Value::from_string("TOA energy imbalance"));
        root.set("formula", json::Value::from_string("N = F + lambda dT"));
        primary.set("value", json::Value::from_number(value));
        primary.set("unit", json::Value::from_string("W/m2"));
        detail("Forcing", forcing, "W/m2");
        detail("Feedback", feedback, "W/m2/K");
        detail("Warming", delta, "K");
    } else if (formula == "forcing-co2") {
        const double ppm = value_number(body, "ppm", 420.0);
        const double reference = value_number(body, "reference_ppm", 280.0);
        const double value = co2_radiative_forcing(ppm, reference);
        root.set("title", json::Value::from_string("CO2 radiative forcing"));
        root.set("formula", json::Value::from_string("dF = 5.35 ln(C/C0)"));
        primary.set("value", json::Value::from_number(value));
        primary.set("unit", json::Value::from_string("W/m2"));
        detail("Concentration", ppm, "ppm");
        detail("Reference", reference, "ppm");
        detail("Ratio", ppm / reference, "");
    } else if (formula == "ecs") {
        const double feedback = value_number(body, "lambda", -1.25);
        const double forcing = value_number(body, "forcing_2x", 3.7);
        const double ecs = climate_sensitivity_2x(feedback, forcing);
        root.set("title", json::Value::from_string("Equilibrium climate sensitivity"));
        root.set("formula", json::Value::from_string("ECS = -F2x / lambda"));
        primary.set("value", json::Value::from_number(ecs));
        primary.set("unit", json::Value::from_string("K"));
        detail("Feedback", feedback, "W/m2/K");
        detail("2xCO2 forcing", forcing, "W/m2");
        detail("Equilibrium warming", equilibrium_temperature_change(forcing, feedback), "K");
    } else if (formula == "mie") {
        const double size = value_number(body, "size", 5.0);
        const double real = value_number(body, "real", 1.33);
        const double imag = value_number(body, "imag", 0.01);
        const MieResult mie = mie_sphere(size, std::complex<double>(real, imag));
        root.set("title", json::Value::from_string("Mie scattering"));
        root.set("formula", json::Value::from_string("homogeneous sphere, size parameter x"));
        primary.set("value", json::Value::from_number(mie.qext));
        primary.set("unit", json::Value::from_string("Qext"));
        detail("Qsca", mie.qsca, "");
        detail("Asymmetry g", mie.g, "");
        detail("Single-scattering albedo", mie.single_scattering_albedo, "");
        detail("Qabs", mie.qext - mie.qsca, "");
    } else if (formula == "wind") {
        const double u = value_number(body, "u", 10.0);
        const double v = value_number(body, "v", 0.0);
        const double rho = value_number(body, "rho", 1.225);
        const double cd = value_number(body, "cd", 1.3e-3);
        const WindStress stress = bulk_wind_stress(rho, cd, u, v);
        root.set("title", json::Value::from_string("Bulk wind stress"));
        root.set("formula", json::Value::from_string("tau = rho_a Cd |U| U"));
        primary.set("value", json::Value::from_number(stress.magnitude()));
        primary.set("unit", json::Value::from_string("N/m2"));
        detail("tau_x", stress.tau_x, "N/m2");
        detail("tau_y", stress.tau_y, "N/m2");
        detail("Wind speed", std::sqrt(u * u + v * v), "m/s");
    } else if (formula == "ekman") {
        const double tau_x = value_number(body, "tau_x", 0.1);
        const double tau_y = value_number(body, "tau_y", 0.0);
        const double latitude = value_number(body, "latitude", 30.0);
        const double rho = value_number(body, "rho", 1025.0);
        const double f = coriolis_parameter(latitude);
        const EkmanTransport transport = ekman_transport(tau_x, tau_y, rho, f);
        root.set("title", json::Value::from_string("Ekman transport"));
        root.set("formula", json::Value::from_string("M = (tau_y, -tau_x)/(rho f)"));
        primary.set("value", json::Value::from_number(transport.magnitude()));
        primary.set("unit", json::Value::from_string("m2/s"));
        detail("M_x", transport.m_x, "m2/s");
        detail("M_y", transport.m_y, "m2/s");
        detail("Coriolis f", f, "1/s");
    } else if (formula == "aod") {
        const double n = value_number(body, "n", 1e8);
        const double radius_um = value_number(body, "radius", 0.1);
        const double qext = value_number(body, "qext", 2.0);
        const double path = value_number(body, "path", 1000.0);
        const double radius = radius_um * 1e-6;
        const double value = aerosol_optical_depth(n, radius, qext, path);
        root.set("title", json::Value::from_string("Aerosol optical depth"));
        root.set("formula", json::Value::from_string("tau = N pi r^2 Qext L"));
        primary.set("value", json::Value::from_number(value));
        primary.set("unit", json::Value::from_string("AOD"));
        detail("Cross section", PI * radius * radius * qext, "m2");
        detail("Column density", n * path, "m-2");
    } else if (formula == "equilibrium") {
        const double forcing = value_number(body, "F", 3.7);
        const double feedback = value_number(body, "lambda", -1.25);
        root.set("title", json::Value::from_string("Equilibrium temperature change"));
        root.set("formula", json::Value::from_string("dT = -F / lambda"));
        primary.set("value", json::Value::from_number(equilibrium_temperature_change(forcing, feedback)));
        primary.set("unit", json::Value::from_string("K"));
        detail("Forcing", forcing, "W/m2");
        detail("Feedback", feedback, "W/m2/K");
    } else if (formula == "coriolis") {
        const double latitude = value_number(body, "latitude", 30.0);
        const double f = coriolis_parameter(latitude);
        root.set("title", json::Value::from_string("Coriolis parameter"));
        root.set("formula", json::Value::from_string("f = 2 Omega sin(phi)"));
        primary.set("value", json::Value::from_number(f));
        primary.set("unit", json::Value::from_string("1/s"));
        detail("Latitude", latitude, "deg");
        detail("Inertial period", f != 0.0 ? 2.0 * PI / std::fabs(f) : 0.0, "s");
    } else if (formula == "ekman-pumping") {
        const double dtau_y_dx = value_number(body, "dtau_y_dx", 1e-8);
        const double dtau_x_dy = value_number(body, "dtau_x_dy", 0.0);
        const double rho = value_number(body, "rho", 1025.0);
        const double f = coriolis_parameter(value_number(body, "latitude", 30.0));
        const double value = ekman_pumping(dtau_y_dx, dtau_x_dy, rho, f);
        root.set("title", json::Value::from_string("Ekman pumping"));
        root.set("formula", json::Value::from_string("w = curl(tau)/(rho f)"));
        primary.set("value", json::Value::from_number(value));
        primary.set("unit", json::Value::from_string("m/s"));
        detail("Coriolis f", f, "1/s");
        detail("Stress curl", dtau_y_dx - dtau_x_dy, "N/m2/m");
    } else if (formula == "geostrophic") {
        const double rho = value_number(body, "rho", 1025.0);
        const double dp_dx = value_number(body, "dp_dx", 1e-5);
        const double dp_dy = value_number(body, "dp_dy", 0.0);
        const double latitude = value_number(body, "latitude", 30.0);
        const GeostrophicVelocity velocity = geostrophic_velocity(rho, dp_dx, dp_dy, latitude);
        root.set("title", json::Value::from_string("Geostrophic velocity"));
        root.set("formula", json::Value::from_string("u = -dp/dy/(rho f), v = dp/dx/(rho f)"));
        primary.set("value", json::Value::from_number(std::sqrt(velocity.u * velocity.u +
                                                                 velocity.v * velocity.v)));
        primary.set("unit", json::Value::from_string("m/s"));
        detail("u", velocity.u, "m/s");
        detail("v", velocity.v, "m/s");
        detail("Coriolis f", coriolis_parameter(latitude), "1/s");
    } else if (formula == "thermal-wind") {
        const double rho = value_number(body, "rho", 1025.0);
        const double drho_dx = value_number(body, "drho_dx", 1e-5);
        const double drho_dy = value_number(body, "drho_dy", 0.0);
        const double latitude = value_number(body, "latitude", 30.0);
        const ThermalWind shear = thermal_wind(rho, drho_dx, drho_dy, latitude);
        root.set("title", json::Value::from_string("Thermal wind shear"));
        root.set("formula", json::Value::from_string("du/dz = g drho/dy/(f rho)"));
        primary.set("value", json::Value::from_number(
                                  std::sqrt(shear.du_dz * shear.du_dz + shear.dv_dz * shear.dv_dz)));
        primary.set("unit", json::Value::from_string("1/s"));
        detail("du/dz", shear.du_dz, "1/s");
        detail("dv/dz", shear.dv_dz, "1/s");
    } else if (formula == "savi") {
        const double red = value_number(body, "red", 0.08);
        const double nir = value_number(body, "nir", 0.45);
        const double soil = value_number(body, "soil", 0.5);
        const double value = savi(red, nir, soil);
        root.set("title", json::Value::from_string("Soil-adjusted vegetation index"));
        root.set("formula", json::Value::from_string("SAVI = (NIR-RED)/(NIR+RED+L)(1+L)"));
        primary.set("value", json::Value::from_number(value));
        primary.set("unit", json::Value::from_string("index"));
        detail("RED", red, "reflectance");
        detail("NIR", nir, "reflectance");
        detail("L", soil, "");
    } else if (formula == "evi") {
        const double red = value_number(body, "red", 0.08);
        const double nir = value_number(body, "nir", 0.45);
        const double blue = value_number(body, "blue", 0.04);
        const double gain = value_number(body, "gain", 2.5);
        const double c1 = value_number(body, "c1", 6.0);
        const double c2 = value_number(body, "c2", 7.5);
        const double canopy = value_number(body, "canopy", 1.0);
        const double value = evi(red, nir, blue, gain, c1, c2, canopy);
        root.set("title", json::Value::from_string("Enhanced vegetation index"));
        root.set("formula", json::Value::from_string("EVI = G(NIR-RED)/(NIR+C1 RED-C2 BLUE+L)"));
        primary.set("value", json::Value::from_number(value));
        primary.set("unit", json::Value::from_string("index"));
        detail("RED", red, "reflectance");
        detail("NIR", nir, "reflectance");
        detail("BLUE", blue, "reflectance");
    } else {
        throw BadRequest("unknown physics formula");
    }

    root.set("primary", primary);
    root.set("details", details);
    root.set("note", json::Value::from_string("Values use the same C++ functions as the command line."));
    return json_response(root);
}

// What kind of change does this variable measure? The label is what the
// challenge asks a judge to compare, so it is stated per variable instead of
// hiding behind a retrieval formula.
std::string change_label(const std::string& id) {
    static const std::map<std::string, std::string> labels = {
        {"sst", "Trend in mean ocean skin temperature (warming / cooling)"},
        {"air_temperature", "Trend in mean surface air temperature"},
        {"co2", "Trend and growth-rate change in atmospheric CO2"},
        {"ch4", "Trend in atmospheric methane concentration"},
        {"ndvi", "Trend in vegetation greenness and growing season"},
        {"sea_ice_extent", "Trend in sea-ice extent (shrinking / expanding)"},
        {"sea_ice_concentration", "Trend in ice concentration (thickening / thinning)"},
        {"precipitation", "Trend in precipitation totals and extremes"},
        {"ocean_heat_content", "Trend in ocean heat uptake"},
        {"toa_energy_imbalance", "Trend in the top-of-atmosphere energy imbalance"},
        {"snow_cover", "Trend in snow-cover extent and duration"},
        {"ozone", "Trend in column ozone (recovery / depletion)"},
        {"aerosol_optical_depth", "Trend in aerosol loading"},
        {"soil_moisture", "Trend in soil moisture (wetting / drying)"},
        {"ocean_color", "Trend in ocean chlorophyll and productivity"},
        {"sea_level", "Trend in sea level (rise / fall)"},
    };
    const auto found = labels.find(id);
    return found == labels.end() ? "Trend in the measured variable" : found->second;
}

Response catalog_request() {
    json::Value list = json::Value::array();
    for (const auto& entry : catalog_entries()) {
        json::Value item = json::Value::object();
        item.set("id", json::Value::from_string(entry.id));
        item.set("name", json::Value::from_string(entry.name));
        item.set("units", json::Value::from_string(entry.units));
        item.set("processes", string_array(entry.processes));
        item.set("regions", string_array(entry.regions));
        item.set("retrieval", json::Value::from_string(entry.retrieval));
        item.set("analysis", json::Value::from_string(entry.analysis_method));
        item.set("description", json::Value::from_string(entry.description));
        item.set("change", json::Value::from_string(change_label(entry.id)));
        list.push(item);
    }
    return json_response(list);
}

Response overview_request() {
    json::Value root = json::Value::object();
    root.set("catalog_count", json::Value::from_number(
                                 static_cast<double>(catalog_entries().size())));
    root.set("processes", string_array(catalog_processes()));
    root.set("regions", string_array(catalog_regions()));

    std::vector<std::string> files;
    try {
        for (const auto& entry : std::filesystem::directory_iterator("data")) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() == ".csv")
                files.push_back(entry.path().filename().string());
        }
    } catch (const std::exception&) {
    }
    std::sort(files.begin(), files.end());
    root.set("files", string_array(files));

    struct Module {
        const char* code;
        const char* name;
        const char* description;
        const char* target;
    };
    const Module modules[] = {
        {"TS", "Time series", "STL decomposition, Mann–Kendall, Sen's slope and Hamed–Rao.", "timeseries"},
        {"CO", "CO2 curve fit", "Polynomial plus harmonics with iterative outlier rejection.", "co2"},
        {"RS", "Retrievals", "NLSST, NDVI/EVI/SAVI, NASA Team sea ice and GPI precipitation.", "retrieval"},
        {"MI", "Information theory", "Histogram MI, Miller–Madow and the KSG k-NN estimator.", "info"},
        {"S", "Uncertainty", "Uncorrelated, correlated and ensemble error components.", "uncertainty"},
        {"PH", "Physics", "Energy balance, forcing, Mie scattering and wind-driven transport.", "physics"},
        {"CA", "Catalog", "Process- and region-tagged climate variables.", "catalog"},
        {"IN", "Intent compiler", "Deterministic natural-language question to variable mapping.", "overview"},
    };
    json::Value module_list = json::Value::array();
    for (const auto& module : modules) {
        json::Value item = json::Value::object();
        item.set("code", json::Value::from_string(module.code));
        item.set("name", json::Value::from_string(module.name));
        item.set("desc", json::Value::from_string(module.description));
        item.set("target", json::Value::from_string(module.target));
        module_list.push(item);
    }
    root.set("modules", module_list);
    return json_response(root);
}

std::vector<std::string> split_command_line(const std::string& text) {
    std::vector<std::string> tokens;
    std::string current;
    char quote = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (quote != 0) {
            if (c == quote) {
                quote = 0;
            } else if (c == '\\' && i + 1 < text.size() && text[i + 1] == quote) {
                current.push_back(text[++i]);
            } else {
                current.push_back(c);
            }
        } else if (c == '"' || c == '\'') {
            quote = c;
        } else if (std::isspace(static_cast<unsigned char>(c))) {
            if (!current.empty()) {
                tokens.push_back(current);
                current.clear();
            }
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty()) tokens.push_back(current);
    return tokens;
}

Response help_request() {
    const CommandResult result = run_command("./climate help");
    json::Value root = json::Value::object();
    root.set("exit_code", json::Value::from_number(result.exit_code));
    root.set("text", json::Value::from_string(result.output));
    return json_response(root);
}

Response run_request(const json::Value& body) {
    const std::string command = value_string(body, "command");
    std::vector<std::string> tokens = split_command_line(command);
    if (tokens.empty()) throw BadRequest("enter a climate command");
    if (tokens[0] == "./climate" || tokens[0] == "climate") tokens.erase(tokens.begin());
    if (tokens.empty()) throw BadRequest("enter a command after 'climate'");
    if (tokens[0] == "ui") throw BadRequest("the ui command cannot be started from the console");

    std::vector<std::string> arguments = {"./climate"};
    arguments.insert(arguments.end(), tokens.begin(), tokens.end());
    const CommandResult result = run_command(build_command(arguments));
    std::string output = result.output;
    if (output.size() > 200000) output = output.substr(0, 200000) + "\n[output truncated]";

    json::Value root = json::Value::object();
    root.set("command", json::Value::from_string(build_command(tokens)));
    root.set("exit_code", json::Value::from_number(result.exit_code));
    root.set("output", json::Value::from_string(output));
    return json_response(root);
}

std::string collapse_whitespace(const std::string& text) {
    std::string out;
    bool space = false;
    for (char c : text) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            space = true;
            continue;
        }
        if (space && !out.empty()) out.push_back(' ');
        space = false;
        out.push_back(c);
    }
    return trim_copy(out);
}

Response formulas_request() {
    std::vector<std::filesystem::path> headers;
    try {
        for (const auto& entry : std::filesystem::directory_iterator("src")) {
            if (entry.is_regular_file() && entry.path().extension() == ".hpp")
                headers.push_back(entry.path().filename());
        }
    } catch (const std::exception&) {
    }
    std::sort(headers.begin(), headers.end());

    json::Value modules = json::Value::array();
    std::size_t total = 0;
    for (const auto& header : headers) {
        const std::string path = "src/" + header.string();
        std::string source = read_text_file(path);
        std::string cleaned;
        cleaned.reserve(source.size());
        bool line_comment = false;
        bool block_comment = false;
        for (std::size_t i = 0; i < source.size(); ++i) {
            const char c = source[i];
            const char next = i + 1 < source.size() ? source[i + 1] : '\0';
            if (line_comment) {
                if (c == '\n') {
                    line_comment = false;
                    cleaned.push_back('\n');
                }
                continue;
            }
            if (block_comment) {
                if (c == '*' && next == '/') {
                    block_comment = false;
                    ++i;
                }
                continue;
            }
            if (c == '/' && next == '/') {
                line_comment = true;
                ++i;
                continue;
            }
            if (c == '/' && next == '*') {
                block_comment = true;
                ++i;
                continue;
            }
            cleaned.push_back(c);
        }

        json::Value functions = json::Value::array();
        std::size_t begin = 0;
        while (begin <= cleaned.size()) {
            const std::size_t end = cleaned.find(';', begin);
            const std::string statement = collapse_whitespace(
                cleaned.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
            begin = end == std::string::npos ? cleaned.size() + 1 : end + 1;
            if (statement.empty() || statement.find('(') == std::string::npos) continue;
            const std::string lowered = lower_copy(statement);
            static const char* forbidden[] = {"#", "using ", "typedef ", "namespace ",
                                              "struct ", "class ", "enum ", "template ",
                                              "friend ", "return ", "if ", "for ",
                                              "while ", "switch ", "operator"};
            bool skip = false;
            for (const char* prefix : forbidden) {
                if (lowered.rfind(prefix, 0) == 0) {
                    skip = true;
                    break;
                }
            }
            if (skip || statement.find('{') != std::string::npos) continue;
            const std::size_t paren = statement.find('(');
            std::size_t name_end = paren;
            while (name_end > 0 && std::isspace(static_cast<unsigned char>(statement[name_end - 1])))
                --name_end;
            std::size_t name_begin = name_end;
            while (name_begin > 0 &&
                   (std::isalnum(static_cast<unsigned char>(statement[name_begin - 1])) ||
                    statement[name_begin - 1] == '_'))
                --name_begin;
            if (name_begin == name_end) continue;
            functions.push(json::Value::from_string(statement));
            ++total;
        }

        json::Value module = json::Value::object();
        module.set("file", json::Value::from_string(header.string()));
        module.set("count", json::Value::from_number(static_cast<double>(functions.items.size())));
        module.set("functions", functions);
        modules.push(module);
    }

    json::Value root = json::Value::object();
    root.set("count", json::Value::from_number(static_cast<double>(total)));
    root.set("modules", modules);
    return json_response(root);
}

Response lab_catalog_request() {
    const std::vector<LabFormula>& formulas = lab_formulas();
    std::vector<std::string> order;
    std::map<std::string, json::Value> groups;
    std::map<std::string, std::string> labels;
    std::map<std::string, std::size_t> counts;
    for (const LabFormula& formula : formulas) {
        if (groups.find(formula.stage) == groups.end()) {
            order.push_back(formula.stage);
            groups[formula.stage] = json::Value::array();
            labels[formula.stage] = formula.stage_label;
            counts[formula.stage] = 0;
        }
        json::Value item = json::Value::object();
        item.set("id", json::Value::from_string(formula.id));
        item.set("name", json::Value::from_string(formula.name));
        item.set("module", json::Value::from_string(formula.module));
        item.set("returns", json::Value::from_string(formula.returns));
        item.set("signature", json::Value::from_string(formula.signature));
        json::Value params = json::Value::array();
        for (const LabParam& param : formula.params) {
            json::Value entry = json::Value::object();
            entry.set("name", json::Value::from_string(param.name));
            entry.set("type", json::Value::from_string(param.type));
            entry.set("kind", json::Value::from_string(param.kind));
            entry.set("default", json::Value::from_number(param.default_value));
            params.push(entry);
        }
        item.set("params", params);
        groups[formula.stage].push(item);
        counts[formula.stage] += 1;
    }

    std::sort(order.begin(), order.end());
    json::Value stages = json::Value::array();
    for (const std::string& stage : order) {
        json::Value group = json::Value::object();
        group.set("id", json::Value::from_string(stage));
        group.set("label", json::Value::from_string(labels[stage]));
        group.set("count", json::Value::from_number(static_cast<double>(counts[stage])));
        group.set("formulas", groups[stage]);
        stages.push(group);
    }

    json::Value root = json::Value::object();
    root.set("count", json::Value::from_number(static_cast<double>(formulas.size())));
    root.set("groups", stages);
    return json_response(root);
}

Response lab_run_request(const json::Value& body) {
    std::string file = value_string(body, "file", "sst_monthly.csv");
    if (file.empty() || file.find('/') != std::string::npos ||
        file.find('\\') != std::string::npos || file.find("..") != std::string::npos)
        throw BadRequest("invalid data file name");
    const Table table = parse_csv_text(read_text_file("data/" + file));
    if (table.rows.size() < 3) throw BadRequest("the series file has too few rows");

    const std::string value_name = value_string(body, "value");
    int value_column = value_name.empty() ? table.column("value") : table.column(value_name);
    if (!value_name.empty() && value_column < 0)
        throw BadRequest("value column '" + value_name + "' was not found");
    if (value_column < 0 && table.header.size() > 1) value_column = 1;
    if (value_column < 0) value_column = 0;

    LabInputs inputs;
    const std::size_t columns = static_cast<std::size_t>(
        std::max(1.0, value_number(body, "matrix_columns", 3.0)));
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        double time = static_cast<double>(i);
        parse_number(table.cell(i, 0), time);
        double value = 0.0;
        if (!parse_number(table.cell(i, value_column), value)) continue;
        inputs.time.push_back(time);
        inputs.series.push_back(value);
    }
    if (inputs.series.size() < 3) throw BadRequest("the series file has too few numeric rows");
    for (std::size_t i = 0; i + columns <= inputs.series.size(); i += columns) {
        inputs.matrix.emplace_back(inputs.series.begin() + static_cast<long>(i),
                                   inputs.series.begin() + static_cast<long>(i + columns));
    }
    const std::size_t side = std::min<std::size_t>(
        32, std::max<std::size_t>(2, static_cast<std::size_t>(
                std::floor(std::sqrt(static_cast<double>(inputs.series.size()))))));
    for (std::size_t r = 0; r < side; ++r) {
        std::vector<double> row;
        row.reserve(side);
        for (std::size_t c = 0; c < side; ++c) {
            const std::size_t index = r * side + c;
            row.push_back(index < inputs.series.size() ? inputs.series[index] : 0.0);
        }
        inputs.square.push_back(row);
    }
    inputs.text = value_string(body, "text");
    if (inputs.text.empty()) {
        // Formulas that take a path should not drop files into the repository.
        const std::filesystem::path scratch =
            std::filesystem::temp_directory_path() / "climate-formula-lab";
        std::error_code ignored;
        std::filesystem::create_directories(scratch, ignored);
        inputs.text = (scratch / "out.csv").string();
    }

    static const std::map<std::string, std::size_t> lookup = [] {
        std::map<std::string, std::size_t> index;
        for (const LabFormula& formula : lab_formulas()) index[formula.id] = formula.index;
        return index;
    }();

    const json::Value* calls = body.find("calls");
    if (calls == nullptr || !calls->is_array() || calls->items.empty())
        throw BadRequest("select at least one formula to run");

    json::Value results = json::Value::array();
    std::size_t succeeded = 0;
    for (const json::Value& call : calls->items) {
        if (!call.is_object()) continue;
        const std::string id = value_string(call, "id");
        const auto found = lookup.find(id);
        json::Value item = json::Value::object();
        item.set("id", json::Value::from_string(id));
        if (found == lookup.end()) {
            item.set("ok", json::Value::from_bool(false));
            item.set("error", json::Value::from_string("unknown formula id"));
            results.push(item);
            continue;
        }
        const LabFormula& formula = lab_formulas()[found->second];
        item.set("name", json::Value::from_string(formula.name));
        item.set("module", json::Value::from_string(formula.module));
        item.set("signature", json::Value::from_string(formula.signature));

        LabInputs call_inputs = inputs;
        const json::Value* args = call.find("args");
        for (const LabParam& param : formula.params) {
            if (param.kind != "number" && param.kind != "complex" && param.kind != "seed" &&
                param.kind != "enum")
                continue;
            double value = param.default_value;
            if (args != nullptr && args->is_object()) {
                const json::Value* supplied = args->find(param.name);
                if (supplied != nullptr && supplied->is_number()) value = supplied->number;
            }
            call_inputs.numbers.push_back(value);
        }

        try {
            item.set("ok", json::Value::from_bool(true));
            item.set("result", json::Value::from_string(run_formula(formula.index, call_inputs)));
            ++succeeded;
        } catch (const std::exception& error) {
            item.set("ok", json::Value::from_bool(false));
            item.set("error", json::Value::from_string(error.what()));
        }
        results.push(item);
    }

    json::Value root = json::Value::object();
    root.set("series", json::Value::from_string(file));
    root.set("samples", json::Value::from_number(static_cast<double>(inputs.series.size())));
    root.set("requested", json::Value::from_number(static_cast<double>(results.items.size())));
    root.set("succeeded", json::Value::from_number(static_cast<double>(succeeded)));
    root.set("results", results);
    return json_response(root);
}

std::string live_interpreter() {
    if (std::filesystem::exists(".venv-test/bin/python")) return ".venv-test/bin/python";
    return "python3";
}

bool safe_token(const std::string& text, const std::string& extra) {
    if (text.empty() || text.size() > 40) return false;
    for (char c : text) {
        const bool ok = std::isalnum(static_cast<unsigned char>(c)) ||
                        extra.find(c) != std::string::npos;
        if (!ok) return false;
    }
    return true;
}

Response live_status_request() {
    const CommandResult result = run_command(live_interpreter() + " tools/nasa_live.py status < /dev/null");
    if (result.exit_code != 0) throw BadRequest("the live-data helper failed: " + result.output);
    return json_response(json::parse(result.output));
}

Response live_pull_request(const json::Value& body) {
    const std::string source = value_string(body, "source");
    static const char* allowed[] = {"gistemp", "power", "seaice", "earthaccess"};
    if (std::find(std::begin(allowed), std::end(allowed), source) == std::end(allowed))
        throw BadRequest("unknown live source");

    std::ostringstream command;
    command << live_interpreter() << " tools/nasa_live.py pull --source " << source;
    const std::string start = value_string(body, "start", "2000-01-01");
    const std::string end = value_string(body, "end");
    if (!safe_token(start, "-") || (!end.empty() && !safe_token(end, "-")))
        throw BadRequest("start and end must be ISO dates");
    if (source == "power" || source == "earthaccess") {
        command << " --start " << start;
        if (!end.empty()) command << " --end " << end;
    }
    if (source == "power") {
        const double lat = value_number(body, "lat", 72.0);
        const double lon = value_number(body, "lon", -40.0);
        if (lat < -90.0 || lat > 90.0 || lon < -180.0 || lon > 180.0)
            throw BadRequest("latitude or longitude is out of range");
        const std::string parameter = value_string(body, "parameter", "T2M");
        if (!safe_token(parameter, "_:")) throw BadRequest("invalid POWER parameter");
        command << " --lat " << lat << " --lon " << lon << " --parameter " << parameter;
    }
    if (source == "seaice") {
        const std::string region = value_string(body, "region", "N");
        if (region != "N" && region != "S") throw BadRequest("region must be N or S");
        command << " --region " << region;
    }
    if (source == "earthaccess") {
        const std::string short_name = value_string(body, "short_name");
        if (!safe_token(short_name, "_-.")) throw BadRequest("invalid CMR short name");
        command << " --short-name " << short_name;
    }

    // Per-request Earthdata credentials: they are handed to the helper through
    // its environment only, never written to disk, logged or echoed back.
    std::vector<std::pair<std::string, std::string>> env;
    auto credential = [&](const char* key, const char* name) {
        const std::string value = value_string(body, key);
        if (value.empty()) return;
        if (value.size() > 400) throw BadRequest(std::string(name) + " is implausibly long");
        env.emplace_back(name, value);
    };
    credential("earthdata_username", "EARTHDATA_USERNAME");
    credential("earthdata_password", "EARTHDATA_PASSWORD");
    credential("earthdata_token", "EARTHDATA_TOKEN");

    const CommandResult result = env.empty()
        ? run_command(command.str() + " < /dev/null")
        : run_command_with_env(command.str() + " < /dev/null", env);
    // A library can print warnings on stderr before its JSON verdict, so the
    // payload is taken from the last line that actually is JSON.
    std::string payload = result.output;
    const std::size_t marker = payload.rfind("{\n  \"ok\"");
    const std::size_t inline_marker = payload.rfind("{\"ok\"");
    if (inline_marker != std::string::npos) payload = payload.substr(inline_marker);
    else if (marker != std::string::npos) payload = payload.substr(marker);
    try {
        return json_response(json::parse(payload));
    } catch (const std::exception&) {
        throw BadRequest("the live-data helper returned something unreadable: " +
                         result.output.substr(result.output.size() > 400 ? result.output.size() - 400 : 0));
    }
}

// ---------------------------------------------------------------------------
// Regional contrast: the same variable, two places, two independent verdicts.
// A measured variable can rise in one region and fall in another even when one
// process drives both; this endpoint makes that opposition measurable.
// ---------------------------------------------------------------------------
struct ContrastRegion {
    std::string label;
    double lat = 0.0;
    double lon = 0.0;
    std::string parameter;
    std::string units;
    std::string file;
    std::size_t samples = 0;
    double start = 0.0;
    double end = 0.0;
    double first_value = 0.0;
    double last_value = 0.0;
    double mean = 0.0;
    double slope_per_year = 0.0;
    double change = 0.0;         // reported as +rise / -fall
    double percent = 0.0;
    double p_value = 1.0;
    std::string direction = "no trend";
    bool significant = false;
};

json::Value contrast_point_json(const ContrastRegion& point) {
    json::Value value = json::Value::object();
    value.set("label", json::Value::from_string(point.label));
    value.set("lat", json::Value::from_number(point.lat));
    value.set("lon", json::Value::from_number(point.lon));
    value.set("parameter", json::Value::from_string(point.parameter));
    value.set("units", json::Value::from_string(point.units));
    value.set("file", json::Value::from_string(point.file));
    value.set("samples", json::Value::from_number(static_cast<double>(point.samples)));
    value.set("start", json::Value::from_number(point.start));
    value.set("end", json::Value::from_number(point.end));
    value.set("first_value", json::Value::from_number(point.first_value));
    value.set("last_value", json::Value::from_number(point.last_value));
    value.set("mean", json::Value::from_number(point.mean));
    value.set("slope_per_year", json::Value::from_number(point.slope_per_year));
    value.set("change", json::Value::from_number(point.change));
    value.set("percent", json::Value::from_number(point.percent));
    value.set("p_value", json::Value::from_number(point.p_value));
    value.set("direction", json::Value::from_string(point.direction));
    value.set("significant", json::Value::from_bool(point.significant));
    return value;
}

ContrastRegion contrast_region(const std::string& label, double lat, double lon,
                               const std::string& parameter, const std::string& start,
                               const std::string& end) {
    std::ostringstream command;
    command << live_interpreter() << " tools/nasa_live.py pull --source power"
            << " --lat " << lat << " --lon " << lon << " --parameter " << parameter
            << " --start " << start << " --end " << end << " < /dev/null";
    const CommandResult result = run_command(command.str());
    json::Value payload;
    try {
        payload = json::parse(result.output);
    } catch (const std::exception&) {
        throw BadRequest("the NASA POWER pull for " + label +
                         " returned something unreadable: " + result.output.substr(0, 300));
    }
    const json::Value* ok = payload.find("ok");
    if (ok == nullptr || !ok->as_bool(false)) {
        const json::Value* error = payload.find("error");
        throw BadRequest("the NASA POWER pull for " + label + " failed: " +
                         (error == nullptr ? std::string("no detail") : error->as_string("")));
    }
    const json::Value* series = payload.find("series");
    if (series == nullptr) throw BadRequest("the NASA POWER pull returned no series for " + label);

    const json::Value* file = series->find("file");
    if (file == nullptr) throw BadRequest("the NASA POWER pull returned no file for " + label);
    const std::string path = file->as_string("");
    if (path.size() < 4 || path.substr(path.size() - 4) != ".csv")
        throw BadRequest("unexpected series path from the live helper");
    if (path.find("..") != std::string::npos) throw BadRequest("unsafe series path");

    const Table table = parse_csv_text(read_text_file(path));
    std::vector<double> times;
    std::vector<double> values;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        double time = 0.0;
        double value = 0.0;
        if (!parse_number(table.cell(i, 0), time)) continue;
        if (!parse_number(table.cell(i, 1), value)) continue;
        times.push_back(time);
        values.push_back(value);
    }
    if (values.size() < 12)
        throw BadRequest("only " + std::to_string(values.size()) +
                         " usable samples were returned for " + label);

    ContrastRegion region;
    region.label = label;
    region.lat = lat;
    region.lon = lon;
    region.parameter = parameter;
    region.file = path;
    region.samples = values.size();
    region.start = times.front();
    region.end = times.back();
    region.first_value = values.front();
    region.last_value = values.back();

    double total = 0.0;
    for (double value : values) total += value;
    region.mean = total / static_cast<double>(values.size());

    const SenResult sen = sens_slope(times, values);
    // Monthly series carry a seasonal cycle, so the seasonal Mann-Kendall test
    // compares months against the same month in other years instead of pooling
    // winter against summer.
    const double spacing = (times.back() - times.front()) /
                           static_cast<double>(values.size() - 1);
    const MKResult mk = (values.size() >= 36 && spacing > 0.06 && spacing < 0.10)
                            ? seasonal_mann_kendall(values, 12)
                            : mann_kendall(values);
    region.slope_per_year = sen.slope;
    region.p_value = mk.p_value;
    region.direction = mk.direction;
    region.significant = mk.significant;
    const double span = region.end - region.start;
    region.change = sen.slope * span;
    region.percent = std::fabs(region.mean) > 1e-9
                         ? 100.0 * region.change / std::fabs(region.mean)
                         : 0.0;

    const json::Value* units = series->find("units");
    region.units = units == nullptr ? std::string() : units->as_string("");
    return region;
}

Response contrast_request(const json::Value& body) {
    const std::string parameter = value_string(body, "parameter", "T2M");
    static const char* allowed[] = {"T2M",   "T2M_MAX", "T2M_MIN", "TS", "PRECTOTCORR",
                                    "WS2M", "RH2M",    "ALLSKY_SFC_SW_DWN"};
    if (std::find(std::begin(allowed), std::end(allowed), parameter) == std::end(allowed))
        throw BadRequest("unsupported POWER parameter for the contrast");

    const int start_year = static_cast<int>(value_number(body, "start", 1990.0));
    const int end_year = static_cast<int>(value_number(body, "end", 2024.0));
    if (start_year < 1981 || end_year > 2030 || end_year - start_year < 10)
        throw BadRequest("the contrast needs at least 10 years inside 1981-2030");
    const std::string start = std::to_string(start_year) + "-01-01";
    const std::string end = std::to_string(end_year) + "-12-31";

    struct PointSpec {
        double lat;
        double lon;
        const char* fallback_label;
    };
    const PointSpec defaults[2] = {{72.0, -40.0, "Arctic (72N, 40W)"},
                                   {-70.0, 0.0, "Antarctic (-70S, 0E)"}};
    std::vector<PointSpec> points;
    const json::Value* list = body.find("points");
    if (list != nullptr && list->is_array() && list->items.size() >= 2) {
        for (std::size_t i = 0; i < 2; ++i) {
            const json::Value& item = list->items[i];
            const json::Value* lat = item.find("lat");
            const json::Value* lon = item.find("lon");
            if (lat == nullptr || lon == nullptr) throw BadRequest("each point needs lat and lon");
            points.push_back({lat->as_number(0.0), lon->as_number(0.0), defaults[i].fallback_label});
        }
    } else {
        points.push_back(defaults[0]);
        points.push_back(defaults[1]);
    }

    std::vector<ContrastRegion> regions;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const PointSpec& point = points[i];
        if (point.lat < -90.0 || point.lat > 90.0 || point.lon < -180.0 || point.lon > 180.0)
            throw BadRequest("a contrast latitude or longitude is out of range");
        std::string label = point.fallback_label;
        if (list != nullptr && list->is_array() && list->items.size() >= 2) {
            const json::Value* name = list->items[i].find("label");
            if (name != nullptr && !name->as_string("").empty()) label = name->as_string("");
        }
        regions.push_back(contrast_region(label, point.lat, point.lon, parameter, start, end));
    }

    const ContrastRegion& first = regions[0];
    const ContrastRegion& second = regions[1];
    const bool opposed = (first.slope_per_year > 0.0) != (second.slope_per_year > 0.0) &&
                         std::fabs(first.slope_per_year) > 1e-12 &&
                         std::fabs(second.slope_per_year) > 1e-12;
    const bool both_significant = first.significant && second.significant;

    std::string verdict;
    if (opposed && both_significant)
        verdict = "OPPOSITE SIGNIFICANT TRENDS: the same variable moves in opposite directions "
                  "in these two regions, so a single global sign would hide one of them.";
    else if (opposed)
        verdict = "OPPOSITE TRENDS, ONE NOT SIGNIFICANT: the signs disagree but at least one "
                  "region cannot be separated from no trend at this alpha.";
    else if (first.significant && second.significant)
        verdict = "SAME DIRECTION, BOTH SIGNIFICANT: both regions change the same way; the "
                  "difference is in magnitude, not sign.";
    else if (first.significant || second.significant)
        verdict = "ONE SIGNIFICANT REGION: only one of the two regions clears the significance "
                  "test at this alpha.";
    else
        verdict = "NO SIGNIFICANT TREND IN EITHER REGION: neither region can be distinguished "
                  "from a flat series at this alpha.";

    json::Value root = json::Value::object();
    root.set("parameter", json::Value::from_string(parameter));
    root.set("units", json::Value::from_string(first.units));
    root.set("start", json::Value::from_number(static_cast<double>(start_year)));
    root.set("end", json::Value::from_number(static_cast<double>(end_year)));
    root.set("alpha", json::Value::from_number(0.05));
    root.set("opposed", json::Value::from_bool(opposed));
    root.set("both_significant", json::Value::from_bool(both_significant));
    root.set("verdict", json::Value::from_string(verdict));
    json::Value out = json::Value::array();
    for (const auto& region : regions) out.push(contrast_point_json(region));
    root.set("regions", out);

    // The two short series are returned too, so the page can draw both curves.
    json::Value curves = json::Value::object();
    for (const auto& region : regions) {
        const Table table = parse_csv_text(read_text_file(region.file));
        json::Value series = json::Value::array();
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            double time = 0.0;
            double value = 0.0;
            if (!parse_number(table.cell(i, 0), time)) continue;
            if (!parse_number(table.cell(i, 1), value)) continue;
            json::Value sample = json::Value::array();
            sample.push(json::Value::from_number(time));
            sample.push(json::Value::from_number(value));
            series.push(sample);
        }
        curves.set(region.label, series);
    }
    root.set("curves", curves);
    return json_response(root);
}

Response visualize_request() {
    std::vector<std::string> files;
    try {
        for (const auto& entry : std::filesystem::directory_iterator("data")) {
            if (!entry.is_regular_file()) continue;
            const std::string name = entry.path().filename().string();
            if (name.size() > 4 && name.substr(name.size() - 4) == ".csv") files.push_back(name);
        }
    } catch (const std::exception&) {
    }
    std::sort(files.begin(), files.end());

    json::Value datasets = json::Value::array();
    for (const std::string& name : files) {
        Table table;
        try {
            table = parse_csv_text(read_text_file("data/" + name));
        } catch (const std::exception&) {
            continue;
        }
        if (table.rows.size() < 4 || table.header.size() < 2) continue;
        int value_column = table.column("value");
        if (value_column < 0) value_column = 1;
        std::vector<double> time, values;
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            double t = static_cast<double>(i), v = 0.0;
            parse_number(table.cell(i, 0), t);
            if (!parse_number(table.cell(i, value_column), v)) continue;
            time.push_back(t);
            values.push_back(v);
        }
        if (values.size() < 4) continue;

        const std::size_t step = std::max<std::size_t>(1, values.size() / 120);
        std::vector<double> thin, thin_time;
        for (std::size_t i = 0; i < values.size(); i += step) {
            thin.push_back(values[i]);
            thin_time.push_back(time[i]);
        }
        if (thin.size() < 2 || thin.back() != values.back()) {
            thin.push_back(values.back());
            thin_time.push_back(time.back());
        }

        const SenResult sen = sens_slope(time, values);
        double mean = 0.0, low = values[0], high = values[0];
        for (double v : values) {
            mean += v;
            low = std::min(low, v);
            high = std::max(high, v);
        }
        mean /= static_cast<double>(values.size());
        const double first = values.front();

        json::Value item = json::Value::object();
        item.set("file", json::Value::from_string(name));
        item.set("name", json::Value::from_string(name));
        item.set("variable", json::Value::from_string(
                                table.header.size() > 1 ? table.header[value_column] : "value"));
        item.set("samples", json::Value::from_number(static_cast<double>(values.size())));
        item.set("mean", json::Value::from_number(mean));
        item.set("minimum", json::Value::from_number(low));
        item.set("maximum", json::Value::from_number(high));
        item.set("slope", json::Value::from_number(sen.slope));
        item.set("percent", json::Value::from_number(
                               std::fabs(first) > 1e-12 ? 100.0 * (values.back() - first) / first
                                                        : 0.0));
        json::Value series = json::Value::array();
        for (double v : thin) series.push(json::Value::from_number(v));
        item.set("values", series);
        datasets.push(item);
    }

    json::Value root = json::Value::object();
    root.set("count", json::Value::from_number(
                           static_cast<double>(datasets.items.size())));
    root.set("datasets", datasets);
    return json_response(root);
}

Response grid_trend_request() {
    // A small deterministic demonstration field. It exists to exercise the
    // per-cell Sen slope and Getis-Ord path without pretending to be a NASA
    // download; the UI labels it as synthetic.
    constexpr int rows = 12;
    constexpr int columns = 24;
    constexpr int months = 120;
    const std::size_t cells = static_cast<std::size_t>(rows * columns);
    std::vector<double> time(months);
    for (int m = 0; m < months; ++m) time[static_cast<std::size_t>(m)] = 2000.0 + m / 12.0;

    std::vector<double> latitudes(rows), longitudes(columns);
    for (int i = 0; i < rows; ++i) latitudes[static_cast<std::size_t>(i)] = -82.5 + 15.0 * i;
    for (int j = 0; j < columns; ++j) longitudes[static_cast<std::size_t>(j)] = -172.5 + 15.0 * j;

    std::vector<std::vector<double>> series(cells, std::vector<double>(months));
    std::vector<double> slopes(cells);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            const std::size_t cell = static_cast<std::size_t>(i * columns + j);
            const double lat = latitudes[static_cast<std::size_t>(i)];
            const double lon = longitudes[static_cast<std::size_t>(j)];
            const double base = 289.0 - 0.52 * std::fabs(lat) + 2.2 * std::sin(deg_to_rad(lon));
            const double arctic = std::exp(-std::pow((lon - 35.0) / 55.0, 2.0)) *
                                  std::exp(-std::pow((lat - 72.0) / 24.0, 2.0));
            const double pacific = std::exp(-std::pow((lon + 150.0) / 60.0, 2.0)) *
                                   std::exp(-std::pow((lat + 10.0) / 30.0, 2.0));
            const double trend = 0.0015 + 0.012 * arctic - 0.004 * pacific;
            for (int m = 0; m < months; ++m) {
                const double seasonal = 0.9 * std::sin(2.0 * PI * m / 12.0 + deg_to_rad(lon) * 0.5);
                const double texture = 0.10 * std::sin(0.37 * m + 0.21 * i - 0.17 * j) +
                                       0.05 * std::cos(0.11 * m + 0.43 * j);
                series[cell][static_cast<std::size_t>(m)] =
                    base + trend * m + seasonal + texture;
            }
            slopes[cell] = sens_slope(time, series[cell]).slope;
        }
    }

    Weights weights(cells, std::vector<double>(cells, 0.0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            const std::size_t cell = static_cast<std::size_t>(i * columns + j);
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    const int ni = i + di;
                    const int nj = j + dj;
                    if (ni < 0 || ni >= rows || nj < 0 || nj >= columns) continue;
                    weights[cell][static_cast<std::size_t>(ni * columns + nj)] = 1.0;
                }
            }
        }
    }

    std::vector<double> gi(cells);
    double gi_mean = 0.0;
    for (std::size_t c = 0; c < cells; ++c) {
        gi[c] = getis_ord_gi(slopes, weights, c);
        gi_mean += gi[c];
    }
    gi_mean /= static_cast<double>(cells);
    double gi_sd = 0.0;
    for (double value : gi) gi_sd += (value - gi_mean) * (value - gi_mean);
    gi_sd = std::sqrt(gi_sd / static_cast<double>(cells));

    json::Value root = json::Value::object();
    root.set("rows", json::Value::from_number(rows));
    root.set("columns", json::Value::from_number(columns));
    root.set("months", json::Value::from_number(months));
    root.set("latitudes", number_array(latitudes));
    root.set("longitudes", number_array(longitudes));
    root.set("slopes", number_array(slopes));
    json::Value z_scores = json::Value::array();
    json::Value hotspots = json::Value::array();
    std::size_t hotspot_count = 0;
    for (double value : gi) {
        const double z = gi_sd > 0.0 ? (value - gi_mean) / gi_sd : 0.0;
        z_scores.push(json::Value::from_number(z));
        const bool hotspot = z >= 1.96;
        hotspots.push(json::Value::from_bool(hotspot));
        if (hotspot) ++hotspot_count;
    }
    root.set("gi_z", z_scores);
    root.set("hotspots", hotspots);
    root.set("hotspot_count", json::Value::from_number(static_cast<double>(hotspot_count)));
    root.set("note", json::Value::from_string(
                          "Synthetic 24 x 12 demonstration field; not a NASA product and not the selected series."));
    return json_response(root);
}

Response nasa_datasets_request() {
    json::Value datasets = json::Value::array();
    for (const auto& dataset : remote_datasets()) {
        json::Value item = json::Value::object();
        item.set("key", json::Value::from_string(dataset.key));
        item.set("catalog_id", json::Value::from_string(dataset.catalog_id));
        item.set("short_name", json::Value::from_string(dataset.short_name));
        item.set("version", json::Value::from_string(dataset.version));
        item.set("provider", json::Value::from_string(dataset.provider));
        item.set("variable", json::Value::from_string(dataset.variable));
        item.set("units", json::Value::from_string(dataset.units));
        item.set("mode", json::Value::from_string(dataset.mode));
        item.set("description", json::Value::from_string(dataset.description));
        datasets.push(item);
    }

    std::vector<std::string> files;
    try {
        const std::filesystem::path directory("data/nasa");
        if (std::filesystem::exists(directory)) {
            for (const auto& entry : std::filesystem::directory_iterator(directory)) {
                if (entry.is_regular_file())
                    files.push_back(entry.path().filename().string());
            }
        }
    } catch (const std::exception&) {
    }
    std::sort(files.begin(), files.end());

    json::Value root = json::Value::object();
    root.set("datasets", datasets);
    root.set("files", string_array(files));
    root.set("python", json::Value::from_string(resolve_python_interpreter()));
    root.set("wrapper", json::Value::from_string(resolve_fetch_wrapper()));
    return json_response(root);
}

Response nasa_plan_request(const json::Value& body) {
    const std::string key = value_string(body, "key");
    const RemoteDataset* dataset = find_remote_dataset(key);
    if (dataset == nullptr) throw BadRequest("unknown NASA dataset key");

    const std::string start = value_string(body, "start", "2024-01-01");
    const std::string end = value_string(body, "end", "2024-12-31");
    const int maximum = std::max(1, static_cast<int>(value_number(body, "max", 4.0)));
    const std::string output_csv = "data/nasa/" + dataset->key + "_monthly.csv";

    std::vector<std::string> arguments = {resolve_python_interpreter(),
                                          resolve_fetch_wrapper(),
                                          "fetch",
                                          "--short-name",
                                          dataset->short_name,
                                          "--start",
                                          start,
                                          "--end",
                                          end,
                                          "--max",
                                          std::to_string(maximum),
                                          "--out",
                                          "data/nasa",
                                          "--out-csv",
                                          output_csv,
                                          "--variable",
                                          dataset->variable,
                                          "--mode",
                                          dataset->mode,
                                          "--units",
                                          dataset->units,
                                          "--dry-run"};
    if (!dataset->version.empty()) {
        arguments.push_back("--version");
        arguments.push_back(dataset->version);
    }

    const std::string command = build_command(arguments);
    const CommandResult result = run_command(command);
    json::Value root = json::Value::object();
    root.set("key", json::Value::from_string(dataset->key));
    root.set("short_name", json::Value::from_string(dataset->short_name));
    root.set("command", json::Value::from_string(command));
    root.set("exit_code", json::Value::from_number(result.exit_code));
    root.set("output", json::Value::from_string(result.output));
    root.set("output_csv", json::Value::from_string(output_csv));
    return json_response(root);
}

Response nasa_search_request(const json::Value& body) {
    const std::string key = value_string(body, "key");
    const RemoteDataset* dataset = find_remote_dataset(key);
    if (dataset == nullptr) throw BadRequest("unknown NASA dataset key");
    const std::string start = value_string(body, "start", "2024-01-01");
    const std::string end = value_string(body, "end", "2024-12-31");
    const int maximum = std::max(1, static_cast<int>(value_number(body, "max", 5.0)));

    std::vector<std::string> arguments = {resolve_python_interpreter(),
                                          "tools/nasa_cmr.py",
                                          "--short-name",
                                          dataset->short_name,
                                          "--start",
                                          start,
                                          "--end",
                                          end,
                                          "--max",
                                          std::to_string(maximum)};
    if (!dataset->version.empty()) {
        arguments.push_back("--version");
        arguments.push_back(dataset->version);
    }
    const std::string command = build_command(arguments);
    const CommandResult result = run_command(command);

    json::Value root = json::Value::object();
    root.set("key", json::Value::from_string(dataset->key));
    root.set("short_name", json::Value::from_string(dataset->short_name));
    root.set("command", json::Value::from_string(command));
    root.set("exit_code", json::Value::from_number(result.exit_code));
    root.set("output", json::Value::from_string(result.output));
    return json_response(root);
}

Response query_request(const json::Value& body) {
    const std::string text = value_string(body, "text");
    if (text.empty()) throw BadRequest("query text is required");
    const Intent intent = parse_intent_text(text, 2025);

    json::Value root = json::Value::object();
    root.set("valid", json::Value::from_bool(intent.valid));
    root.set("error", json::Value::from_string(intent.error));
    json::Value intent_json = json::Value::object();
    intent_json.set("process", json::Value::from_string(intent.process));
    intent_json.set("region", json::Value::from_string(intent.region));
    intent_json.set("start_year", json::Value::from_number(intent.start_year));
    intent_json.set("end_year", json::Value::from_number(intent.end_year));
    intent_json.set("variable", json::Value::from_string(intent.variable));
    root.set("intent", intent_json);

    json::Value variables = json::Value::array();
    json::Value reasons = json::Value::array();
    if (intent.valid) {
        const Selection selection = select_variables(intent);
        for (const auto& entry : selection.entries) {
            json::Value item = json::Value::object();
            item.set("id", json::Value::from_string(entry.id));
            item.set("name", json::Value::from_string(entry.name));
            item.set("units", json::Value::from_string(entry.units));
            variables.push(item);
        }
        for (const auto& reason : selection.reasons) reasons.push(json::Value::from_string(reason));
    }
    root.set("variables", variables);
    root.set("reasons", reasons);
    return json_response(root);
}

Response route(const Request& request) {
    if (request.method == "GET") {
        if (request.path == "/" || request.path == "/index.html") {
            Response response;
            response.content_type = "text/html; charset=utf-8";
            response.body = climate_ui_html();
            return response;
        }
        if (request.path == "/favicon.ico") {
            Response response;
            response.status = 204;
            response.body.clear();
            return response;
        }
        if (request.path == "/api/overview") return overview_request();
        if (request.path == "/api/catalog") return catalog_request();
        if (request.path == "/api/help") return help_request();
        if (request.path == "/api/formulas") return formulas_request();
        if (request.path == "/api/lab") return lab_catalog_request();
        if (request.path == "/api/visualize") return visualize_request();
        if (request.path == "/api/live/status") return live_status_request();
        if (request.path == "/api/grid-trend") return grid_trend_request();
        if (request.path == "/api/nasa/datasets") return nasa_datasets_request();
        if (request.path == "/api/data/file") {
            const std::string name = query_parameter(request.query, "name");
            if (name.empty() || name.find('/') != std::string::npos ||
                name.find('\\') != std::string::npos || name.find("..") != std::string::npos ||
                name.size() < 4 || name.substr(name.size() - 4) != ".csv")
                throw BadRequest("invalid data file name");
            Response response;
            response.content_type = "text/csv; charset=utf-8";
            response.body = read_text_file("data/" + name);
            return response;
        }
        throw BadRequest("unknown endpoint");
    }

    if (request.method == "POST") {
        json::Value body;
        try {
            body = json::parse(request.body.empty() ? "{}" : request.body);
        } catch (const std::exception&) {
            throw BadRequest("request body is not valid JSON");
        }
        if (!body.is_object()) throw BadRequest("request body must be a JSON object");
        if (request.path == "/api/analyze") return analyze_series_request(body);
        if (request.path == "/api/co2") return analyze_co2_request(body);
        if (request.path == "/api/retrieve/sst") return sst_request(body);
        if (request.path == "/api/retrieve/index") return vegetation_index_request(body);
        if (request.path == "/api/retrieve/seaice") return seaice_request(body);
        if (request.path == "/api/retrieve/precip") return precip_request(body);
        if (request.path == "/api/info") return information_request(body);
        if (request.path == "/api/sstunc") return uncertainty_request(body);
        if (request.path == "/api/physics") return physics_request(body);
        if (request.path == "/api/lab/run") return lab_run_request(body);
        if (request.path == "/api/live/pull") return live_pull_request(body);
        if (request.path == "/api/contrast") return contrast_request(body);
        if (request.path == "/api/query") return query_request(body);
        if (request.path == "/api/run") return run_request(body);
        if (request.path == "/api/nasa/plan") return nasa_plan_request(body);
        if (request.path == "/api/nasa/search") return nasa_search_request(body);
        throw BadRequest("unknown endpoint");
    }

    Response response;
    response.status = 405;
    response.body = "{\"error\":\"method not allowed\"}";
    return response;
}

bool send_all(int socket, const char* data, std::size_t size) {
    std::size_t sent = 0;
    while (sent < size) {
        const ssize_t n = send(socket, data + sent, size - sent, 0);
        if (n <= 0) return false;
        sent += static_cast<std::size_t>(n);
    }
    return true;
}

void handle_connection(int socket) {
    std::string buffer;
    char chunk[8192];
    std::size_t header_end = std::string::npos;
    while (header_end == std::string::npos) {
        const ssize_t n = recv(socket, chunk, sizeof(chunk), 0);
        if (n <= 0) {
            close(socket);
            return;
        }
        buffer.append(chunk, static_cast<std::size_t>(n));
        header_end = buffer.find("\r\n\r\n");
        if (buffer.size() > 2 * 1024 * 1024) {
            close(socket);
            return;
        }
    }

    Request request;
    const std::string header = buffer.substr(0, header_end);
    std::istringstream lines(header);
    std::string line;
    std::getline(lines, line);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    {
        std::istringstream first(line);
        std::string target;
        first >> request.method >> target;
        const std::size_t question = target.find('?');
        request.path = url_decode(target.substr(0, question));
        if (question != std::string::npos) request.query = target.substr(question + 1);
    }
    while (std::getline(lines, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::size_t colon = line.find(':');
        if (colon == std::string::npos) continue;
        request.headers[lower_copy(trim_copy(line.substr(0, colon)))] =
            trim_copy(line.substr(colon + 1));
    }

    const auto expect = request.headers.find("expect");
    if (expect != request.headers.end() &&
        lower_copy(expect->second).find("100-continue") != std::string::npos) {
        const std::string cont = "HTTP/1.1 100 Continue\r\n\r\n";
        send_all(socket, cont.data(), cont.size());
    }

    const std::size_t body_start = header_end + 4;
    const auto transfer = request.headers.find("transfer-encoding");
    const bool chunked = transfer != request.headers.end() &&
                         lower_copy(transfer->second).find("chunked") != std::string::npos;

    const auto receive_more = [&]() -> bool {
        const ssize_t n = recv(socket, chunk, sizeof(chunk), 0);
        if (n <= 0) return false;
        buffer.append(chunk, static_cast<std::size_t>(n));
        return true;
    };

    if (chunked) {
        std::size_t position = body_start;
        while (true) {
            std::size_t line_end = buffer.find("\r\n", position);
            while (line_end == std::string::npos) {
                if (!receive_more() || buffer.size() > 40 * 1024 * 1024) break;
                line_end = buffer.find("\r\n", position);
            }
            if (line_end == std::string::npos) break;
            const std::string size_text = buffer.substr(position, line_end - position);
            const std::size_t extension = size_text.find(';');
            std::size_t size = 0;
            try {
                size = static_cast<std::size_t>(
                    std::stoull(size_text.substr(0, extension), nullptr, 16));
            } catch (const std::exception&) {
                request.body.clear();
                break;
            }
            position = line_end + 2;
            while (buffer.size() < position + size + 2) {
                if (!receive_more() || buffer.size() > 40 * 1024 * 1024) break;
            }
            if (buffer.size() < position + size) break;
            request.body.append(buffer, position, size);
            position += size;
            if (position + 2 <= buffer.size()) position += 2;
            if (size == 0) break;
            if (request.body.size() > 32 * 1024 * 1024) {
                request.body.clear();
                break;
            }
        }
    } else {
        std::size_t content_length = 0;
        const auto length = request.headers.find("content-length");
        if (length != request.headers.end()) {
            try {
                content_length = static_cast<std::size_t>(std::stoull(length->second));
            } catch (const std::exception&) {
                content_length = 0;
            }
        }
        if (content_length > 32 * 1024 * 1024) {
            const std::string response =
                "HTTP/1.1 413 Payload Too Large\r\nContent-Length: 0\r\nConnection: close\r\n\r\n";
            send_all(socket, response.data(), response.size());
            close(socket);
            return;
        }
        request.body = buffer.substr(body_start);
        while (request.body.size() < content_length) {
            if (!receive_more()) break;
            request.body = buffer.substr(body_start);
        }
        if (request.body.size() > content_length) request.body.resize(content_length);
    }

    Response response;
    try {
        response = route(request);
    } catch (const BadRequest& exception) {
        response.status = 400;
        response.body = "{\"error\":" + json::escape(exception.what()) + "}";
    } catch (const std::exception& exception) {
        response.status = 500;
        response.body = "{\"error\":" + json::escape(exception.what()) + "}";
    }

    std::ostringstream head;
    head << "HTTP/1.1 " << response.status << ' ' << status_text(response.status) << "\r\n"
         << "Content-Type: " << response.content_type << "\r\n"
         << "Content-Length: " << response.body.size() << "\r\n"
         << "Cache-Control: no-store\r\n"
         << "X-Content-Type-Options: nosniff\r\n"
         << "Connection: close\r\n\r\n";
    const std::string head_text = head.str();
    send_all(socket, head_text.data(), head_text.size());
    if (!response.body.empty()) send_all(socket, response.body.data(), response.body.size());
    close(socket);
}

}  // namespace

int run_ui_server(const UiOptions& options) {
    if (options.port <= 0 || options.port > 65535) {
        std::cerr << "error: --port must be between 1 and 65535\n";
        return 1;
    }

    signal(SIGPIPE, SIG_IGN);

    const int server = socket(AF_INET, SOCK_STREAM, 0);
    if (server < 0) {
        std::cerr << "error: could not create a socket\n";
        return 1;
    }
    int reuse = 1;
    setsockopt(server, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse));

    sockaddr_in address{};
    address.sin_family = AF_INET;
    address.sin_port = htons(static_cast<unsigned short>(options.port));
    if (options.host == "0.0.0.0") {
        address.sin_addr.s_addr = INADDR_ANY;
    } else if (inet_pton(AF_INET, options.host.c_str(), &address.sin_addr) != 1) {
        std::cerr << "error: invalid IPv4 address " << options.host << "\n";
        close(server);
        return 1;
    }

    if (bind(server, reinterpret_cast<sockaddr*>(&address), sizeof(address)) != 0) {
        std::cerr << "error: could not bind " << options.host << ":" << options.port
                  << " (is the port already in use?)\n";
        close(server);
        return 1;
    }
    if (listen(server, 32) != 0) {
        std::cerr << "error: could not listen on the socket\n";
        close(server);
        return 1;
    }

    std::cout << "Climate Workbench running at http://" << options.host << ':' << options.port
              << "\nPress Ctrl+C to stop.\n";
    std::cout.flush();

    while (true) {
        sockaddr_in client{};
        socklen_t length = sizeof(client);
        const int connection =
            accept(server, reinterpret_cast<sockaddr*>(&client), &length);
        if (connection < 0) continue;
        std::thread([connection]() { handle_connection(connection); }).detach();
    }
}

}  // namespace climate
