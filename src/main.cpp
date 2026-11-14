#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "analysis.hpp"
#include "catalog.hpp"
#include "common.hpp"
#include "csv.hpp"
#include "harmonic.hpp"
#include "fetch.hpp"
#include "intent.hpp"
#include "json.hpp"
#include "netcdf.hpp"
#include "physics.hpp"
#include "report.hpp"
#include "retrieval.hpp"
#include "series.hpp"
#include "stats.hpp"
#include "trend.hpp"

using namespace climate;

namespace {

struct Arguments {
    std::map<std::string, std::string> options;
    std::vector<std::string> positional;

    bool has(const std::string& key) const { return options.find(key) != options.end(); }

    std::string get(const std::string& key, const std::string& fallback = "") const {
        const auto it = options.find(key);
        return it == options.end() ? fallback : it->second;
    }

    bool flag(const std::string& key) const {
        const auto it = options.find(key);
        if (it == options.end()) return false;
        return it->second != "false" && it->second != "0";
    }

    double number(const std::string& key, double fallback) const {
        const auto it = options.find(key);
        if (it == options.end()) return fallback;
        try {
            return std::stod(it->second);
        } catch (const std::exception&) {
            return fallback;
        }
    }

    int integer(const std::string& key, int fallback) const {
        const auto it = options.find(key);
        if (it == options.end()) return fallback;
        try {
            return std::stoi(it->second);
        } catch (const std::exception&) {
            return fallback;
        }
    }

    std::string text() const {
        std::ostringstream out;
        for (std::size_t i = 0; i < positional.size(); ++i) {
            if (i > 0) out << ' ';
            out << positional[i];
        }
        return out.str();
    }
};

Arguments parse_arguments(int argc, char** argv, int start) {
    Arguments args;
    for (int i = start; i < argc; ++i) {
        const std::string arg = argv[i];
        const bool looks_like_option = arg.size() > 2 && arg[0] == '-' && arg[1] == '-';
        if (looks_like_option) {
            std::string key = arg.substr(2);
            const std::size_t equals = key.find('=');
            if (equals != std::string::npos) {
                args.options[key.substr(0, equals)] = key.substr(equals + 1);
            } else {
                const bool next_is_option =
                    i + 1 < argc && std::string(argv[i + 1]).size() > 2 &&
                    std::string(argv[i + 1])[0] == '-' && std::string(argv[i + 1])[1] == '-';
                if (i + 1 < argc && !next_is_option)
                    args.options[key] = argv[++i];
                else
                    args.options[key] = "true";
            }
        } else {
            args.positional.push_back(arg);
        }
    }
    return args;
}

std::string lowered(std::string text) {
    for (auto& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return text;
}

std::string read_text_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open file: " + path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::string format_name(const Arguments& args) {
    std::string format = lowered(args.get("format", "text"));
    if (format == "md") format = "markdown";
    if (format != "text" && format != "json" && format != "markdown")
        throw std::runtime_error("unknown format '" + format + "'; use text, json or markdown");
    return format;
}

void emit(const std::string& content, const std::string& path) {
    if (path.empty()) {
        std::cout << content;
        return;
    }
    write_text_file(path, content);
    std::cout << "wrote " << content.size() << " bytes to " << path << "\n";
}

bool load_series(const std::string& path, const std::string& time_column,
                 const std::string& value_column, std::vector<double>& time,
                 std::vector<double>& values, std::string& error) {
    try {
        const Table table = read_csv(path);
        int time_index = table.column(time_column.empty() ? "time" : time_column);
        if (time_index < 0) {
            if (time_column.empty() && !table.header.empty()) {
                time_index = 0;
            } else {
                error = "time column '" + (time_column.empty() ? std::string("time") : time_column) +
                        "' not found in " + path;
                return false;
            }
        }

        int value_index = -1;
        if (!value_column.empty()) {
            value_index = table.column(value_column);
            if (value_index < 0) {
                error = "value column '" + value_column + "' not found in " + path;
                return false;
            }
        } else {
            value_index = table.column("value");
            if (value_index < 0 && table.header.size() > 1) value_index = 1;
            if (value_index < 0) value_index = 0;
        }

        time.clear();
        values.clear();
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            double t = 0.0;
            if (!parse_number(table.cell(i, time_index), t)) continue;
            double v = std::numeric_limits<double>::quiet_NaN();
            parse_number(table.cell(i, value_index), v);
            time.push_back(t);
            values.push_back(v);
        }
        if (time.size() < 3) {
            error = "need at least three numeric observations in " + path;
            return false;
        }
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

SeriesOptions options_from(const Arguments& args) {
    SeriesOptions options;
    options.period = args.integer("period", 12);
    options.seasonal_window = args.integer("seasonal-window", std::max(7, options.period + 1));
    options.trend_window = args.integer("trend-window", std::max(21, 2 * options.period + 1));
    options.low_pass_window = args.integer("low-pass-window", std::max(7, options.period + 1));
    options.seasonal_degree = args.integer("seasonal-degree", 1);
    options.trend_degree = args.integer("trend-degree", 1);
    options.inner_iterations = args.integer("inner", 2);
    options.outer_iterations = args.integer("outer", 5);
    options.alpha = args.number("alpha", 0.05);
    options.hamed_rao = !args.flag("no-hamed-rao");
    options.prewhiten = args.flag("prewhiten");
    options.test_deseasonalised = !args.flag("no-detrend");
    return options;
}

std::string join(const std::vector<std::string>& values, const std::string& separator) {
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) out << separator;
        out << values[i];
    }
    return out.str();
}

void print_usage() {
    std::cout
        << "climate - NASA-style climate variable retrieval and trend analysis\n\n"
        << "Usage: climate <command> [options]\n\n"
        << "Commands\n"
        << "  catalog [--process P] [--region R] [--format text|json]\n"
        << "      List the variable catalog.\n\n"
        << "  query \"<text>\" [--intent FILE] [--data DIR] [--format text|json]\n"
        << "      Parse a free-text question into the intent IR, compile it to catalog\n"
        << "      variables, and optionally run the analysis on DIR/<id>_monthly.csv.\n\n"
        << "  intent --file IR.json [--format text|json]\n"
        << "      Validate an intent JSON produced by an external language model.\n\n"
        << "  analyze --file DATA.csv [--time COL] [--value COL] [--period N]\n"
        << "          [--alpha A] [--seasonal-window N] [--trend-window N]\n"
        << "          [--low-pass-window N] [--prewhiten] [--no-hamed-rao] [--no-detrend]\n"
        << "          [--format text|json|markdown] [--out FILE] [--label NAME]\n"
        << "      STL decomposition, Mann-Kendall test and Sen's slope.\n\n"
        << "  co2 --file DATA.csv [--time COL] [--value COL] [--poly N] [--harmonics N]\n"
        << "          [--format text|json] [--out FILE]\n"
        << "      NOAA/CMDL polynomial plus harmonic curve fit and CO2 growth rate.\n\n"
        << "  retrieve sst --file MATCHUPS.csv [--out MODEL.csv]\n"
        << "  retrieve sst --coeffs MODEL.csv --bt11 K --bt12 K --first-guess K --zenith DEG\n"
        << "          [--mirror M] [--theta-star DEG]\n"
        << "      Calibrate or apply the NLSST sea surface temperature retrieval.\n\n"
        << "  retrieve ndvi --file DATA.csv --red COL --nir COL [--out FILE]\n"
        << "      Compute NDVI = (NIR - RED) / (NIR + RED).\n\n"
        << "  retrieve seaice --file GRID.csv [--threshold 0.15] [--format text|json]\n"
        << "      Sea ice extent (15% rule) and area from a concentration grid.\n\n"
        << "  retrieve seaice-pixel --tb19v K --tb19h K --tb37v K [--tie-points FILE]\n"
        << "      NASA Team concentration for a single passive-microwave footprint.\n\n"
        << "  retrieve precip --file ESTIMATES.csv [--cold-fraction F] [--water-vapour W]\n"
        << "          [--a1 A] [--a2 B]\n"
        << "      GPI cold-cloud index, water-vapour adjustment and GPCP merge.\n\n"
        << "  calibrate nlsst --file MATCHUPS.csv --out MODEL.csv\n"
        << "  calibrate seaice [--tie-points FILE] --out MODEL.csv\n"
        << "      Write calibrated retrieval coefficients.\n\n"
        << "  physics toa --forcing F --feedback L --delta-t T\n"
        << "  physics surface --albedo A --shortwave S --emissivity-air EA --longwave L\n"
        << "          --emissivity-surface ES --temperature T\n"
        << "  physics carbon --start-ppm P1 --end-ppm P2 --start-year Y1 --end-year Y2\n"
        << "  physics amoc --file VELOCITY.csv [--depth M]\n"
        << "  physics feedback [--planck L] [--water-vapour L] [--lapse-rate L] [--albedo L]\n"
        << "          [--cloud L] [--forcing F]\n"
        << "  physics forcing-co2 [--ppm C] [--reference-ppm C0]\n"
        << "  physics mie --size-parameter X [--real N] [--imag K]\n"
        << "  physics aod --number-concentration N --radius-um R --path-length L [--qext Q]\n"
        << "  physics wind-stress --u10 U --v10 V [--air-density RHO] [--drag CD]\n"
        << "  physics ekman --tau-x TX --tau-y TY --latitude LAT [--density RHO]\n"
        << "  physics geostrophic --density RHO --dp-dx DX --dp-dy DY --latitude LAT\n"
        << "  physics thermal-wind --density RHO --drho-dx DX --drho-dy DY --latitude LAT\n"
        << "      Coupled-system diagnostics and radiative-transfer formulas.\n\n"
        << "  retrieve index --file DATA.csv --red COL --nir COL [--blue COL]\n"
        << "          [--index ndvi|evi|savi] [--soil L] [--out FILE]\n"
        << "      NDVI, EVI or SAVI from reflectance pairs.\n\n"
        << "  data list [--format text|json]\n"
        << "  data search --dataset KEY --start YYYY-MM-DD --end YYYY-MM-DD [--bbox W,S,E,N]\n"
        << "  data fetch --dataset KEY --start ... --end ... [--out DIR] [--max N] [--dry-run]\n"
        << "  data convert --dataset KEY --out CSV FILE.nc [FILE.nc ...]\n"
        << "  data inspect --file FILE.nc\n"
        << "  data export --file FILE.nc --variable V --out CSV [--mode timeseries|grid]\n"
        << "      NASA Earthdata acquisition wrapper and NetCDF-3 inspection.\n\n"
        << "  help\n"
        << "      Show this message.\n";
}

int command_catalog(const Arguments& args) {
    const std::string process = args.get("process");
    const std::string region = args.get("region");
    const std::string format = format_name(args);

    std::vector<const CatalogEntry*> selected;
    for (const auto& entry : catalog_entries()) {
        if (!process.empty() &&
            std::find(entry.processes.begin(), entry.processes.end(), process) == entry.processes.end())
            continue;
        if (!region.empty() &&
            std::find(entry.regions.begin(), entry.regions.end(), region) == entry.regions.end())
            continue;
        selected.push_back(&entry);
    }

    if (format == "json") {
        json::Value list = json::Value::array();
        for (const auto* entry : selected) {
            json::Value item = json::Value::object();
            item.set("id", json::Value::from_string(entry->id));
            item.set("name", json::Value::from_string(entry->name));
            item.set("units", json::Value::from_string(entry->units));
            json::Value processes = json::Value::array();
            for (const auto& value : entry->processes) processes.push(json::Value::from_string(value));
            json::Value regions = json::Value::array();
            for (const auto& value : entry->regions) regions.push(json::Value::from_string(value));
            item.set("processes", processes);
            item.set("regions", regions);
            item.set("retrieval", json::Value::from_string(entry->retrieval));
            item.set("analysis", json::Value::from_string(entry->analysis_method));
            item.set("description", json::Value::from_string(entry->description));
            list.push(item);
        }
        std::cout << json::dump(list, 2) << "\n";
        return 0;
    }

    if (selected.empty()) {
        std::cerr << "no catalog variables match the filter\n";
        return 1;
    }
    for (const auto* entry : selected) {
        std::cout << entry->id << " - " << entry->name << " [" << entry->units << "]\n";
        std::cout << "  processes: " << join(entry->processes, ", ") << "\n";
        std::cout << "  regions  : " << join(entry->regions, ", ") << "\n";
        std::cout << "  retrieval: " << entry->retrieval << "\n";
        std::cout << "  analysis : " << entry->analysis_method << "\n";
        std::cout << "  purpose  : " << entry->description << "\n\n";
    }
    return 0;
}

void run_selection_on_data(const CatalogEntry& entry, const std::string& data_dir) {
    const std::string path = data_dir + "/" + entry.id + "_monthly.csv";
    std::ifstream probe(path);
    if (!probe) {
        std::cout << "\n[" << entry.id << "] no data file at " << path << "; selection only\n";
        return;
    }
    probe.close();

    std::vector<double> time;
    std::vector<double> values;
    std::string error;
    if (!load_series(path, "", "", time, values, error)) {
        std::cout << "\n[" << entry.id << "] " << error << "\n";
        return;
    }
    std::cout << "\n";
    if (entry.id == "co2" || entry.id == "ch4") {
        CO2Options options;
        options.harmonic_count = entry.id == "ch4" ? 2 : 4;
        CO2Analysis analysis = analyze_co2(time, values, options);
        std::cout << co2_report_text(analysis, entry.id, entry.units, path);
    } else {
        SeriesAnalysis analysis = analyze_series(time, values, SeriesOptions{});
        std::cout << series_report_text(analysis, entry.id, entry.units, path);
    }
}

int command_query(const Arguments& args) {
    const int reference_year = args.integer("reference-year", 2025);
    const std::string format = format_name(args);

    Intent intent;
    if (args.has("intent")) {
        intent = parse_intent_json(read_text_file(args.get("intent")), reference_year);
    } else {
        const std::string text = args.text();
        if (text.empty()) {
            std::cerr << "error: query text is required, for example: climate query \"is the arctic "
                         "warming?\"\n";
            return 1;
        }
        intent = parse_intent_text(text, reference_year);
    }
    if (!intent.valid) {
        std::cerr << "error: " << intent.error << "\n";
        return 1;
    }

    const Selection selection = select_variables(intent);

    if (format == "json") {
        json::Value root = json::Value::object();
        root.set("intent", json::parse(intent_to_json(intent)));
        json::Value entries = json::Value::array();
        for (const auto& entry : selection.entries) {
            json::Value item = json::Value::object();
            item.set("id", json::Value::from_string(entry.id));
            item.set("name", json::Value::from_string(entry.name));
            item.set("units", json::Value::from_string(entry.units));
            item.set("retrieval", json::Value::from_string(entry.retrieval));
            item.set("analysis", json::Value::from_string(entry.analysis_method));
            entries.push(item);
        }
        root.set("selection", entries);
        json::Value reasons = json::Value::array();
        for (const auto& reason : selection.reasons) reasons.push(json::Value::from_string(reason));
        root.set("reasons", reasons);
        std::cout << json::dump(root, 2) << "\n";
    } else {
        std::cout << "Intent IR\n" << intent_to_json(intent) << "\n\n";
        if (selection.entries.empty()) {
            std::cout << "No catalog variables matched.\n";
        } else {
            std::cout << "Selected variables (" << selection.entries.size() << ")\n";
            for (const auto& entry : selection.entries)
                std::cout << "  " << entry.id << " - " << entry.name << " [" << entry.units << "]\n";
            std::cout << "\nWhy\n";
            for (const auto& reason : selection.reasons) std::cout << "  - " << reason << "\n";
        }
    }

    const std::string data_dir = args.get("data");
    if (!data_dir.empty() && format != "json") {
        for (const auto& entry : selection.entries) run_selection_on_data(entry, data_dir);
    }
    return 0;
}

int command_intent(const Arguments& args) {
    const int reference_year = args.integer("reference-year", 2025);
    const std::string format = format_name(args);
    std::string text;
    if (args.has("file"))
        text = read_text_file(args.get("file"));
    else
        text = args.text();
    if (text.empty()) {
        std::cerr << "error: provide --file IR.json or the JSON text as an argument\n";
        return 1;
    }

    const Intent intent = parse_intent_json(text, reference_year);
    if (!intent.valid) {
        std::cerr << "error: " << intent.error << "\n";
        return 1;
    }
    const Selection selection = select_variables(intent);

    if (format == "json") {
        json::Value root = json::Value::object();
        root.set("intent", json::parse(intent_to_json(intent)));
        json::Value entries = json::Value::array();
        for (const auto& entry : selection.entries)
            entries.push(json::Value::from_string(entry.id));
        root.set("variables", entries);
        std::cout << json::dump(root, 2) << "\n";
        return 0;
    }

    std::cout << intent_to_json(intent) << "\n\nVariables: ";
    std::vector<std::string> ids;
    for (const auto& entry : selection.entries) ids.push_back(entry.id);
    std::cout << (ids.empty() ? std::string("(none)") : join(ids, ", ")) << "\n";
    return 0;
}

int command_analyze(const Arguments& args) {
    const std::string path = args.get("file");
    if (path.empty()) {
        std::cerr << "error: --file is required\n";
        return 1;
    }
    const std::string format = format_name(args);
    const std::string time_column = args.get("time");
    const std::string value_column = args.get("value");

    std::vector<double> time;
    std::vector<double> values;
    std::string error;
    if (!load_series(path, time_column, value_column, time, values, error)) {
        std::cerr << "error: " << error << "\n";
        return 1;
    }

    SeriesOptions options = options_from(args);
    SeriesAnalysis analysis = analyze_series(time, values, options);
    analysis.label = args.get("label");

    const std::string variable = args.get("variable", args.get("label", "series"));
    const std::string units = args.get("units", "");

    std::string content;
    if (format == "json")
        content = series_report_json(analysis, variable, units, path);
    else if (format == "markdown")
        content = series_report_markdown(analysis, variable, units, path);
    else
        content = series_report_text(analysis, variable, units, path);
    if (format != "json") content += "\n";
    emit(content, args.get("out"));
    return 0;
}

int command_co2(const Arguments& args) {
    const std::string path = args.get("file");
    if (path.empty()) {
        std::cerr << "error: --file is required\n";
        return 1;
    }
    const std::string format = format_name(args);

    std::vector<double> time;
    std::vector<double> values;
    std::string error;
    if (!load_series(path, args.get("time"), args.get("value"), time, values, error)) {
        std::cerr << "error: " << error << "\n";
        return 1;
    }

    CO2Options options;
    options.polynomial_order = args.integer("poly", 2);
    options.harmonic_count = args.integer("harmonics", 4);
    options.reject_iterations = args.integer("reject-iterations", 4);
    options.reject_threshold = args.number("reject-threshold", 3.0);
    options.alpha = args.number("alpha", 0.05);

    CO2Analysis analysis = analyze_co2(time, values, options);
    analysis.label = args.get("label");

    const std::string variable = args.get("variable", "co2");
    const std::string units = args.get("units", "ppm");
    if (format == "json") {
        emit(co2_report_json(analysis, variable, units, path), args.get("out"));
    } else {
        emit(co2_report_text(analysis, variable, units, path) + "\n", args.get("out"));
    }
    return 0;
}

int retrieve_sst(const Arguments& args) {
    NLSSTModel model;
    bool have_model = false;

    if (args.has("file")) {
        std::vector<NLSSTSample> samples;
        std::string error;
        if (!read_nlsst_samples(args.get("file"), samples, error)) {
            std::cerr << "error: " << error << "\n";
            return 1;
        }
        model = calibrate_nlsst(samples);
        have_model = true;
        std::cout << "NLSST calibration from " << samples.size() << " matchups\n";
        const char* names[7] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};
        for (std::size_t i = 0; i < model.coefficients.size(); ++i)
            std::cout << "  " << names[i] << " = " << model.coefficients[i] << "\n";
        std::cout << "  RMSE = " << model.rmse << " K\n";
        std::cout << "  bias = " << model.bias << " K\n";
    }

    if (args.has("coeffs")) {
        std::string error;
        if (!read_nlsst_model(args.get("coeffs"), model, error)) {
            std::cerr << "error: " << error << "\n";
            return 1;
        }
        have_model = true;
        std::cout << "Loaded NLSST coefficients from " << args.get("coeffs") << "\n";
        const char* names[7] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};
        for (std::size_t i = 0; i < model.coefficients.size(); ++i)
            std::cout << "  " << names[i] << " = " << model.coefficients[i] << "\n";
        std::cout << "  RMSE = " << model.rmse << " K\n";
        std::cout << "  bias = " << model.bias << " K\n";
    }

    if (args.has("bt11")) {
        if (!have_model) {
            std::cerr << "error: provide --file MATCHUPS.csv or --coeffs MODEL.csv\n";
            return 1;
        }
        const double bt11 = args.number("bt11", std::nan(""));
        const double bt12 = args.number("bt12", std::nan(""));
        const double first_guess = args.number("first-guess", std::nan(""));
        const double zenith = args.number("zenith", std::nan(""));
        if (!std::isfinite(bt11) || !std::isfinite(bt12) || !std::isfinite(first_guess) ||
            !std::isfinite(zenith)) {
            std::cerr << "error: --bt11, --bt12, --first-guess and --zenith are required\n";
            return 1;
        }
        const double mirror = args.number("mirror", 0.0);
        const double theta_star = deg_to_rad(args.number("theta-star", 0.0));
        const double retrieved =
            nlsst_retrieve(model, bt11, bt12, first_guess, deg_to_rad(zenith), mirror, theta_star);
        std::cout << "NLSST SST = " << retrieved << " K (" << (retrieved - 273.15) << " degC)\n";
    }

    if (args.has("out")) {
        if (!have_model) {
            std::cerr << "error: nothing to write; calibrate or load a model first\n";
            return 1;
        }
        if (!write_nlsst_model(args.get("out"), model)) {
            std::cerr << "error: could not write " << args.get("out") << "\n";
            return 1;
        }
        std::cout << "wrote model to " << args.get("out") << "\n";
    }

    if (!have_model && !args.has("bt11")) {
        std::cerr << "error: nothing to do; use --file, --coeffs or --bt11\n";
        return 1;
    }
    return 0;
}

int retrieve_ndvi(const Arguments& args) {
    const std::string path = args.get("file");
    if (path.empty()) {
        std::cerr << "error: --file is required\n";
        return 1;
    }
    const Table table = read_csv(path);
    int red = table.column(args.get("red", "red"));
    int nir = table.column(args.get("nir", "nir"));
    if (red < 0 || nir < 0) {
        std::cerr << "error: red/nir columns not found; use --red and --nir\n";
        return 1;
    }

    std::vector<double> results;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        double r = 0.0;
        double n = 0.0;
        if (!parse_number(table.cell(i, red), r) || !parse_number(table.cell(i, nir), n)) continue;
        results.push_back(ndvi(r, n));
        if (std::isfinite(results.back())) ++valid;
    }
    if (results.empty()) {
        std::cerr << "error: no valid reflectance pairs in " << path << "\n";
        return 1;
    }

    double total = 0.0;
    double lowest = std::numeric_limits<double>::infinity();
    double highest = -std::numeric_limits<double>::infinity();
    for (double value : results) {
        if (!std::isfinite(value)) continue;
        total += value;
        lowest = std::min(lowest, value);
        highest = std::max(highest, value);
    }

    std::cout << "NDVI from " << path << "\n";
    std::cout << "  pixels   : " << results.size() << " (" << valid << " valid)\n";
    std::cout << "  mean     : " << (valid > 0 ? total / static_cast<double>(valid) : 0.0) << "\n";
    std::cout << "  minimum  : " << lowest << "\n";
    std::cout << "  maximum  : " << highest << "\n";

    if (args.has("out")) {
        std::ostringstream out;
        out << "index,ndvi\n";
        for (std::size_t i = 0; i < results.size(); ++i)
            out << i << "," << (std::isfinite(results[i]) ? std::to_string(results[i]) : "NA") << "\n";
        write_text_file(args.get("out"), out.str());
        std::cout << "wrote " << args.get("out") << "\n";
    }
    return 0;
}

int retrieve_seaice_grid(const Arguments& args) {
    const std::string path = args.get("file");
    if (path.empty()) {
        std::cerr << "error: --file is required\n";
        return 1;
    }
    std::vector<SeaIceCell> cells;
    std::string error;
    if (!read_sea_ice_grid(path, cells, error)) {
        std::cerr << "error: " << error << "\n";
        return 1;
    }
    const double threshold = args.number("threshold", SIC_EXTENT_THRESHOLD);
    const SeaIceAreaResult result = sea_ice_extent_area(cells, threshold);
    const std::string format = format_name(args);

    if (format == "json") {
        json::Value root = json::Value::object();
        root.set("source", json::Value::from_string(path));
        root.set("threshold", json::Value::from_number(threshold));
        root.set("cells_total", json::Value::from_number(static_cast<double>(result.cells_total)));
        root.set("cells_ice_covered",
                 json::Value::from_number(static_cast<double>(result.cells_ice_covered)));
        root.set("extent_km2", json::Value::from_number(result.extent_km2));
        root.set("extent_million_km2", json::Value::from_number(result.extent_km2 / 1.0e6));
        root.set("area_km2", json::Value::from_number(result.area_km2));
        root.set("area_million_km2", json::Value::from_number(result.area_km2 / 1.0e6));
        root.set("grid_area_million_km2", json::Value::from_number(result.total_grid_km2 / 1.0e6));
        std::cout << json::dump(root, 2) << "\n";
        return 0;
    }

    std::cout << "Sea ice grid " << path << "\n";
    std::cout << "  cells               : " << result.cells_total << "\n";
    std::cout << "  cells >= threshold  : " << result.cells_ice_covered << "\n";
    std::cout << "  threshold           : " << threshold << "\n";
    std::cout << "  extent              : " << result.extent_km2 / 1.0e6 << " million km2\n";
    std::cout << "  area                : " << result.area_km2 / 1.0e6 << " million km2\n";
    std::cout << "  grid area covered   : " << result.total_grid_km2 / 1.0e6 << " million km2\n";
    return 0;
}

int retrieve_seaice_pixel(const Arguments& args) {
    const double tb19v = args.number("tb19v", std::nan(""));
    const double tb19h = args.number("tb19h", std::nan(""));
    const double tb37v = args.number("tb37v", std::nan(""));
    if (!std::isfinite(tb19v) || !std::isfinite(tb19h) || !std::isfinite(tb37v)) {
        std::cerr << "error: --tb19v, --tb19h and --tb37v are required (Kelvin)\n";
        return 1;
    }

    SeaIceTiePoints tie_points;
    if (args.has("tie-points")) {
        std::string error;
        if (!read_sea_ice_tie_points(args.get("tie-points"), tie_points, error)) {
            std::cerr << "error: " << error << "\n";
            return 1;
        }
    }

    const SeaIceConcentrations mixture =
        nasa_team_mixture_invert(tie_points, tb19v, tb19h, tb37v);
    const PolarGradientModel model = calibrate_polar_gradient(tie_points);
    const SeaIceConcentrations empirical =
        nasa_team_polar_gradient(model, tb19v, tb19h, tb37v);

    std::cout << "Passive-microwave footprint\n";
    std::cout << "  PR(19)          : " << polarization_ratio(tb19v, tb19h) << "\n";
    std::cout << "  GR(37V/19V)     : " << gradient_ratio(tb37v, tb19v) << "\n";
    std::cout << "Mixture inversion (physical tie points)\n";
    std::cout << "  open water      : " << mixture.open_water << "\n";
    std::cout << "  first-year ice  : " << mixture.first_year << "\n";
    std::cout << "  multiyear ice   : " << mixture.multiyear << "\n";
    std::cout << "  total ice       : " << mixture.total() << "\n";
    std::cout << "Polarization/gradient retrieval (calibrated)\n";
    std::cout << "  first-year ice  : " << empirical.first_year << "\n";
    std::cout << "  multiyear ice   : " << empirical.multiyear << "\n";
    std::cout << "  total ice       : " << empirical.total() << "\n";
    std::cout << "  extent indicator: " << empirical.extent_indicator() << "\n";
    std::cout << "  model RMSE      : " << model.rmse_cf << " (first-year), " << model.rmse_cm
              << " (multiyear)\n";
    return 0;
}

int retrieve_precip(const Arguments& args) {
    bool printed = false;
    if (args.has("cold-fraction")) {
        const double fraction = args.number("cold-fraction", std::nan(""));
        if (!std::isfinite(fraction)) {
            std::cerr << "error: --cold-fraction must be numeric\n";
            return 1;
        }
        double value = gpi(fraction);
        std::cout << "GPI cold-cloud fraction " << fraction << " -> " << value << " mm/day\n";
        if (args.has("water-vapour")) {
            value = gpi_with_water_vapour(value, args.number("water-vapour", 0.0),
                                          args.number("a1", 0.1), args.number("a2", 3.0));
            std::cout << "water-vapour adjusted -> " << value << " mm/day\n";
        }
        printed = true;
    }

    if (args.has("file")) {
        std::vector<PrecipEstimate> estimates;
        std::string error;
        if (!read_precip_estimates(args.get("file"), estimates, error)) {
            std::cerr << "error: " << error << "\n";
            return 1;
        }
        const MergedPrecip merged = gpcp_merge(estimates);
        std::cout << "GPCP inverse-variance merge from " << args.get("file") << "\n";
        for (std::size_t i = 0; i < estimates.size(); ++i) {
            std::cout << "  " << (estimates[i].source.empty() ? "estimate" : estimates[i].source)
                      << ": value " << estimates[i].value << ", error SD "
                      << estimates[i].error_sd << ", weight " << merged.weights[i] << "\n";
        }
        if (merged.valid)
            std::cout << "  merged = " << merged.value << " mm/day from " << merged.used
                      << " sources\n";
        else
            std::cout << "  no valid estimates\n";
        printed = true;
    }

    if (!printed) {
        std::cerr << "error: provide --file ESTIMATES.csv and/or --cold-fraction F\n";
        return 1;
    }
    return 0;
}

int retrieve_index(const Arguments& args) {
    const std::string path = args.get("file");
    if (path.empty()) {
        std::cerr << "error: --file is required\n";
        return 1;
    }
    const std::string index_name = lowered(args.get("index", "ndvi"));
    if (index_name != "ndvi" && index_name != "evi" && index_name != "savi") {
        std::cerr << "error: --index must be ndvi, evi or savi\n";
        return 1;
    }

    const Table table = read_csv(path);
    const int red = table.column(args.get("red", "red"));
    const int nir = table.column(args.get("nir", "nir"));
    const int blue = table.column(args.get("blue", "blue"));
    if (red < 0 || nir < 0) {
        std::cerr << "error: red/nir columns not found; use --red and --nir\n";
        return 1;
    }
    if (index_name == "evi" && blue < 0) {
        std::cerr << "error: EVI needs a blue reflectance column; use --blue\n";
        return 1;
    }

    const double soil = args.number("soil", 0.5);
    std::vector<double> results;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        double r = 0.0;
        double n = 0.0;
        if (!parse_number(table.cell(i, red), r) || !parse_number(table.cell(i, nir), n)) continue;
        if (index_name == "ndvi") {
            results.push_back(ndvi(r, n));
        } else if (index_name == "savi") {
            results.push_back(savi(r, n, soil));
        } else {
            double b = 0.0;
            if (!parse_number(table.cell(i, blue), b)) continue;
            results.push_back(evi(r, n, b));
        }
    }
    if (results.empty()) {
        std::cerr << "error: no valid reflectance rows in " << path << "\n";
        return 1;
    }

    double total = 0.0;
    double lowest = std::numeric_limits<double>::infinity();
    double highest = -std::numeric_limits<double>::infinity();
    std::size_t valid = 0;
    for (double value : results) {
        if (!std::isfinite(value)) continue;
        total += value;
        ++valid;
        lowest = std::min(lowest, value);
        highest = std::max(highest, value);
    }

    std::cout << index_name << " from " << path << "\n";
    std::cout << "  pixels   : " << results.size() << " (" << valid << " valid)\n";
    std::cout << "  mean     : " << (valid > 0 ? total / static_cast<double>(valid) : 0.0)
              << "\n";
    std::cout << "  minimum  : " << lowest << "\n";
    std::cout << "  maximum  : " << highest << "\n";

    if (args.has("out")) {
        std::ostringstream out;
        out << "index,value\n";
        for (std::size_t i = 0; i < results.size(); ++i)
            out << i << "," << (std::isfinite(results[i]) ? std::to_string(results[i]) : "NA")
                << "\n";
        write_text_file(args.get("out"), out.str());
        std::cout << "wrote " << args.get("out") << "\n";
    }
    return 0;
}

int command_data(const std::string& sub, const Arguments& args) {
    if (sub == "list") {
        const std::string format = format_name(args);
        if (format == "json") {
            json::Value list = json::Value::array();
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
                list.push(item);
            }
            std::cout << json::dump(list, 2) << "\n";
            return 0;
        }
        for (const auto& dataset : remote_datasets()) {
            std::cout << dataset.key << " -> " << dataset.short_name;
            if (!dataset.version.empty()) std::cout << " v" << dataset.version;
            std::cout << " (" << dataset.provider << ")\n";
            std::cout << "  catalog variable : " << dataset.catalog_id << "\n";
            std::cout << "  science variable : " << dataset.variable << " [" << dataset.units
                      << "]\n";
            std::cout << "  conversion mode  : " << dataset.mode << "\n";
            std::cout << "  " << dataset.description << "\n\n";
        }
        return 0;
    }

    if (sub == "inspect") {
        const std::string path = args.get("file");
        if (path.empty()) {
            std::cerr << "error: --file is required\n";
            return 1;
        }
        if (is_hdf5_file(path)) {
            std::cerr << "error: " << path
                      << " is NetCDF-4/HDF5; use the wrapper (data convert) to produce CSV\n";
            return 1;
        }
        const NetCDFFile file = NetCDFFile::open(path);
        const auto print_attributes = [](const std::vector<NcAttribute>& attributes, const std::string& indent) {
            for (const auto& attribute : attributes) {
                std::cout << indent << attribute.name << " = ";
                if (attribute.type == 2) {
                    std::cout << "\"" << attribute.text << "\"";
                } else {
                    const std::size_t limit = std::min<std::size_t>(attribute.numbers.size(), 8);
                    for (std::size_t i = 0; i < limit; ++i) {
                        if (i > 0) std::cout << ", ";
                        std::cout << attribute.numbers[i];
                    }
                    if (attribute.numbers.size() > limit)
                        std::cout << ", ... (" << attribute.numbers.size() << " values)";
                }
                std::cout << "\n";
            }
        };

        std::cout << "NetCDF file : " << path << "\n";
        std::cout << "Format      : classic CDF-" << file.version() << "\n";
        std::cout << "Records     : " << file.record_count() << " (stride "
                  << file.record_stride() << " bytes)\n\nDimensions\n";
        for (const auto& dimension : file.dimensions())
            std::cout << "  " << dimension.name << " = " << dimension.length
                      << (dimension.unlimited ? " (unlimited)" : "") << "\n";
        std::cout << "\nGlobal attributes\n";
        if (file.global_attributes().empty()) std::cout << "  (none)\n";
        print_attributes(file.global_attributes(), "  ");
        std::cout << "\nVariables\n";
        for (const auto& variable : file.variables()) {
            std::cout << "  " << variable.name << ": " << netcdf_type_name(variable.type) << "(";
            for (std::size_t i = 0; i < variable.dimension_ids.size(); ++i) {
                if (i > 0) std::cout << ", ";
                const int id = variable.dimension_ids[i];
                std::cout << file.dimensions()[static_cast<std::size_t>(id)].name;
            }
            std::cout << ")\n";
            print_attributes(variable.attributes, "      ");
        }
        return 0;
    }

    if (sub == "export") {
        const std::string path = args.get("file");
        const std::string variable = args.get("variable");
        if (path.empty() || variable.empty()) {
            std::cerr << "error: --file and --variable are required\n";
            return 1;
        }
        const std::string mode = args.get("mode", "timeseries");
        const std::string output = args.get("out", variable + ".csv");
        const NetCDFExportResult result = export_netcdf_to_csv(path, variable, mode, output);
        if (!result.ok) {
            std::cerr << "error: " << result.error << "\n";
            return 1;
        }
        std::cout << result.summary_json << "\n";
        return 0;
    }

    if (sub != "search" && sub != "fetch" && sub != "convert") {
        std::cerr << "error: unknown data command '" << sub
                  << "'; expected list, search, fetch, convert, inspect or export\n";
        return 1;
    }

    const std::string key = args.get("dataset");
    const RemoteDataset* spec = key.empty() ? nullptr : find_remote_dataset(key);
    if (!key.empty() && spec == nullptr) {
        std::cerr << "error: unknown dataset '" << key << "'; known keys:";
        for (const auto& dataset : remote_datasets()) std::cout << " " << dataset.key;
        std::cerr << "\n";
        return 1;
    }

    std::vector<std::string> command = {resolve_python_interpreter(), resolve_fetch_wrapper(), sub};
    const auto add = [&](const std::string& flag, const std::string& value) {
        if (value.empty()) return;
        command.push_back(flag);
        command.push_back(value);
    };

    if (sub == "convert") {
        const std::string output = args.get("out");
        if (output.empty()) {
            std::cerr << "error: --out CSV is required\n";
            return 1;
        }
        add("--variable", args.get("variable", spec != nullptr ? spec->variable : ""));
        add("--mode", args.get("mode", spec != nullptr ? spec->mode : "timeseries"));
        add("--units", args.get("units", spec != nullptr ? spec->units : ""));
        add("--threshold", args.get("threshold"));
        add("--out", output);
        for (const auto& file : args.positional) command.push_back(file);
        if (args.positional.empty()) {
            std::cerr << "error: provide one or more NetCDF files\n";
            return 1;
        }
    } else {
        std::string short_name = args.get("short-name");
        if (short_name.empty() && spec != nullptr) short_name = spec->short_name;
        if (short_name.empty()) {
            std::cerr << "error: --dataset KEY or --short-name is required\n";
            return 1;
        }
        add("--short-name", short_name);
        add("--version", args.get("version", spec != nullptr ? spec->version : ""));
        add("--concept-id", args.get("concept-id"));
        add("--start", args.get("start"));
        add("--end", args.get("end"));
        add("--bbox", args.get("bbox"));
        add("--max", args.get("max"));
        const std::string out_dir = args.get("out", "data/nasa");
        if (sub == "fetch") {
            std::string csv = args.get("out-csv");
            if (csv.empty() && spec != nullptr) {
                csv = out_dir + "/" + spec->key +
                      (spec->mode == "grid" ? "_grid.csv" : "_monthly.csv");
            }
            if (csv.empty()) {
                std::cerr << "error: --out-csv is required when --short-name is used directly\n";
                return 1;
            }
            add("--out", out_dir);
            add("--out-csv", csv);
            add("--variable", args.get("variable", spec != nullptr ? spec->variable : ""));
            add("--mode", args.get("mode", spec != nullptr ? spec->mode : "timeseries"));
            add("--units", args.get("units", spec != nullptr ? spec->units : ""));
            add("--threshold", args.get("threshold"));
            if (args.flag("dry-run")) command.push_back("--dry-run");
        }
    }

    const std::string invocation = build_command(command);
    if (args.flag("print-command")) std::cout << invocation << "\n";
    const CommandResult result = run_command(invocation);
    std::cout << result.output;
    if (result.exit_code != 0) {
        std::cerr << "error: data wrapper exited with code " << result.exit_code << "\n";
        return result.exit_code == 0 ? 1 : result.exit_code;
    }
    return 0;
}

int command_retrieve(const std::string& sub, const Arguments& args) {
    if (sub == "sst") return retrieve_sst(args);
    if (sub == "ndvi") return retrieve_ndvi(args);
    if (sub == "seaice") return retrieve_seaice_grid(args);
    if (sub == "seaice-pixel") return retrieve_seaice_pixel(args);
    if (sub == "precip") return retrieve_precip(args);
    if (sub == "index") return retrieve_index(args);
    std::cerr << "error: unknown retrieval '" << sub
              << "'; expected sst, ndvi, index, seaice, seaice-pixel or precip\n";
    return 1;
}

int command_calibrate(const std::string& sub, const Arguments& args) {
    if (sub == "nlsst") {
        const std::string path = args.get("file");
        if (path.empty()) {
            std::cerr << "error: --file MATCHUPS.csv is required\n";
            return 1;
        }
        std::vector<NLSSTSample> samples;
        std::string error;
        if (!read_nlsst_samples(path, samples, error)) {
            std::cerr << "error: " << error << "\n";
            return 1;
        }
        const NLSSTModel model = calibrate_nlsst(samples);
        const std::string out = args.get("out", "nlsst_model.csv");
        if (!write_nlsst_model(out, model)) {
            std::cerr << "error: could not write " << out << "\n";
            return 1;
        }
        std::cout << "calibrated NLSST from " << model.samples << " matchups\n";
        std::cout << "  RMSE = " << model.rmse << " K, bias = " << model.bias << " K\n";
        std::cout << "  wrote " << out << "\n";
        return 0;
    }

    if (sub == "seaice") {
        SeaIceTiePoints tie_points;
        if (args.has("tie-points")) {
            std::string error;
            if (!read_sea_ice_tie_points(args.get("tie-points"), tie_points, error)) {
                std::cerr << "error: " << error << "\n";
                return 1;
            }
        }
        const PolarGradientModel model = calibrate_polar_gradient(tie_points);
        const std::string out = args.get("out", "nasa_team_model.csv");
        if (!write_polar_gradient_model(out, model)) {
            std::cerr << "error: could not write " << out << "\n";
            return 1;
        }
        std::cout << "calibrated NASA Team polarization/gradient model from " << model.samples
                  << " mixture samples\n";
        std::cout << "  RMSE (first-year) = " << model.rmse_cf << "\n";
        std::cout << "  RMSE (multiyear)  = " << model.rmse_cm << "\n";
        std::cout << "  wrote " << out << "\n";
        return 0;
    }

    std::cerr << "error: unknown calibration '" << sub << "'; expected nlsst or seaice\n";
    return 1;
}

int command_physics(const std::string& sub, const Arguments& args) {
    if (sub == "toa") {
        const double forcing = args.number("forcing", std::nan(""));
        const double feedback = args.number("feedback", std::nan(""));
        const double delta_t = args.number("delta-t", std::nan(""));
        if (!std::isfinite(forcing) || !std::isfinite(feedback) || !std::isfinite(delta_t)) {
            std::cerr << "error: --forcing, --feedback and --delta-t are required\n";
            return 1;
        }
        std::cout << "TOA energy imbalance N = F + lambda dT = "
                  << toa_energy_imbalance(forcing, feedback, delta_t) << " W/m2\n";
        return 0;
    }

    if (sub == "surface") {
        const double albedo = args.number("albedo", std::nan(""));
        const double shortwave = args.number("shortwave", std::nan(""));
        const double emissivity_air = args.number("emissivity-air", std::nan(""));
        const double longwave = args.number("longwave", std::nan(""));
        const double emissivity_surface = args.number("emissivity-surface", std::nan(""));
        const double temperature = args.number("temperature", std::nan(""));
        if (!std::isfinite(albedo) || !std::isfinite(shortwave) || !std::isfinite(emissivity_air) ||
            !std::isfinite(longwave) || !std::isfinite(emissivity_surface) ||
            !std::isfinite(temperature)) {
            std::cerr << "error: --albedo, --shortwave, --emissivity-air, --longwave, "
                         "--emissivity-surface and --temperature are required\n";
            return 1;
        }
        const double net =
            net_surface_radiation(albedo, shortwave, emissivity_air, longwave,
                                  emissivity_surface, temperature);
        std::cout << "net surface radiation Rn = " << net << " W/m2\n";
        if (args.has("latent") && args.has("sensible") && args.has("ground")) {
            const double residual =
                surface_energy_residual(net, args.number("latent", 0.0),
                                        args.number("sensible", 0.0), args.number("ground", 0.0));
            std::cout << "energy balance residual Rn - (LE + H + G) = " << residual << " W/m2\n";
        }
        return 0;
    }

    if (sub == "carbon") {
        const double start_ppm = args.number("start-ppm", std::nan(""));
        const double end_ppm = args.number("end-ppm", std::nan(""));
        const double start_year = args.number("start-year", std::nan(""));
        const double end_year = args.number("end-year", std::nan(""));
        if (!std::isfinite(start_ppm) || !std::isfinite(end_ppm) || !std::isfinite(start_year) ||
            !std::isfinite(end_year)) {
            std::cerr << "error: --start-ppm, --end-ppm, --start-year and --end-year are required\n";
            return 1;
        }
        std::cout << "atmospheric carbon growth = "
                  << carbon_growth_rate_pgc_per_year(start_ppm, end_ppm, start_year, end_year)
                  << " PgC/year\n";
        return 0;
    }

    if (sub == "amoc") {
        const std::string path = args.get("file");
        if (path.empty()) {
            std::cerr << "error: --file VELOCITY.csv is required\n";
            return 1;
        }
        std::vector<VelocitySample> samples;
        std::string error;
        if (!read_velocity_samples(path, samples, error)) {
            std::cerr << "error: " << error << "\n";
            return 1;
        }
        const double depth = args.number("depth", 0.0);
        std::cout << "AMOC streamfunction Q(z = " << depth << " m) = "
                  << amoc_streamfunction_sv(samples, depth) << " Sv\n";
        return 0;
    }

    if (sub == "feedback") {
        FeedbackSet feedback;
        feedback.planck = args.number("planck", feedback.planck);
        feedback.water_vapour = args.number("water-vapour", feedback.water_vapour);
        feedback.lapse_rate = args.number("lapse-rate", feedback.lapse_rate);
        feedback.albedo = args.number("albedo", feedback.albedo);
        feedback.cloud = args.number("cloud", feedback.cloud);
        const double forcing = args.number("forcing", 3.7);
        const double total = feedback.total();
        std::cout << "feedback components (W/m2/K): Planck " << feedback.planck
                  << ", water vapour " << feedback.water_vapour << ", lapse rate "
                  << feedback.lapse_rate << ", albedo " << feedback.albedo << ", cloud "
                  << feedback.cloud << "\n";
        std::cout << "net feedback lambda = " << total << " W/m2/K\n";
        if (total != 0.0) {
            std::cout << "equilibrium warming for F = " << forcing
                      << " W/m2: " << equilibrium_temperature_change(forcing, total) << " K\n";
            std::cout << "ECS (2xCO2): " << climate_sensitivity_2x(total) << " K\n";
        }
        return 0;
    }

    if (sub == "forcing-co2") {
        const double ppm = args.number("ppm", std::nan(""));
        const double reference = args.number("reference-ppm", 280.0);
        if (!std::isfinite(ppm) || ppm <= 0.0 || reference <= 0.0) {
            std::cerr << "error: --ppm (positive) is required; --reference-ppm defaults to 280\n";
            return 1;
        }
        std::cout << "CO2 radiative forcing " << ppm << " vs " << reference << " ppm = "
                  << co2_radiative_forcing(ppm, reference) << " W/m2\n";
        return 0;
    }

    if (sub == "mie") {
        const double size_parameter = args.number("size-parameter", std::nan(""));
        const double real = args.number("real", 1.33);
        const double imag = args.number("imag", 0.0);
        if (!std::isfinite(size_parameter) || size_parameter <= 0.0) {
            std::cerr << "error: --size-parameter (2 pi r / lambda) is required\n";
            return 1;
        }
        const MieResult result =
            mie_sphere(size_parameter, std::complex<double>(real, imag));
        std::cout << "Mie sphere, x = " << size_parameter << ", m = " << real;
        if (imag != 0.0) std::cout << " + " << imag << "i";
        std::cout << "\n  Qext = " << result.qext << "\n  Qsca = " << result.qsca
                  << "\n  asymmetry g = " << result.g
                  << "\n  single-scattering albedo = " << result.single_scattering_albedo
                  << "\n";
        return 0;
    }

    if (sub == "aod") {
        const double number_concentration = args.number("number-concentration", std::nan(""));
        const double radius_um = args.number("radius-um", std::nan(""));
        const double path_length = args.number("path-length", std::nan(""));
        const double qext = args.number("qext", 2.0);
        if (!std::isfinite(number_concentration) || !std::isfinite(radius_um) ||
            !std::isfinite(path_length)) {
            std::cerr << "error: --number-concentration, --radius-um and --path-length are required\n";
            return 1;
        }
        std::cout << "aerosol optical depth = "
                  << aerosol_optical_depth(number_concentration, radius_um * 1e-6, qext,
                                           path_length)
                  << "\n";
        return 0;
    }

    if (sub == "wind-stress") {
        const double u10 = args.number("u10", std::nan(""));
        const double v10 = args.number("v10", std::nan(""));
        if (!std::isfinite(u10) || !std::isfinite(v10)) {
            std::cerr << "error: --u10 and --v10 are required (m/s)\n";
            return 1;
        }
        const WindStress stress = bulk_wind_stress(args.number("air-density", 1.225),
                                                   args.number("drag", 1.3e-3), u10, v10);
        std::cout << "wind stress tau_x = " << stress.tau_x << " N/m2, tau_y = " << stress.tau_y
                  << " N/m2, magnitude = " << stress.magnitude() << " N/m2\n";
        return 0;
    }

    if (sub == "ekman") {
        const double latitude = args.number("latitude", std::nan(""));
        const double tau_x = args.number("tau-x", std::nan(""));
        const double tau_y = args.number("tau-y", std::nan(""));
        if (!std::isfinite(latitude) || !std::isfinite(tau_x) || !std::isfinite(tau_y)) {
            std::cerr << "error: --tau-x, --tau-y and --latitude are required\n";
            return 1;
        }
        const double f = coriolis_parameter(latitude);
        if (f == 0.0) {
            std::cerr << "error: Ekman transport is undefined at the equator\n";
            return 1;
        }
        const EkmanTransport transport =
            ekman_transport(tau_x, tau_y, args.number("density", 1027.0), f);
        std::cout << "Coriolis parameter f = " << f << " rad/s\n";
        std::cout << "Ekman transport Mx = " << transport.m_x << " m2/s, My = " << transport.m_y
                  << " m2/s, magnitude = " << transport.magnitude() << " m2/s\n";
        return 0;
    }

    if (sub == "geostrophic") {
        const double latitude = args.number("latitude", std::nan(""));
        const double density = args.number("density", std::nan(""));
        const double dp_dx = args.number("dp-dx", std::nan(""));
        const double dp_dy = args.number("dp-dy", std::nan(""));
        if (!std::isfinite(latitude) || !std::isfinite(density) || !std::isfinite(dp_dx) ||
            !std::isfinite(dp_dy)) {
            std::cerr << "error: --density, --dp-dx, --dp-dy and --latitude are required\n";
            return 1;
        }
        const GeostrophicVelocity velocity =
            geostrophic_velocity(density, dp_dx, dp_dy, latitude);
        std::cout << "geostrophic velocity u = " << velocity.u << " m/s, v = " << velocity.v
                  << " m/s\n";
        return 0;
    }

    if (sub == "thermal-wind") {
        const double latitude = args.number("latitude", std::nan(""));
        const double density = args.number("density", std::nan(""));
        const double drho_dx = args.number("drho-dx", std::nan(""));
        const double drho_dy = args.number("drho-dy", std::nan(""));
        if (!std::isfinite(latitude) || !std::isfinite(density) || !std::isfinite(drho_dx) ||
            !std::isfinite(drho_dy)) {
            std::cerr << "error: --density, --drho-dx, --drho-dy and --latitude are required\n";
            return 1;
        }
        const ThermalWind shear = thermal_wind(density, drho_dx, drho_dy, latitude);
        std::cout << "thermal wind shear du/dz = " << shear.du_dz << " 1/s, dv/dz = "
                  << shear.dv_dz << " 1/s\n";
        return 0;
    }

    std::cerr << "error: unknown physics diagnostic '" << sub
              << "'; expected toa, surface, carbon, amoc, feedback, forcing-co2, mie, aod, "
                 "wind-stress, ekman, geostrophic or thermal-wind\n";
    return 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    try {
        const std::string command = argv[1];
        if (command == "help" || command == "--help" || command == "-h") {
            print_usage();
            return 0;
        }

        if (command == "catalog") {
            const Arguments args = parse_arguments(argc, argv, 2);
            return command_catalog(args);
        }
        if (command == "query") {
            const Arguments args = parse_arguments(argc, argv, 2);
            return command_query(args);
        }
        if (command == "intent") {
            const Arguments args = parse_arguments(argc, argv, 2);
            return command_intent(args);
        }
        if (command == "analyze") {
            const Arguments args = parse_arguments(argc, argv, 2);
            return command_analyze(args);
        }
        if (command == "co2") {
            const Arguments args = parse_arguments(argc, argv, 2);
            return command_co2(args);
        }
        if (command == "data") {
            if (argc < 3) {
                std::cerr << "error: data needs a subcommand; run climate help\n";
                return 1;
            }
            const std::string sub = argv[2];
            const Arguments args = parse_arguments(argc, argv, 3);
            return command_data(sub, args);
        }
        if (command == "retrieve" || command == "calibrate" || command == "physics") {
            if (argc < 3) {
                std::cerr << "error: " << command << " needs a subcommand; run climate help\n";
                return 1;
            }
            const std::string sub = argv[2];
            const Arguments args = parse_arguments(argc, argv, 3);
            if (command == "retrieve") return command_retrieve(sub, args);
            if (command == "calibrate") return command_calibrate(sub, args);
            return command_physics(sub, args);
        }

        std::cerr << "error: unknown command '" << command << "'; run climate help\n";
        return 1;
    } catch (const std::exception& exception) {
        std::cerr << "error: " << exception.what() << "\n";
        return 1;
    }
}
