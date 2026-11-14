#include "retrieval.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>

#include "common.hpp"
#include "csv.hpp"
#include "linalg.hpp"
#include "stats.hpp"

namespace climate {
namespace {

double safe_divide(double numerator, double denominator) {
    if (denominator == 0.0 || !std::isfinite(denominator)) return std::nan("");
    return numerator / denominator;
}

int find_column(const Table& table, const std::vector<std::string>& names) {
    for (const auto& name : names) {
        const int index = table.column(name);
        if (index >= 0) return index;
    }
    return -1;
}

}  // namespace

std::vector<double> nlsst_features(double bt11_k, double bt12_k, double first_guess_k,
                                   double sat_zenith_rad, double mirror, double theta_star_rad) {
    const double difference = bt11_k - bt12_k;
    const double cosine = std::cos(sat_zenith_rad);
    if (std::fabs(cosine) < 1e-6)
        throw std::runtime_error("nlsst_features: satellite zenith angle is too close to 90 degrees");
    const double secant = 1.0 / cosine;

    return {
        1.0,
        bt11_k,
        difference * first_guess_k,
        (secant - 1.0) * difference,
        mirror,
        theta_star_rad,
        sat_zenith_rad * sat_zenith_rad,
    };
}

NLSSTModel calibrate_nlsst(const std::vector<NLSSTSample>& samples) {
    if (samples.size() < 8)
        throw std::runtime_error("calibrate_nlsst: at least 8 matchups are required");

    std::vector<std::vector<double>> design;
    std::vector<double> response;
    design.reserve(samples.size());
    response.reserve(samples.size());

    for (const auto& sample : samples) {
        if (!std::isfinite(sample.buoy_sst_k) || !std::isfinite(sample.bt11_k) ||
            !std::isfinite(sample.bt12_k) || !std::isfinite(sample.first_guess_k))
            continue;
        design.push_back(nlsst_features(sample.bt11_k, sample.bt12_k, sample.first_guess_k,
                                        sample.sat_zenith_rad, sample.mirror,
                                        sample.theta_star_rad));
        response.push_back(sample.buoy_sst_k);
    }

    if (design.size() < 8)
        throw std::runtime_error("calibrate_nlsst: too few valid matchups after filtering");

    NLSSTModel model;
    if (!least_squares(design, response, model.coefficients))
        throw std::runtime_error("calibrate_nlsst: rank-deficient matchup set");

    double squared = 0.0;
    double bias = 0.0;
    for (std::size_t i = 0; i < design.size(); ++i) {
        double fitted = 0.0;
        for (std::size_t j = 0; j < model.coefficients.size(); ++j)
            fitted += model.coefficients[j] * design[i][j];
        const double residual = response[i] - fitted;
        squared += residual * residual;
        bias += residual;
    }
    model.rmse = std::sqrt(squared / static_cast<double>(design.size()));
    model.bias = bias / static_cast<double>(design.size());
    model.samples = design.size();
    model.calibrated = true;
    return model;
}

double nlsst_retrieve(const NLSSTModel& model, double bt11_k, double bt12_k, double first_guess_k,
                      double sat_zenith_rad, double mirror, double theta_star_rad) {
    if (model.coefficients.size() != 7)
        throw std::runtime_error("nlsst_retrieve: model does not contain 7 coefficients");
    const std::vector<double> features =
        nlsst_features(bt11_k, bt12_k, first_guess_k, sat_zenith_rad, mirror, theta_star_rad);
    return dot_product(model.coefficients, features);
}

bool read_nlsst_samples(const std::string& path, std::vector<NLSSTSample>& out, std::string& error) {
    try {
        const Table table = read_csv(path);
        const int bt11 = find_column(table, {"bt11", "bt_11", "t11", "brightness_temperature_11"});
        const int bt12 = find_column(table, {"bt12", "bt_12", "t12", "brightness_temperature_12"});
        const int guess = find_column(table, {"first_guess", "tsfc", "t_sfc", "climatology"});
        const int zenith = find_column(table, {"sat_zenith_deg", "zenith_deg", "satellite_zenith_deg", "theta_deg"});
        const int mirror = find_column(table, {"mirror", "mirror_side"});
        const int theta_star = find_column(table, {"theta_star_deg", "relative_azimuth_deg", "azimuth_deg"});
        const int buoy = find_column(table, {"buoy_sst", "buoy_sst_k", "truth_sst", "sst_buoy"});
        if (bt11 < 0 || bt12 < 0 || guess < 0 || zenith < 0 || buoy < 0)
            throw std::runtime_error(
                "required columns are bt11, bt12, first_guess, sat_zenith_deg and buoy_sst");

        out.clear();
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            NLSSTSample sample;
            double value = 0.0;
            if (!parse_number(table.cell(i, bt11), value)) continue;
            sample.bt11_k = value;
            if (!parse_number(table.cell(i, bt12), value)) continue;
            sample.bt12_k = value;
            if (!parse_number(table.cell(i, guess), value)) continue;
            sample.first_guess_k = value;
            double zenith_deg = 0.0;
            if (!parse_number(table.cell(i, zenith), zenith_deg)) continue;
            sample.sat_zenith_rad = deg_to_rad(zenith_deg);
            if (mirror >= 0 && parse_number(table.cell(i, mirror), value))
                sample.mirror = value;
            if (theta_star >= 0 && parse_number(table.cell(i, theta_star), value))
                sample.theta_star_rad = deg_to_rad(value);
            if (!parse_number(table.cell(i, buoy), value)) continue;
            sample.buoy_sst_k = value;
            out.push_back(sample);
        }
        if (out.empty()) throw std::runtime_error("no valid matchup rows found");
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

bool write_nlsst_model(const std::string& path, const NLSSTModel& model) {
    try {
        std::string content = "coefficient,value\n";
        const char* names[7] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};
        for (std::size_t i = 0; i < model.coefficients.size() && i < 7; ++i) {
            content += names[i];
            content += ",";
            content += std::to_string(model.coefficients[i]);
            content += "\n";
        }
        content += "rmse," + std::to_string(model.rmse) + "\n";
        content += "bias," + std::to_string(model.bias) + "\n";
        content += "samples," + std::to_string(model.samples) + "\n";
        write_text_file(path, content);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool read_nlsst_model(const std::string& path, NLSSTModel& model, std::string& error) {
    try {
        const Table table = read_csv(path);
        const int name_column = find_column(table, {"coefficient", "name", "parameter"});
        const int value_column = find_column(table, {"value", "coefficient_value"});
        if (name_column < 0 || value_column < 0)
            throw std::runtime_error("coefficient file needs 'coefficient' and 'value' columns");

        std::map<std::string, double> values;
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            double value = 0.0;
            if (!parse_number(table.cell(i, value_column), value)) continue;
            values[trim_copy(table.cell(i, name_column))] = value;
        }

        model.coefficients.clear();
        for (int i = 0; i < 7; ++i) {
            const std::string key = "a" + std::to_string(i);
            const auto found = values.find(key);
            if (found == values.end())
                throw std::runtime_error("coefficient file is missing '" + key + "'");
            model.coefficients.push_back(found->second);
        }
        if (values.count("rmse")) model.rmse = values["rmse"];
        if (values.count("bias")) model.bias = values["bias"];
        if (values.count("samples")) model.samples = static_cast<std::size_t>(values["samples"]);
        model.calibrated = true;
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

double ndvi(double red_reflectance, double nir_reflectance) {
    const double denominator = nir_reflectance + red_reflectance;
    return safe_divide(nir_reflectance - red_reflectance, denominator);
}

double polarization_ratio(double tb19v, double tb19h) {
    return safe_divide(tb19v - tb19h, tb19v + tb19h);
}

double gradient_ratio(double tb37v, double tb19v) {
    return safe_divide(tb37v - tb19v, tb37v + tb19v);
}

SeaIceConcentrations nasa_team_mixture_invert(const SeaIceTiePoints& tie_points, double tb19v,
                                              double tb19h, double tb37v) {
    // Three channel equations plus the closure C_OW + C_FY + C_MY = 1.
    std::vector<std::vector<double>> matrix(4, std::vector<double>(4, 1.0));
    for (int channel = 0; channel < 3; ++channel)
        for (int constituent = 0; constituent < 3; ++constituent)
            matrix[static_cast<std::size_t>(channel)][static_cast<std::size_t>(constituent)] =
                tie_points.tb[static_cast<std::size_t>(constituent)][static_cast<std::size_t>(channel)];

    std::vector<double> rhs = {tb19v, tb19h, tb37v, 1.0};
    std::vector<double> solution;
    if (!solve_linear_system(matrix, rhs, solution))
        throw std::runtime_error("nasa_team_mixture_invert: singular tie-point matrix");

    SeaIceConcentrations result;
    result.open_water = clamp_value(solution[0], 0.0, 1.0);
    result.first_year = clamp_value(solution[1], 0.0, 1.0);
    result.multiyear = clamp_value(solution[2], 0.0, 1.0);

    const double sum = result.open_water + result.first_year + result.multiyear;
    if (sum > 1.0) {
        result.open_water /= sum;
        result.first_year /= sum;
        result.multiyear /= sum;
    }
    return result;
}

PolarGradientModel calibrate_polar_gradient(const SeaIceTiePoints& tie_points, int steps) {
    if (steps < 2) throw std::runtime_error("calibrate_polar_gradient: steps must be >= 2");

    std::vector<std::vector<double>> design_cf;
    std::vector<double> response_cf;
    std::vector<std::vector<double>> design_cm;
    std::vector<double> response_cm;

    for (int i = 0; i <= steps; ++i) {
        for (int j = 0; j <= steps - i; ++j) {
            const int k = steps - i - j;
            const double ow = static_cast<double>(i) / static_cast<double>(steps);
            const double fy = static_cast<double>(j) / static_cast<double>(steps);
            const double my = static_cast<double>(k) / static_cast<double>(steps);

            double tb[3] = {0.0, 0.0, 0.0};
            for (int channel = 0; channel < 3; ++channel) {
                tb[channel] = ow * tie_points.tb[0][static_cast<std::size_t>(channel)] +
                              fy * tie_points.tb[1][static_cast<std::size_t>(channel)] +
                              my * tie_points.tb[2][static_cast<std::size_t>(channel)];
            }
            const double pr = polarization_ratio(tb[0], tb[1]);
            const double gr = gradient_ratio(tb[2], tb[0]);
            if (!std::isfinite(pr) || !std::isfinite(gr)) continue;

            // CF (1 + d1 PR + d2 GR + d3 PR GR) = a0 + a1 PR + a2 GR + a3 PR GR
            // rearranges to a linear regression of CF on the seven features below.
            design_cf.push_back({1.0, pr, gr, pr * gr, -fy * pr, -fy * gr, -fy * pr * gr});
            response_cf.push_back(fy);
            design_cm.push_back({1.0, pr, gr, pr * gr, -my * pr, -my * gr, -my * pr * gr});
            response_cm.push_back(my);
        }
    }

    if (design_cf.size() < 8)
        throw std::runtime_error("calibrate_polar_gradient: not enough mixture samples");

    std::vector<double> cf_solution;
    std::vector<double> cm_solution;
    if (!least_squares(design_cf, response_cf, cf_solution) ||
        !least_squares(design_cm, response_cm, cm_solution))
        throw std::runtime_error("calibrate_polar_gradient: rank-deficient mixture set");

    PolarGradientModel model;
    model.cf.assign(cf_solution.begin(), cf_solution.begin() + 4);
    model.cm.assign(cm_solution.begin(), cm_solution.begin() + 4);
    model.d = {1.0, cf_solution[4], cf_solution[5], cf_solution[6]};
    model.samples = design_cf.size();

    double squared_cf = 0.0;
    double squared_cm = 0.0;
    for (std::size_t s = 0; s < design_cf.size(); ++s) {
        double numerator_cf = 0.0;
        double numerator_cm = 0.0;
        for (std::size_t j = 0; j < 4; ++j) {
            numerator_cf += model.cf[j] * design_cf[s][j];
            numerator_cm += model.cm[j] * design_cf[s][j];
        }
        double denominator = model.d[0];
        for (std::size_t j = 1; j < 4; ++j) denominator += model.d[j] * design_cf[s][j];
        if (std::fabs(denominator) < 1e-9) continue;
        const double residual_cf = response_cf[s] - numerator_cf / denominator;
        const double residual_cm = response_cm[s] - numerator_cm / denominator;
        squared_cf += residual_cf * residual_cf;
        squared_cm += residual_cm * residual_cm;
    }
    model.rmse_cf = std::sqrt(squared_cf / static_cast<double>(design_cf.size()));
    model.rmse_cm = std::sqrt(squared_cm / static_cast<double>(design_cf.size()));
    model.calibrated = true;
    return model;
}

SeaIceConcentrations nasa_team_polar_gradient(const PolarGradientModel& model, double tb19v,
                                              double tb19h, double tb37v) {
    if (!model.calibrated || model.cf.size() != 4 || model.cm.size() != 4 || model.d.size() != 4)
        throw std::runtime_error("nasa_team_polar_gradient: model is not calibrated");

    const double pr = polarization_ratio(tb19v, tb19h);
    const double gr = gradient_ratio(tb37v, tb19v);
    if (!std::isfinite(pr) || !std::isfinite(gr)) return SeaIceConcentrations{};

    const std::vector<double> features = {1.0, pr, gr, pr * gr};
    const double numerator_cf = dot_product(model.cf, features);
    const double numerator_cm = dot_product(model.cm, features);
    const double denominator =
        model.d[0] + model.d[1] * pr + model.d[2] * gr + model.d[3] * pr * gr;
    if (std::fabs(denominator) < 1e-9) return SeaIceConcentrations{};

    SeaIceConcentrations result;
    result.first_year = clamp_value(numerator_cf / denominator, 0.0, 1.0);
    result.multiyear = clamp_value(numerator_cm / denominator, 0.0, 1.0);
    if (result.total() > 1.0) {
        const double scale = 1.0 / result.total();
        result.first_year *= scale;
        result.multiyear *= scale;
    }
    result.open_water = 1.0 - result.total();
    return result;
}

bool read_sea_ice_tie_points(const std::string& path, SeaIceTiePoints& tie_points,
                             std::string& error) {
    try {
        const Table table = read_csv(path);
        const int constituent = find_column(table, {"constituent", "type", "class"});
        const int channel19v = find_column(table, {"channel_19v", "tb19v", "19v"});
        const int channel19h = find_column(table, {"channel_19h", "tb19h", "19h"});
        const int channel37v = find_column(table, {"channel_37v", "tb37v", "37v"});
        if (constituent < 0 || channel19v < 0 || channel19h < 0 || channel37v < 0)
            throw std::runtime_error(
                "tie-point file needs constituent, channel_19v, channel_19h and channel_37v");

        bool found[3] = {false, false, false};
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            std::string name = trim_copy(table.cell(i, constituent));
            for (auto& c : name) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            int index = -1;
            if (name == "open_water" || name == "openwater" || name == "water")
                index = 0;
            else if (name == "first_year" || name == "firstyear" || name == "fy")
                index = 1;
            else if (name == "multiyear" || name == "multi_year" || name == "my")
                index = 2;
            if (index < 0) continue;

            double v19 = 0.0;
            double h19 = 0.0;
            double v37 = 0.0;
            if (!parse_number(table.cell(i, channel19v), v19) ||
                !parse_number(table.cell(i, channel19h), h19) ||
                !parse_number(table.cell(i, channel37v), v37))
                throw std::runtime_error("invalid brightness temperature in tie-point file");
            tie_points.tb[static_cast<std::size_t>(index)][0] = v19;
            tie_points.tb[static_cast<std::size_t>(index)][1] = h19;
            tie_points.tb[static_cast<std::size_t>(index)][2] = v37;
            found[index] = true;
        }
        if (!found[0] || !found[1] || !found[2])
            throw std::runtime_error(
                "tie-point file must define open_water, first_year and multiyear rows");
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

bool write_polar_gradient_model(const std::string& path, const PolarGradientModel& model) {
    try {
        std::string content = "term,cf,cm,d\n";
        const char* names[4] = {"intercept", "pr", "gr", "pr_gr"};
        for (std::size_t i = 0; i < 4; ++i) {
            content += names[i];
            content += ",";
            content += std::to_string(model.cf[i]);
            content += ",";
            content += std::to_string(model.cm[i]);
            content += ",";
            content += std::to_string(model.d[i]);
            content += "\n";
        }
        content += "rmse_cf," + std::to_string(model.rmse_cf) + ",,\n";
        content += "rmse_cm," + std::to_string(model.rmse_cm) + ",,\n";
        content += "samples," + std::to_string(model.samples) + ",,\n";
        write_text_file(path, content);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool read_polar_gradient_model(const std::string& path, PolarGradientModel& model,
                               std::string& error) {
    try {
        const Table table = read_csv(path);
        const int term = find_column(table, {"term", "parameter", "coefficient"});
        const int cf = find_column(table, {"cf"});
        const int cm = find_column(table, {"cm"});
        const int d = find_column(table, {"d"});
        if (term < 0 || cf < 0 || cm < 0 || d < 0)
            throw std::runtime_error("model file needs term, cf, cm and d columns");

        std::vector<double> cf_values(4, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> cm_values(4, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> d_values(4, std::numeric_limits<double>::quiet_NaN());

        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            const std::string name = trim_copy(table.cell(i, term));
            int index = -1;
            if (name == "intercept")
                index = 0;
            else if (name == "pr")
                index = 1;
            else if (name == "gr")
                index = 2;
            else if (name == "pr_gr")
                index = 3;
            if (index < 0) continue;
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
            if (!parse_number(table.cell(i, cf), a) || !parse_number(table.cell(i, cm), b) ||
                !parse_number(table.cell(i, d), c))
                throw std::runtime_error("invalid model coefficient row");
            cf_values[static_cast<std::size_t>(index)] = a;
            cm_values[static_cast<std::size_t>(index)] = b;
            d_values[static_cast<std::size_t>(index)] = c;
        }
        for (int i = 0; i < 4; ++i) {
            if (!std::isfinite(cf_values[static_cast<std::size_t>(i)]) ||
                !std::isfinite(cm_values[static_cast<std::size_t>(i)]) ||
                !std::isfinite(d_values[static_cast<std::size_t>(i)]))
                throw std::runtime_error("model file is missing a required coefficient");
        }
        model.cf = cf_values;
        model.cm = cm_values;
        model.d = d_values;
        model.calibrated = true;
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

double spherical_cell_area_km2(double lat_min_deg, double lat_max_deg, double lon_min_deg,
                               double lon_max_deg) {
    const double lat_min = deg_to_rad(std::min(lat_min_deg, lat_max_deg));
    const double lat_max = deg_to_rad(std::max(lat_min_deg, lat_max_deg));
    double delta_lon = std::fabs(lon_max_deg - lon_min_deg);
    if (delta_lon > 360.0) delta_lon = 360.0;
    const double delta_lon_rad = deg_to_rad(delta_lon);
    const double area = EARTH_RADIUS_KM * EARTH_RADIUS_KM * delta_lon_rad *
                        (std::sin(lat_max) - std::sin(lat_min));
    return std::fabs(area);
}

SeaIceAreaResult sea_ice_extent_area(const std::vector<SeaIceCell>& cells, double threshold) {
    SeaIceAreaResult result;
    result.cells_total = cells.size();
    for (const auto& cell : cells) {
        const double area = spherical_cell_area_km2(cell.lat_min_deg, cell.lat_max_deg,
                                                    cell.lon_min_deg, cell.lon_max_deg);
        result.total_grid_km2 += area;
        const double concentration = clamp_value(cell.concentration, 0.0, 1.0);
        if (concentration >= threshold) {
            result.extent_km2 += area;
            ++result.cells_ice_covered;
        }
        if (concentration > 0.0) result.area_km2 += concentration * area;
    }
    return result;
}

bool read_sea_ice_grid(const std::string& path, std::vector<SeaIceCell>& cells,
                       std::string& error) {
    try {
        const Table table = read_csv(path);
        const int lat_min = find_column(table, {"lat_min", "lat_min_deg", "latitude_min"});
        const int lat_max = find_column(table, {"lat_max", "lat_max_deg", "latitude_max"});
        const int lon_min = find_column(table, {"lon_min", "lon_min_deg", "longitude_min"});
        const int lon_max = find_column(table, {"lon_max", "lon_max_deg", "longitude_max"});
        const int concentration =
            find_column(table, {"concentration", "sic", "ice_concentration", "sea_ice_concentration"});
        if (lat_min < 0 || lat_max < 0 || lon_min < 0 || lon_max < 0 || concentration < 0)
            throw std::runtime_error(
                "grid file needs lat_min, lat_max, lon_min, lon_max and concentration columns");

        cells.clear();
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            SeaIceCell cell;
            double value = 0.0;
            if (!parse_number(table.cell(i, lat_min), value)) continue;
            cell.lat_min_deg = value;
            if (!parse_number(table.cell(i, lat_max), value)) continue;
            cell.lat_max_deg = value;
            if (!parse_number(table.cell(i, lon_min), value)) continue;
            cell.lon_min_deg = value;
            if (!parse_number(table.cell(i, lon_max), value)) continue;
            cell.lon_max_deg = value;
            if (!parse_number(table.cell(i, concentration), value)) continue;
            cell.concentration = value;
            cells.push_back(cell);
        }
        if (cells.empty()) throw std::runtime_error("no valid grid cells found");
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

double gpi(double cold_cloud_fraction) {
    return GPI_MM_PER_HOUR * clamp_value(cold_cloud_fraction, 0.0, 1.0);
}

double gpi_with_water_vapour(double gpi_value, double precipitable_water_cm, double a1, double a2) {
    return std::exp(a1 * (precipitable_water_cm - a2)) * gpi_value;
}

MergedPrecip gpcp_merge(const std::vector<PrecipEstimate>& estimates) {
    MergedPrecip result;
    result.weights.assign(estimates.size(), 0.0);

    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    for (std::size_t i = 0; i < estimates.size(); ++i) {
        const PrecipEstimate& estimate = estimates[i];
        if (!std::isfinite(estimate.value) || !std::isfinite(estimate.error_sd) ||
            estimate.error_sd <= 0.0)
            continue;
        const double weight = 1.0 / (estimate.error_sd * estimate.error_sd);
        result.weights[i] = weight;
        weighted_sum += weight * estimate.value;
        weight_sum += weight;
        ++result.used;
    }
    if (weight_sum <= 0.0) return result;
    result.value = weighted_sum / weight_sum;
    result.valid = true;
    return result;
}

bool read_precip_estimates(const std::string& path, std::vector<PrecipEstimate>& out,
                           std::string& error) {
    try {
        const Table table = read_csv(path);
        const int source = find_column(table, {"source", "product", "name"});
        const int value = find_column(table, {"value", "precipitation", "estimate"});
        const int error_column = find_column(table, {"error_sd", "sd", "sigma", "error"});
        if (value < 0 || error_column < 0)
            throw std::runtime_error("estimate file needs value and error_sd columns");

        out.clear();
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            PrecipEstimate estimate;
            if (!parse_number(table.cell(i, value), estimate.value)) continue;
            if (!parse_number(table.cell(i, error_column), estimate.error_sd)) continue;
            if (source >= 0) estimate.source = trim_copy(table.cell(i, source));
            out.push_back(estimate);
        }
        if (out.empty()) throw std::runtime_error("no valid precipitation estimates found");
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

double toa_energy_imbalance(double effective_forcing, double feedback_parameter,
                            double temperature_change) {
    return effective_forcing + feedback_parameter * temperature_change;
}

double net_surface_radiation(double albedo, double shortwave_down, double emissivity_air,
                             double longwave_down, double emissivity_surface,
                             double surface_temperature) {
    const double emitted = STEFAN_BOLTZMANN * emissivity_surface *
                           std::pow(surface_temperature, 4.0);
    return (1.0 - albedo) * shortwave_down + emissivity_air * longwave_down - emitted;
}

double surface_energy_residual(double net_radiation, double latent_heat, double sensible_heat,
                               double ground_heat) {
    return net_radiation - (latent_heat + sensible_heat + ground_heat);
}

double carbon_growth_rate_pgc_per_year(double ppm_start, double ppm_end, double year_start,
                                       double year_end) {
    if (year_end == year_start)
        throw std::runtime_error("carbon_growth_rate_pgc_per_year: zero time interval");
    return CARBON_PPM_TO_PGC * (ppm_end - ppm_start) / (year_end - year_start);
}

double amoc_streamfunction_sv(const std::vector<VelocitySample>& samples, double depth_m) {
    if (samples.empty()) return 0.0;

    // Group by depth, integrate zonally at each level, then integrate in depth.
    std::map<double, std::vector<VelocitySample>> by_depth;
    for (const auto& sample : samples) {
        if (!std::isfinite(sample.velocity_m_s) || !std::isfinite(sample.depth_m)) continue;
        by_depth[sample.depth_m].push_back(sample);
    }
    if (by_depth.empty()) return 0.0;

    std::vector<double> depths;
    std::vector<double> zonal_transport;  // m^2/s
    for (auto& entry : by_depth) {
        std::vector<VelocitySample>& level = entry.second;
        std::sort(level.begin(), level.end(),
                  [](const VelocitySample& a, const VelocitySample& b) {
                      return a.longitude_deg < b.longitude_deg;
                  });
        double transport = 0.0;
        for (std::size_t i = 0; i + 1 < level.size(); ++i) {
            const double delta_lon_rad = deg_to_rad(level[i + 1].longitude_deg - level[i].longitude_deg);
            const double latitude_rad = deg_to_rad(0.5 * (level[i].latitude_deg + level[i + 1].latitude_deg));
            const double dx = delta_lon_rad * EARTH_RADIUS_KM * 1000.0 * std::cos(latitude_rad);
            transport += 0.5 * (level[i].velocity_m_s + level[i + 1].velocity_m_s) * dx;
        }
        depths.push_back(entry.first);
        zonal_transport.push_back(transport);
    }

    // Q(z) = integral from z to the surface (depth 0).
    double integral = 0.0;
    for (std::size_t i = 0; i + 1 < depths.size(); ++i) {
        const double top = depths[i];
        const double bottom = depths[i + 1];
        // Only the part of the layer shallower than or equal to depth_m.
        const double layer_top = std::max(top, depth_m);
        if (layer_top >= bottom) continue;
        const double thickness = bottom - layer_top;
        const double fraction = thickness / (bottom - top);
        const double top_transport = zonal_transport[i];
        const double bottom_transport =
            zonal_transport[i] + fraction * (zonal_transport[i + 1] - zonal_transport[i]);
        integral += 0.5 * (top_transport + bottom_transport) * thickness;
    }
    return integral / 1.0e6;
}

bool read_velocity_samples(const std::string& path, std::vector<VelocitySample>& out,
                           std::string& error) {
    try {
        const Table table = read_csv(path);
        const int latitude = find_column(table, {"latitude", "lat", "latitude_deg"});
        const int longitude = find_column(table, {"longitude", "lon", "longitude_deg"});
        const int depth = find_column(table, {"depth_m", "depth", "z_m"});
        const int velocity = find_column(table, {"velocity_m_s", "v", "velocity", "v_m_s"});
        if (latitude < 0 || longitude < 0 || depth < 0 || velocity < 0)
            throw std::runtime_error(
                "velocity file needs latitude, longitude, depth_m and velocity_m_s columns");

        out.clear();
        for (std::size_t i = 0; i < table.rows.size(); ++i) {
            VelocitySample sample;
            double value = 0.0;
            if (!parse_number(table.cell(i, latitude), value)) continue;
            sample.latitude_deg = value;
            if (!parse_number(table.cell(i, longitude), value)) continue;
            sample.longitude_deg = value;
            if (!parse_number(table.cell(i, depth), value)) continue;
            sample.depth_m = value;
            if (!parse_number(table.cell(i, velocity), value)) continue;
            sample.velocity_m_s = value;
            out.push_back(sample);
        }
        if (out.empty()) throw std::runtime_error("no valid velocity samples found");
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
}

}  // namespace climate
