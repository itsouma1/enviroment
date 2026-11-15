#include "precip.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "ensemble.hpp"
#include "stats.hpp"

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double wet_day_probability(double p01, double p11) {
    require(p01 >= 0.0 && p01 <= 1.0 && p11 >= 0.0 && p11 <= 1.0,
            "wet_day_probability: probabilities");
    const double denominator = 1.0 - p11 + p01;
    require(denominator > 0.0, "wet_day_probability: degenerate chain");
    return p01 / denominator;
}

Matrix markov_occurrence_field(std::size_t rows, std::size_t columns, std::size_t days,
                               double p01, double p11, double neighbour_correlation,
                               unsigned seed) {
    const std::size_t sites = rows * columns;
    return markov_occurrence_field_sites(rows, columns, days,
                                         std::vector<double>(sites, p01),
                                         std::vector<double>(sites, p11),
                                         neighbour_correlation, seed);
}

Matrix markov_occurrence_field_sites(std::size_t rows, std::size_t columns,
                                     std::size_t days, const std::vector<double>& p01,
                                     const std::vector<double>& p11,
                                     double neighbour_correlation, unsigned seed) {
    const std::size_t sites = rows * columns;
    require(sites > 0 && days > 0 && p01.size() == sites && p11.size() == sites,
            "markov_occurrence_field_sites: size");
    for (std::size_t i = 0; i < sites; ++i)
        require(p01[i] >= 0.0 && p01[i] <= 1.0 && p11[i] >= 0.0 && p11[i] <= 1.0,
                "markov_occurrence_field_sites: probability");
    const Matrix correlation = spatial_correlation_matrix(rows, columns, neighbour_correlation);
    Matrix factor;
    require(cholesky_factor(correlation, factor),
            "markov_occurrence_field_sites: invalid correlation");
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::vector<double> uniform(sites, 0.0), previous(sites, 0.0);
    Matrix occurrence(days, std::vector<double>(sites, 0.0));
    for (std::size_t day = 0; day < days; ++day) {
        // Build the spatial field from one innovation vector.
        std::vector<double> innovation(sites, 0.0);
        for (double& value : innovation) value = normal(generator);
        for (std::size_t site = 0; site < sites; ++site) {
            double sum = 0.0;
            for (std::size_t k = 0; k <= site; ++k) sum += factor[site][k] * innovation[k];
            uniform[site] = normal_cdf(sum);
        }
        for (std::size_t site = 0; site < sites; ++site) {
            const double probability =
                day == 0 ? wet_day_probability(p01[site], p11[site])
                         : (previous[site] > 0.5 ? p11[site] : p01[site]);
            occurrence[day][site] = uniform[site] < probability ? 1.0 : 0.0;
            previous[site] = occurrence[day][site];
        }
    }
    return occurrence;
}

PrecipitationField simulate_precipitation(std::size_t rows, std::size_t columns,
                                         std::size_t days, double p01, double p11,
                                         double neighbour_correlation, double shape,
                                         double rate, unsigned seed) {
    require(shape > 0.0 && rate > 0.0, "simulate_precipitation: gamma parameters");
    PrecipitationField result;
    result.occurrence = markov_occurrence_field(rows, columns, days, p01, p11,
                                                neighbour_correlation, seed);
    const std::size_t sites = rows * columns;
    result.amount.assign(days, std::vector<double>(sites, 0.0));
    std::mt19937 generator(seed + 104729u);
    std::gamma_distribution<double> gamma(shape, 1.0 / rate);
    for (std::size_t day = 0; day < days; ++day)
        for (std::size_t site = 0; site < sites; ++site)
            if (result.occurrence[day][site] > 0.5) result.amount[day][site] = gamma(generator);
    return result;
}

double wilks_correlation_distance_derivative(double distance, double d1, double d2) {
    require(distance >= 0.0, "wilks_correlation_distance_derivative: distance");
    const double correlation = wilks_correlation(distance, d1, d2);
    return correlation * std::pow(distance, d2);
}

double wilks_correlation_shape_derivative(double distance, double d1, double d2) {
    require(distance > 0.0, "wilks_correlation_shape_derivative: distance");
    const double correlation = wilks_correlation(distance, d1, d2);
    return correlation * d1 * std::pow(distance, d2) * std::log(distance);
}

double temperature_correlation(double distance, double elevation_difference, double b0,
                               double b1, double b2) {
    require(distance >= 0.0, "temperature_correlation: distance");
    return std::exp(b0 + b1 * distance + b2 * std::fabs(elevation_difference));
}

double locally_weighted_parameter(const std::vector<double>& site_x,
                                  const std::vector<double>& site_y,
                                  const std::vector<double>& values, double grid_x,
                                  double grid_y, double bandwidth) {
    require(site_x.size() == site_y.size() && site_x.size() == values.size() &&
                !site_x.empty() && bandwidth > 0.0,
            "locally_weighted_parameter: inputs");
    double weighted = 0.0, total = 0.0;
    for (std::size_t s = 0; s < site_x.size(); ++s) {
        const double dx = site_x[s] - grid_x, dy = site_y[s] - grid_y;
        const double weight = std::exp(-(dx * dx + dy * dy) / (2.0 * bandwidth * bandwidth));
        weighted += weight * values[s];
        total += weight;
    }
    require(total > 0.0, "locally_weighted_parameter: zero weight");
    return weighted / total;
}

double locally_weighted_parameter_elevation(
    const std::vector<double>& site_x, const std::vector<double>& site_y,
    const std::vector<double>& site_elevation, const std::vector<double>& values,
    double grid_x, double grid_y, double grid_elevation, double bandwidth,
    double elevation_bandwidth) {
    require(site_x.size() == site_y.size() && site_x.size() == values.size() &&
                site_x.size() == site_elevation.size() && !site_x.empty() &&
                bandwidth > 0.0 && elevation_bandwidth > 0.0,
            "locally_weighted_parameter_elevation: inputs");
    double weighted = 0.0, total = 0.0;
    for (std::size_t s = 0; s < site_x.size(); ++s) {
        const double dx = site_x[s] - grid_x, dy = site_y[s] - grid_y;
        const double dh = site_elevation[s] - grid_elevation;
        const double weight =
            std::exp(-(dx * dx + dy * dy) / (2.0 * bandwidth * bandwidth) -
                     std::fabs(dh) / (2.0 * elevation_bandwidth * elevation_bandwidth));
        weighted += weight * values[s];
        total += weight;
    }
    require(total > 0.0, "locally_weighted_parameter_elevation: zero weight");
    return weighted / total;
}

}  // namespace climate
