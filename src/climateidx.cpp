#include "climateidx.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

double area_mean_anomaly(const std::vector<double>& anomalies,
                         const std::vector<double>& weights) {
    require(!anomalies.empty() && anomalies.size() == weights.size(),
            "area_mean_anomaly: size mismatch");
    double weighted = 0.0, total = 0.0;
    for (std::size_t i = 0; i < anomalies.size(); ++i) {
        require(weights[i] >= 0.0, "area_mean_anomaly: negative weight");
        weighted += anomalies[i] * weights[i];
        total += weights[i];
    }
    require(total > 0.0, "area_mean_anomaly: zero weight");
    return weighted / total;
}

double relative_nino34(double nino34, double nino34_climatology, double tropical_mean,
                       double tropical_climatology, double scale) {
    return scale * ((nino34 - nino34_climatology) - (tropical_mean - tropical_climatology));
}

double amo_relative(double north_atlantic_anomaly, double global_anomaly) {
    return north_atlantic_anomaly - global_anomaly;
}

double pdo_projection(const std::vector<double>& leading_eof,
                      const std::vector<double>& anomaly_field) {
    require(leading_eof.size() == anomaly_field.size(), "pdo_projection: size mismatch");
    double projection = 0.0;
    for (std::size_t i = 0; i < leading_eof.size(); ++i)
        projection += leading_eof[i] * anomaly_field[i];
    return projection;
}

double amoc_index(const std::vector<std::vector<double>>& streamfunction) {
    require(!streamfunction.empty(), "amoc_index: empty field");
    bool found = false;
    double maximum = 0.0;
    for (const auto& row : streamfunction) {
        for (double value : row) {
            if (!found || value > maximum) {
                maximum = value;
                found = true;
            }
        }
    }
    return maximum;
}

double ocean_heat_content(double density, double specific_heat,
                          const std::vector<double>& temperatures,
                          const std::vector<double>& climatology,
                          const std::vector<double>& layer_thickness) {
    require(temperatures.size() == climatology.size() &&
                temperatures.size() == layer_thickness.size(),
            "ocean_heat_content: size mismatch");
    double content = 0.0;
    for (std::size_t i = 0; i < temperatures.size(); ++i)
        content += density * specific_heat * (temperatures[i] - climatology[i]) *
                   layer_thickness[i];
    return content;
}

double layer_heat_content_anomaly(double density, double specific_heat,
                                  double temperature_anomaly, double layer_thickness) {
    return density * specific_heat * temperature_anomaly * layer_thickness;
}

double earth_energy_imbalance(double incoming, double outgoing_longwave,
                              double outgoing_shortwave) {
    return incoming - (outgoing_longwave + outgoing_shortwave);
}

double energy_imbalance_decomposition(double forcing, double feedback_parameter,
                                      double temperature_anomaly) {
    return forcing - feedback_parameter * temperature_anomaly;
}

double sea_level_budget(double sterodynamic, double gravity_rotation_deformation,
                        double glacial_isostatic_adjustment, double residual) {
    return sterodynamic + gravity_rotation_deformation + glacial_isostatic_adjustment +
           residual;
}

double barystatic_sea_level(double greenland, double antarctic, double glaciers,
                            double terrestrial_water, double atmospheric_water) {
    return greenland + antarctic + glaciers + terrestrial_water + atmospheric_water;
}

double quadratic_sea_level(double offset, double rate, double acceleration, double time) {
    return offset + rate * time + acceleration * time * time;
}

double quadratic_sea_level_acceleration(double acceleration) {
    return 2.0 * acceleration;
}

double broken_linear_sea_level(double offset, double rate, double slope_change,
                               double time, double change_point) {
    const double stepped = time >= change_point ? time - change_point : 0.0;
    return offset + rate * time + slope_change * stepped;
}

void DynamicLinearState::update(double observation, double observation_variance,
                                double level_noise, double rate_noise) {
    require(observation_variance > 0.0, "DynamicLinearState: observation variance");

    // Transition F = [[1, 1], [0, 1]]: level += rate.
    const double predicted_level = level + rate;
    const double predicted_rate = rate;
    const double predicted_00 = covariance_00 + 2.0 * covariance_01 + covariance_11 + level_noise;
    const double predicted_01 = covariance_01 + covariance_11;
    const double predicted_11 = covariance_11 + rate_noise;

    const double innovation = observation - predicted_level;
    const double innovation_variance = predicted_00 + observation_variance;
    const double gain_level = predicted_00 / innovation_variance;
    const double gain_rate = predicted_01 / innovation_variance;

    level = predicted_level + gain_level * innovation;
    rate = predicted_rate + gain_rate * innovation;

    covariance_00 = (1.0 - gain_level) * predicted_00;
    covariance_01 = (1.0 - gain_level) * predicted_01;
    covariance_11 = predicted_11 - gain_rate * predicted_01;
}

double groundwater_storage_anomaly(double total_storage, double snow_water_equivalent,
                                   double surface_water, double soil_moisture) {
    return total_storage - snow_water_equivalent - surface_water - soil_moisture;
}

double recharge_convolution(double initial_head, const std::vector<double>& recharge,
                            double characteristic_time) {
    require(characteristic_time > 0.0, "recharge_convolution: time constant");
    double head = initial_head;
    const std::size_t last = recharge.size() - 1;
    for (std::size_t tau = 0; tau <= last; ++tau) {
        const std::size_t index = last - tau;
        head += recharge[index] * std::exp(-static_cast<double>(tau) / characteristic_time);
    }
    return head;
}

double storage_depletion(const std::vector<double>& drawdown,
                         const std::vector<double>& cell_area,
                         double storage_coefficient) {
    require(drawdown.size() == cell_area.size(), "storage_depletion: size mismatch");
    double volume = 0.0;
    for (std::size_t i = 0; i < drawdown.size(); ++i)
        volume += drawdown[i] * cell_area[i] * storage_coefficient;
    return volume;
}

Matrix hovmoller_average(const std::vector<Matrix>& fields) {
    require(!fields.empty() && !fields[0].empty(), "hovmoller_average: empty field");
    const std::size_t columns = fields[0][0].size();
    require(columns > 0, "hovmoller_average: no columns");
    Matrix result(fields.size(), std::vector<double>(columns, 0.0));
    for (std::size_t time = 0; time < fields.size(); ++time) {
        require(!fields[time].empty(), "hovmoller_average: empty time slice");
        for (const auto& row : fields[time])
            require(row.size() == columns, "hovmoller_average: ragged field");
        for (std::size_t column = 0; column < columns; ++column) {
            double sum = 0.0;
            for (const auto& row : fields[time]) sum += row[column];
            result[time][column] = sum / static_cast<double>(fields[time].size());
        }
    }
    return result;
}

std::vector<double> phase_space_derivative(const std::vector<double>& values,
                                          double time_step) {
    require(values.size() >= 2, "phase_space_derivative: need at least two values");
    require(time_step > 0.0, "phase_space_derivative: time step must be positive");
    std::vector<double> derivative(values.size(), 0.0);
    for (std::size_t i = 1; i + 1 < values.size(); ++i)
        derivative[i] = (values[i + 1] - values[i - 1]) / (2.0 * time_step);
    derivative.front() = (values[1] - values[0]) / time_step;
    derivative.back() = (values.back() - values[values.size() - 2]) / time_step;
    return derivative;
}

double heston_variance_step(double variance, double mean_reversion, double long_run_variance,
                            double volatility_of_volatility, double time_step, double shock) {
    require(variance >= 0.0, "heston_variance_step: negative variance");
    return variance + mean_reversion * (long_run_variance - variance) * time_step +
           volatility_of_volatility * std::sqrt(variance) * shock;
}

std::vector<double> markov_stationary_distribution(
    const std::vector<std::vector<double>>& transition) {
    const std::size_t n = transition.size();
    require(n > 0, "markov_stationary_distribution: empty matrix");
    std::vector<double> distribution(n, 1.0 / static_cast<double>(n));
    for (int iteration = 0; iteration < 10000; ++iteration) {
        std::vector<double> next(n, 0.0);
        for (std::size_t from = 0; from < n; ++from) {
            require(transition[from].size() == n, "markov_stationary_distribution: ragged");
            for (std::size_t to = 0; to < n; ++to)
                next[to] += distribution[from] * transition[from][to];
        }
        double change = 0.0;
        for (std::size_t i = 0; i < n; ++i) change += std::fabs(next[i] - distribution[i]);
        distribution.swap(next);
        if (change < 1e-14) break;
    }
    return distribution;
}

double regime_expected_value(const std::vector<double>& regime_means,
                             const std::vector<double>& regime_probabilities) {
    require(regime_means.size() == regime_probabilities.size(),
            "regime_expected_value: size mismatch");
    double expectation = 0.0;
    for (std::size_t i = 0; i < regime_means.size(); ++i)
        expectation += regime_means[i] * regime_probabilities[i];
    return expectation;
}

}  // namespace climate
