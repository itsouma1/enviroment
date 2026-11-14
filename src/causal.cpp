#include "causal.hpp"

#include <map>
#include <stdexcept>

namespace climate {

Matrix sem_interaction_design(const std::vector<double>& temperature,
                              const std::vector<double>& precipitation) {
    if (temperature.size() != precipitation.size())
        throw std::runtime_error("sem_interaction_design: size mismatch");
    Matrix design;
    design.reserve(temperature.size());
    for (std::size_t i = 0; i < temperature.size(); ++i)
        design.push_back({temperature[i], precipitation[i], temperature[i] * precipitation[i]});
    return design;
}

std::vector<double> two_way_demean(const std::vector<double>& values,
                                   const std::vector<int>& units,
                                   const std::vector<int>& times) {
    if (values.size() != units.size() || values.size() != times.size())
        throw std::runtime_error("two_way_demean: size mismatch");
    if (values.empty()) return {};

    std::map<int, double> unit_sum, time_sum;
    std::map<int, std::size_t> unit_count, time_count;
    double grand_sum = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        unit_sum[units[i]] += values[i];
        time_sum[times[i]] += values[i];
        ++unit_count[units[i]];
        ++time_count[times[i]];
        grand_sum += values[i];
    }
    const double grand_mean = grand_sum / static_cast<double>(values.size());

    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double unit_mean = unit_sum[units[i]] / static_cast<double>(unit_count[units[i]]);
        const double time_mean = time_sum[times[i]] / static_cast<double>(time_count[times[i]]);
        result[i] = values[i] - unit_mean - time_mean + grand_mean;
    }
    return result;
}

}  // namespace climate
