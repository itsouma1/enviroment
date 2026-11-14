#include "series.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "stats.hpp"

namespace climate {

std::vector<double> interpolate_missing(const std::vector<double>& values, std::size_t* filled) {
    const std::size_t n = values.size();
    std::vector<double> result(values);
    std::size_t substitutions = 0;

    std::size_t first_valid = n;
    for (std::size_t i = 0; i < n; ++i) {
        if (std::isfinite(values[i])) {
            first_valid = i;
            break;
        }
    }
    if (first_valid == n) {
        if (filled != nullptr) *filled = n;
        return result;
    }

    for (std::size_t i = 0; i < first_valid; ++i) {
        result[i] = values[first_valid];
        ++substitutions;
    }

    std::size_t previous = first_valid;
    for (std::size_t i = first_valid + 1; i < n; ++i) {
        if (!std::isfinite(values[i])) continue;
        const std::size_t gap = i - previous - 1;
        if (gap > 0) {
            const double start = values[previous];
            const double end = values[i];
            for (std::size_t k = 1; k <= gap; ++k) {
                const double fraction = static_cast<double>(k) / static_cast<double>(gap + 1);
                result[previous + k] = start + (end - start) * fraction;
                ++substitutions;
            }
        }
        previous = i;
    }

    for (std::size_t i = previous + 1; i < n; ++i) {
        result[i] = values[previous];
        ++substitutions;
    }

    if (filled != nullptr) *filled = substitutions;
    return result;
}

Climatology phase_climatology(const std::vector<double>& values, std::size_t phase_count) {
    if (phase_count == 0) throw std::runtime_error("phase_climatology: phase count must be > 0");
    Climatology result;
    result.mean.assign(phase_count, 0.0);
    result.sd.assign(phase_count, 0.0);
    result.count.assign(phase_count, 0);

    std::vector<std::vector<double>> buckets(phase_count);
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i])) continue;
        buckets[i % phase_count].push_back(values[i]);
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        result.count[phase] = buckets[phase].size();
        if (buckets[phase].empty()) continue;
        result.mean[phase] = mean_value(buckets[phase]);
        result.sd[phase] = stdev_value(buckets[phase]);
    }
    return result;
}

std::vector<double> standardized_anomalies(const std::vector<double>& values,
                                           const Climatology& climatology) {
    std::vector<double> result(values.size(), 0.0);
    if (climatology.mean.empty()) return result;
    const std::size_t phases = climatology.mean.size();
    for (std::size_t i = 0; i < values.size(); ++i) {
        const std::size_t phase = i % phases;
        const double deviation = values[i] - climatology.mean[phase];
        result[i] = climatology.sd[phase] > 0.0 ? deviation / climatology.sd[phase] : 0.0;
    }
    return result;
}

std::vector<double> anomalies(const std::vector<double>& values, const Climatology& climatology) {
    std::vector<double> result(values.size(), 0.0);
    if (climatology.mean.empty()) return result;
    const std::size_t phases = climatology.mean.size();
    for (std::size_t i = 0; i < values.size(); ++i)
        result[i] = values[i] - climatology.mean[i % phases];
    return result;
}

std::vector<double> remove_seasonal_cycle(const std::vector<double>& values,
                                          std::size_t phase_count) {
    return anomalies(values, phase_climatology(values, phase_count));
}

double mean_spacing(const std::vector<double>& time) {
    if (time.size() < 2) return 1.0;
    const double span = time.back() - time.front();
    if (span <= 0.0) return 1.0;
    return span / static_cast<double>(time.size() - 1);
}

}  // namespace climate
