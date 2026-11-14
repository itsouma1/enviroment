#pragma once

#include <cstddef>
#include <vector>

namespace climate {

// Fills missing (non-finite) entries by linear interpolation between the nearest
// valid neighbours. Leading and trailing gaps are filled with the nearest valid
// value. filled receives the number of substituted entries.
std::vector<double> interpolate_missing(const std::vector<double>& values,
                                        std::size_t* filled = nullptr);

struct Climatology {
    std::vector<double> mean;
    std::vector<double> sd;
    std::vector<std::size_t> count;
};

// Calendar-phase climatology: phase i is the set of samples with index % phase_count == i.
Climatology phase_climatology(const std::vector<double>& values, std::size_t phase_count);

// (value - mean_phase) / sd_phase; zero deviation when the phase has no spread.
std::vector<double> standardized_anomalies(const std::vector<double>& values,
                                           const Climatology& climatology);

// (value - mean_phase) without scaling.
std::vector<double> anomalies(const std::vector<double>& values, const Climatology& climatology);

// Deseasonalises by subtracting the phase climatology.
std::vector<double> remove_seasonal_cycle(const std::vector<double>& values,
                                          std::size_t phase_count);

// Mean spacing of the time vector (used to report rates per decade).
double mean_spacing(const std::vector<double>& time);

}  // namespace climate
