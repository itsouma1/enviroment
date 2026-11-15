#ifndef CLIMATE_SSA_HPP
#define CLIMATE_SSA_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

struct SsaFeatures {
    double trend = 0.0;      // slope of the SSA trend component
    double amplitude = 0.0;  // seasonal amplitude
    double phase = 0.0;      // seasonal phase in radians
    double residual = 0.0;   // residual standard deviation
};

Matrix ssa_trajectory(const std::vector<double>& series, std::size_t window);

std::vector<double> ssa_trend(const std::vector<double>& series, std::size_t window);

std::vector<double> sinusoid_fit(const std::vector<double>& series, double period,
                                 double& amplitude, double& phase);

SsaFeatures ssa_features(const std::vector<double>& series, std::size_t window,
                         double period);

}  // namespace climate

#endif
