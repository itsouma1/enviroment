#ifndef CLIMATE_CAUSAL_HPP
#define CLIMATE_CAUSAL_HPP

#include <vector>

#include "ml.hpp"

namespace climate {

// Design columns [Temp, Precip, Temp * Precip] for the panel SEM.
Matrix sem_interaction_design(const std::vector<double>& temperature,
                              const std::vector<double>& precipitation);

// y - unit mean - time mean + grand mean.
std::vector<double> two_way_demean(const std::vector<double>& values,
                                   const std::vector<int>& units,
                                   const std::vector<int>& times);

}  // namespace climate

#endif
