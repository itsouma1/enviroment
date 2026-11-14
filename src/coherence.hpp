#ifndef CLIMATE_COHERENCE_HPP
#define CLIMATE_COHERENCE_HPP

#include <complex>
#include <vector>

namespace climate {

using Complex = std::complex<double>;
using ComplexField = std::vector<std::vector<Complex>>;  // [scale][time]

// Continuous Morlet wavelet transform, psi(eta) = pi^-1/4 exp(i omega0 eta)
// exp(-eta^2 / 2), sampled at the given scales.
ComplexField morlet_cwt(const std::vector<double>& values, double time_step,
                        const std::vector<double>& scales, double omega0 = 6.0);

// Wavelet coherence R^2 at one time-scale point, smoothed over a box.
double wavelet_coherence(const ComplexField& first, const ComplexField& second,
                         std::size_t time_index, std::size_t scale_index, int time_window,
                         int scale_window);

}  // namespace climate

#endif
