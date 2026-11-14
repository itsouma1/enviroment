#pragma once

#include <complex>
#include <vector>

namespace climate {

using Complex = std::complex<double>;

// In-place discrete Fourier transform of arbitrary length. Power-of-two sizes
// use the radix-2 Cooley-Tukey algorithm; other sizes use Bluestein's chirp-z
// algorithm so that no zero padding or approximation is required.
void fft(std::vector<Complex>& data, bool inverse);

// Thoning et al. (1989) style low-pass filter: the series is reflected about its
// ends to reduce edge effects, transformed, all harmonics with a period shorter
// than cutoff (in the same units as the sample spacing) are removed, and the
// series is transformed back.
std::vector<double> lowpass_filter(const std::vector<double>& values, double cutoff,
                                   double sample_spacing = 1.0);

// High-pass complement of lowpass_filter.
std::vector<double> highpass_filter(const std::vector<double>& values, double cutoff,
                                    double sample_spacing = 1.0);

}  // namespace climate
