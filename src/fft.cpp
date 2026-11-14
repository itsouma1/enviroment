#include "fft.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "common.hpp"

namespace climate {
namespace {

void fft_radix2(std::vector<Complex>& data, bool inverse) {
    const std::size_t n = data.size();
    if (n <= 1) return;

    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }

    for (std::size_t length = 2; length <= n; length <<= 1) {
        const double angle = 2.0 * PI / static_cast<double>(length) * (inverse ? 1.0 : -1.0);
        const Complex root(std::cos(angle), std::sin(angle));
        for (std::size_t start = 0; start < n; start += length) {
            Complex factor(1.0, 0.0);
            const std::size_t half = length >> 1;
            for (std::size_t k = 0; k < half; ++k) {
                const Complex even = data[start + k];
                const Complex odd = data[start + k + half] * factor;
                data[start + k] = even + odd;
                data[start + k + half] = even - odd;
                factor *= root;
            }
        }
    }

    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (auto& value : data) value *= scale;
    }
}

}  // namespace

void fft(std::vector<Complex>& data, bool inverse) {
    const std::size_t n = data.size();
    if (n <= 1) return;
    if ((n & (n - 1)) == 0) {
        fft_radix2(data, inverse);
        return;
    }

    // Bluestein: X_k = w_k * sum_j (x_j w_j) conj(w_{k-j}), where
    // w_j = exp(s * i * pi * j^2 / n) and s = -1 for the forward transform.
    std::size_t m = 1;
    while (m < 2 * n - 1) m <<= 1;

    const double sign = inverse ? 1.0 : -1.0;
    std::vector<Complex> a(m, Complex(0.0, 0.0));
    std::vector<Complex> b(m, Complex(0.0, 0.0));
    std::vector<Complex> chirp(n);

    for (std::size_t j = 0; j < n; ++j) {
        const double angle = sign * PI * (static_cast<double>(j) * static_cast<double>(j)) /
                             static_cast<double>(n);
        chirp[j] = Complex(std::cos(angle), std::sin(angle));
        a[j] = data[j] * chirp[j];
        b[j] = std::conj(chirp[j]);
        if (j != 0) b[m - j] = std::conj(chirp[j]);
    }

    fft_radix2(a, false);
    fft_radix2(b, false);
    for (std::size_t i = 0; i < m; ++i) a[i] *= b[i];
    fft_radix2(a, true);

    for (std::size_t k = 0; k < n; ++k) data[k] = a[k] * chirp[k];
    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (auto& value : data) value *= scale;
    }
}

namespace {

// Reflects a series about its ends until it reaches the target length.
std::vector<double> reflect_pad(const std::vector<double>& values, std::size_t target) {
    const std::size_t n = values.size();
    std::vector<double> padded;
    padded.reserve(target);
    const std::size_t period = 2 * n;
    for (std::size_t i = 0; i < target; ++i) {
        const std::size_t folded = i % period;
        const std::size_t index = folded < n ? folded : period - 1 - folded;
        padded.push_back(values[index]);
    }
    return padded;
}

}  // namespace

std::vector<double> lowpass_filter(const std::vector<double>& values, double cutoff,
                                   double sample_spacing) {
    if (values.empty()) return values;
    if (cutoff <= 0.0 || sample_spacing <= 0.0)
        throw std::runtime_error("lowpass_filter: cutoff and spacing must be positive");

    const std::size_t n = values.size();
    if (n < 4) return values;

    std::size_t length = 1;
    while (length < 2 * n) length <<= 1;

    const double offset = values.front();
    std::vector<Complex> spectrum(length, Complex(0.0, 0.0));
    const std::vector<double> padded = reflect_pad(values, length);
    for (std::size_t i = 0; i < length; ++i) spectrum[i] = Complex(padded[i] - offset, 0.0);

    fft(spectrum, false);

    const double cutoff_frequency = 1.0 / cutoff;
    for (std::size_t k = 1; k <= length / 2; ++k) {
        const double frequency = static_cast<double>(k) / (static_cast<double>(length) * sample_spacing);
        if (frequency > cutoff_frequency) {
            spectrum[k] = Complex(0.0, 0.0);
            if (k != length - k) spectrum[length - k] = Complex(0.0, 0.0);
        }
    }

    fft(spectrum, true);

    std::vector<double> result(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) result[i] = spectrum[i].real() + offset;
    return result;
}

std::vector<double> highpass_filter(const std::vector<double>& values, double cutoff,
                                    double sample_spacing) {
    const std::vector<double> smooth = lowpass_filter(values, cutoff, sample_spacing);
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) result[i] = values[i] - smooth[i];
    return result;
}

}  // namespace climate
