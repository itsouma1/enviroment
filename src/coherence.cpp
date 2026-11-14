#include "coherence.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

ComplexField morlet_cwt(const std::vector<double>& values, double time_step,
                        const std::vector<double>& scales, double omega0) {
    require(!values.empty(), "morlet_cwt: empty series");
    require(time_step > 0.0, "morlet_cwt: time step must be positive");
    require(!scales.empty(), "morlet_cwt: no scales");
    for (double scale : scales) require(scale > 0.0, "morlet_cwt: scale must be positive");

    const double normalisation = std::pow(kPi, -0.25);
    ComplexField field(scales.size(), std::vector<Complex>(values.size(), Complex(0.0, 0.0)));
    for (std::size_t s = 0; s < scales.size(); ++s) {
        for (std::size_t n = 0; n < values.size(); ++n) {
            Complex sum(0.0, 0.0);
            for (std::size_t m = 0; m < values.size(); ++m) {
                const double eta = (static_cast<double>(m) - static_cast<double>(n)) *
                                   time_step / scales[s];
                const double envelope = normalisation * std::exp(-0.5 * eta * eta);
                const Complex wavelet = envelope * std::exp(Complex(0.0, -omega0 * eta));
                sum += values[m] * wavelet;
            }
            field[s][n] = sum;
        }
    }
    return field;
}

double wavelet_coherence(const ComplexField& first, const ComplexField& second,
                         std::size_t time_index, std::size_t scale_index, int time_window,
                         int scale_window) {
    require(first.size() == second.size(), "wavelet_coherence: scale mismatch");
    require(scale_index < first.size() && time_index < first[scale_index].size(),
            "wavelet_coherence: index out of range");
    require(time_window >= 0 && scale_window >= 0, "wavelet_coherence: negative window");

    const int time_begin = std::max(0, static_cast<int>(time_index) - time_window);
    const int time_end =
        std::min(static_cast<int>(first[scale_index].size()) - 1,
                 static_cast<int>(time_index) + time_window);
    const int scale_begin = std::max(0, static_cast<int>(scale_index) - scale_window);
    const int scale_end = std::min(static_cast<int>(first.size()) - 1,
                                   static_cast<int>(scale_index) + scale_window);

    Complex cross(0.0, 0.0);
    double power_first = 0.0, power_second = 0.0;
    int count = 0;
    for (int s = scale_begin; s <= scale_end; ++s) {
        require(first[static_cast<std::size_t>(s)].size() ==
                    first[scale_index].size() &&
                    second[static_cast<std::size_t>(s)].size() == first[scale_index].size(),
                "wavelet_coherence: time size mismatch");
        for (int t = time_begin; t <= time_end; ++t) {
            const Complex& a = first[static_cast<std::size_t>(s)][static_cast<std::size_t>(t)];
            const Complex& b = second[static_cast<std::size_t>(s)][static_cast<std::size_t>(t)];
            cross += a * std::conj(b);
            power_first += std::norm(a);
            power_second += std::norm(b);
            ++count;
        }
    }
    require(count > 0, "wavelet_coherence: empty smoothing window");
    cross /= static_cast<double>(count);
    power_first /= static_cast<double>(count);
    power_second /= static_cast<double>(count);
    const double denominator = power_first * power_second;
    if (denominator <= 0.0) return 0.0;
    return std::norm(cross) / denominator;
}

}  // namespace climate
