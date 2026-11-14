#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace climate {

// Polynomial-plus-harmonics model used for the NOAA/CMDL CO2 curve fitting
// method (Thoning, Tans & Komhyr, 1989):
//   f(t) = a0 + a1 t + ... + ak t^k
//          + sum_n [A_n cos(2 pi n t / period) + B_n sin(2 pi n t / period)]
struct HarmonicModel {
    int polynomial_order = 2;
    int harmonic_count = 4;
    double period = 1.0;       // in the units of t (1.0 = one year for fractional years)
    double time_mean = 0.0;    // centring used internally for conditioning
    double time_scale = 1.0;
    std::vector<double> coefficients;
    std::size_t used_points = 0;
    std::size_t rejected_points = 0;
    double residual_sd = 0.0;
    double r_squared = 0.0;
    bool converged = true;
};

struct HarmonicComponent {
    int number = 1;
    double amplitude = 0.0;
    double phase_radians = 0.0;
    double phase_degrees = 0.0;
};

// Fits the model with iteratively reweighted outlier rejection: points whose
// residual exceeds reject_threshold residual standard deviations are removed
// and the model is refitted. Set reject_iterations to 0 to disable.
HarmonicModel fit_harmonic(const std::vector<double>& time, const std::vector<double>& values,
                           int polynomial_order, int harmonic_count, double period,
                           int reject_iterations = 4, double reject_threshold = 3.0);

double harmonic_evaluate(const HarmonicModel& model, double t);
double harmonic_derivative(const HarmonicModel& model, double t);

// Derivative of the polynomial (trend) part only. The CO2 growth rate is defined
// on the long-term trend, not on the seasonal harmonics.
double harmonic_trend_derivative(const HarmonicModel& model, double t);
double harmonic_seasonal_value(const HarmonicModel& model, double t);
double harmonic_seasonal_amplitude(const HarmonicModel& model);
std::vector<HarmonicComponent> harmonic_components(const HarmonicModel& model);

// Mean of the model derivative over [start, end], computed with adaptive-like
// fine trapezoidal sampling; used for CO2 growth rates.
double harmonic_mean_derivative(const HarmonicModel& model, double start, double end, int samples = 512);

// Mean of the polynomial-part derivative over [start, end].
double harmonic_mean_trend_derivative(const HarmonicModel& model, double start, double end,
                                      int samples = 512);

}  // namespace climate
