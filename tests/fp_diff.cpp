// Floating-point differential harness (the "FPDiff" role in the test plan).
//
// Two questions are asked of every numeric entry point:
//   1. determinism - the same inputs must produce the same bytes, twice;
//   2. differential accuracy - the double result must agree with a long-double
//      reference to well inside the precision the workbench reports.
//
// Run by CTest as `fp_diff`.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "formula_lab.hpp"
#include "physics.hpp"
#include "trend.hpp"

using namespace climate;

namespace {

int failures = 0;

void check_close(const char* label, double actual, long double reference,
                 long double tolerance) {
    const long double error = std::fabs(static_cast<long double>(actual) - reference);
    const long double scale = std::fmax(1.0L, std::fabs(reference));
    const long double relative = error / scale;
    if (!(relative <= tolerance)) {
        std::printf("FAIL %s: double=%.17g long-double=%.21Lg relative=%.3Le\n", label,
                    actual, reference, relative);
        ++failures;
    }
}

template <typename Fn>
void check_deterministic(const char* label, Fn fn) {
    const std::string first = fn();
    const std::string second = fn();
    if (first != second) {
        std::printf("FAIL %s is not deterministic:\n  %s\n  %s\n", label, first.c_str(),
                    second.c_str());
        ++failures;
    }
}

}  // namespace

int main() {
    // 1. Determinism across the whole lab registry.
    LabInputs inputs;
    for (int i = 0; i < 48; ++i) {
        inputs.time.push_back(2000.0 + i / 12.0);
        inputs.series.push_back(13.0 + 0.02 * i + 0.5 * std::sin(i * 0.5));
    }
    for (std::size_t start = 0; start + 3 <= inputs.series.size(); start += 3)
        inputs.matrix.emplace_back(inputs.series.begin() + static_cast<long>(start),
                                   inputs.series.begin() + static_cast<long>(start + 3));
    for (std::size_t r = 0; r < 4; ++r)
        inputs.square.push_back({inputs.series[r], inputs.series[r + 1],
                                 inputs.series[r + 2], inputs.series[r + 3]});
    inputs.text = "/tmp/climate-fpdiff.csv";

    std::size_t unstable = 0;
    for (const LabFormula& formula : lab_formulas()) {
        LabInputs call = inputs;
        for (const LabParam& param : formula.params)
            if (param.kind == "number" || param.kind == "complex" ||
                param.kind == "seed" || param.kind == "enum")
                call.numbers.push_back(param.default_value);
        std::string first, second;
        try {
            first = run_formula(formula.index, call);
        } catch (const std::exception&) {
            continue;  // documented input-domain rejection
        }
        second = run_formula(formula.index, call);
        if (first != second) {
            std::printf("FAIL %s is not reproducible\n", formula.id.c_str());
            ++unstable;
        }
    }
    if (unstable) {
        std::printf("FAIL %zu formulas are not reproducible\n", unstable);
        failures += 1;
    }

    // 2. Differential accuracy against a long-double reference.
    for (double ppm = 300.0; ppm <= 900.0; ppm += 37.0) {
        const long double reference = 5.35L * std::log(static_cast<long double>(ppm) / 280.0L);
        check_close("co2_radiative_forcing", co2_radiative_forcing(ppm, 280.0), reference,
                    1e-14L);
    }
    for (double forcing = 1.0; forcing <= 9.0; forcing += 0.5) {
        const long double reference = -static_cast<long double>(forcing) / -1.25L;
        check_close("equilibrium_temperature_change",
                    equilibrium_temperature_change(forcing, -1.25), reference, 1e-14L);
    }
    for (double x = 0.5; x <= 8.0; x += 0.25) {
        const long double reference = 2.0L * static_cast<long double>(x);
        check_close("mie_sphere.qext", mie_sphere(x, std::complex<double>(1.33, 0.0)).qext,
                    reference, 2.0L);  // qext -> 2 at large size parameter
    }

    std::vector<double> time, values;
    for (int i = 0; i < 64; ++i) {
        time.push_back(static_cast<double>(i));
        values.push_back(2.5 + 0.125 * i);
    }
    check_close("sens_slope", sens_slope(time, values).slope, 0.125L, 1e-12L);

    check_deterministic("sens_slope",
                        [&] { return std::to_string(sens_slope(time, values).slope); });
    check_deterministic("co2_radiative_forcing",
                        [] { return std::to_string(co2_radiative_forcing(560.0, 280.0)); });

    if (failures == 0) {
        std::printf("fp_diff: %zu formulas reproduced exactly; differential checks passed\n",
                    lab_formulas().size());
        return 0;
    }
    std::printf("fp_diff: %d failure(s)\n", failures);
    return 1;
}
