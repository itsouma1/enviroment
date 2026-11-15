// libFuzzer target for the formula lab.
//
// Fuzzed bytes become the lab inputs: the registry index, the scalar arguments,
// the shared series and the text argument. The target asserts only that a
// formula either returns text or raises - never that it crashes. Build with:
//
//   cmake -DBUILD_FUZZERS=ON -DCMAKE_CXX_COMPILER=clang++ ..
//   ./fuzz_lab -runs=100000
//
// GCC has no -fsanitize=fuzzer, so this target is clang-only; the CTest run of
// tests/test_all.cpp and the RapidCheck property already exercise the same
// entry point with generated inputs on every compiler.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "formula_lab.hpp"

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    if (size < 8) return 0;
    const std::size_t count = climate::lab_formulas().size();
    if (count == 0) return 0;
    const std::size_t index = *reinterpret_cast<const std::uint32_t*>(data) % count;

    climate::LabInputs inputs;
    const std::size_t numbers = size / 8;
    for (std::size_t i = 0; i < numbers; ++i) {
        const std::uint8_t* chunk = data + (i * 8) % (size - 7);
        double value = 0.0;
        for (int b = 0; b < 8; ++b)
            value = value * 256.0 + static_cast<double>(chunk[b]);
        inputs.numbers.push_back(value / 1e6);
    }
    const std::size_t points = 4 + (size % 32);
    for (std::size_t i = 0; i < points; ++i) {
        inputs.time.push_back(static_cast<double>(i));
        inputs.series.push_back(static_cast<double>(data[(i * 3) % size]) / 8.0);
    }
    for (std::size_t r = 0; r + 3 < points; r += 2)
        inputs.matrix.emplace_back(inputs.series.begin() + static_cast<long>(r),
                                   inputs.series.begin() + static_cast<long>(r + 3));
    for (std::size_t r = 0; r + 3 < points; ++r)
        inputs.square.push_back({inputs.series[r], inputs.series[r + 1],
                                 inputs.series[r + 2], inputs.series[r + 3]});
    inputs.text.assign(reinterpret_cast<const char*>(data), size < 64 ? size : 64);

    try {
        const std::string out = climate::run_formula(index, inputs);
        (void)out;
    } catch (const std::exception&) {
        // A documented input-domain rejection is a valid outcome.
    }
    return 0;
}
