// RapidCheck property tests. Run by CTest as the `rapidcheck` test.
//
// Each property states a contract that must hold for every generated input, not
// just for the hand-picked examples in the unit suite. `property()` propagates
// rc::check's verdict so a falsified property fails the test.

#include <rapidcheck.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "common.hpp"
#include "formula_lab.hpp"
#include "json.hpp"
#include "physics.hpp"
#include "trend.hpp"

using namespace climate;

namespace {

int failures = 0;

template <typename Fn>
void property(const std::string& name, Fn fn) {
    if (!rc::check(name, fn)) ++failures;
}

}  // namespace

int main() {
    property("Sen slope recovers any straight line", [] {
        const double slope = *rc::gen::inRange(-50, 50) / 100.0;
        const double intercept = *rc::gen::inRange(-100, 100) / 10.0;
        const int n = *rc::gen::inRange(8, 60);
        std::vector<double> time, values;
        for (int i = 0; i < n; ++i) {
            time.push_back(static_cast<double>(i));
            values.push_back(intercept + slope * static_cast<double>(i));
        }
        const SenResult result = sens_slope(time, values);
        RC_ASSERT(std::fabs(result.slope - slope) < 1e-9);
    });

    property("increasing series never has a negative Sen slope", [] {
        const int n = *rc::gen::inRange(6, 40);
        std::vector<double> time, values;
        double value = 0.0;
        for (int i = 0; i < n; ++i) {
            time.push_back(static_cast<double>(i));
            value += *rc::gen::inRange(1, 100) / 100.0;
            values.push_back(value);
        }
        RC_ASSERT(sens_slope(time, values).slope > 0.0);
    });

    property("CO2 forcing reverses when the ratio is inverted", [] {
        const double reference = *rc::gen::inRange(200, 400);
        const double a = *rc::gen::inRange(200, 1000);
        RC_ASSERT(std::fabs(co2_radiative_forcing(a, reference) +
                            co2_radiative_forcing(reference, a)) < 1e-12);
        RC_ASSERT(std::fabs(co2_radiative_forcing(reference, reference)) < 1e-12);
    });

    property("CO2 forcing increases with concentration", [] {
        const double reference = *rc::gen::inRange(180, 300);
        const double low = *rc::gen::inRange(300, 600);
        const double high = low + *rc::gen::inRange(1, 400);
        RC_ASSERT(co2_radiative_forcing(low, reference) <
                  co2_radiative_forcing(high, reference));
    });

    property("vegetation indices stay bounded for physical reflectances", [] {
        // Physical ordering: blue <= red <= nir, all inside [0, 1].
        const int nir_hundredths = *rc::gen::inRange(5, 100);
        const int red_hundredths = *rc::gen::inRange(0, nir_hundredths);
        const int blue_hundredths = *rc::gen::inRange(0, red_hundredths);
        const double nir = nir_hundredths / 100.0;
        const double red = red_hundredths / 100.0;
        const double blue = blue_hundredths / 100.0;
        const double denominator = nir + 6.0 * red - 7.5 * blue + 1.0;
        RC_ASSERT(denominator > 0.0);
        // For blue <= red <= nir the denominator stays above ~0.97, so EVI's
        // mathematical bound is 2.5 / 0.97; it only approaches [-1, 1] for
        // mid-range reflectances, which the sign check below captures.
        const double value = evi(red, nir, blue);
        RC_ASSERT(std::isfinite(value));
        RC_ASSERT(value >= -2.6 && value <= 2.6);
        RC_ASSERT((value >= 0.0) == (nir >= red));
        const double savi_value = savi(red, nir, 0.5);
        RC_ASSERT(savi_value >= -1.0 && savi_value <= 1.0);
        RC_ASSERT((savi_value >= 0.0) == (nir >= red));
    });

    property("JSON round-trip preserves every number", [] {
        const double a = *rc::gen::inRange(-1000000, 1000000) / 1000.0;
        json::Value root = json::Value::object();
        root.set("a", json::Value::from_number(a));
        json::Value list = json::Value::array();
        list.push(json::Value::from_number(a));
        list.push(json::Value::from_number(-a));
        root.set("list", list);
        const json::Value parsed = json::parse(json::dump(root));
        RC_ASSERT(parsed.find("a")->number == a);
        RC_ASSERT(parsed.find("list")->items.size() == 2u);
        RC_ASSERT(parsed.find("list")->items[1].number == -a);
    });

    property("the lab always returns text or raises, never crashes", [] {
        LabInputs inputs;
        const double a = *rc::gen::inRange(-100, 100);
        inputs.numbers = {a, a + 1.0, a * 0.5};
        inputs.time = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        inputs.series = {1.0, 2.0, 1.5, 3.0, 2.5, 4.0};
        inputs.matrix = {{1.0, 2.0}, {3.0, 4.0}};
        inputs.square = {{1.0, 0.0}, {0.0, 1.0}};
        inputs.text = "property";
        const std::size_t count = lab_formulas().size();
        const std::size_t index = *rc::gen::inRange<std::size_t>(0, count);
        try {
            const std::string out = run_formula(index, inputs);
            RC_ASSERT(!out.empty());
        } catch (const std::exception&) {
            RC_SUCCEED("formula rejected the generated inputs, which is allowed");
        }
    });

    property("equilibrium warming is positive for stabilising feedbacks", [] {
        const double forcing = *rc::gen::inRange(1, 100) / 10.0;
        const double feedback = -*rc::gen::inRange(1, 50) / 10.0;
        RC_ASSERT(equilibrium_temperature_change(forcing, feedback) >= 0.0);
    });

    if (failures == 0) {
        std::printf("rapidcheck: all properties hold\n");
        return 0;
    }
    std::printf("rapidcheck: %d propert%s failed\n", failures, failures == 1 ? "y" : "ies");
    return 1;
}
