#pragma once

#include <string>
#include <vector>

namespace climate {

// How a formula in the plan is fed.
//   series       one vector parameter: the values themselves
//   time         one vector parameter: the time axis
//   time_series  two vector parameters: time, then values
//   derived      scalar parameters only, taken from the record's own statistics
//   sensor       needs radiances, reflectances or a fitted model, which a CSV
//                series cannot supply; listed as part of the variable's method
//                but not run on data it does not have
enum class PlanInput { Series, Time, TimeSeries, Derived, Sensor };

struct PlanStep {
    std::string id;    // registry id, module::name
    std::string role;  // the part of the analysis this formula supplies
    PlanInput input = PlanInput::Series;

    bool executable() const { return input != PlanInput::Sensor; }
};

struct FormulaPlan {
    std::string variable;
    std::string method;              // "harmonic" or "stl"
    std::vector<PlanStep> steps;     // retrieval first, then the trend chain
    std::vector<std::string> notes;  // why this variable gets this chain

    std::size_t executable_count() const;
};

// Chooses the formulas that belong to a variable's analysis from the whole
// registry: the variable's own retrieval formulas, then the seasonality, trend,
// significance, power and uncertainty chain every trend verdict needs. The
// method comes from the catalog, not from the caller. Ids the registry does not
// contain are dropped, so a plan can only name real formulas.
FormulaPlan plan_for_variable(const std::string& variable_id, int seasonal_period = 12);

}  // namespace climate
