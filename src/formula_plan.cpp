#include "formula_plan.hpp"

#include <algorithm>

#include "catalog.hpp"
#include "formula_lab.hpp"

namespace climate {
namespace {

bool registry_has(const std::string& id) {
    for (const auto& formula : lab_formulas())
        if (formula.id == id) return true;
    return false;
}

void append(FormulaPlan& plan, const std::vector<PlanStep>& steps) {
    for (const PlanStep& step : steps)
        if (registry_has(step.id)) plan.steps.push_back(step);
}

// Formulas that describe how the variable is retrieved. A CSV series carries no
// radiance or reflectance columns, so these are part of the variable's method but
// are not executed on a series; the Formula Lab runs them given real inputs.
std::vector<PlanStep> method_steps(const std::string& id) {
    if (id == "sst")
        return {{"retrieval::nlsst_features",
                 "split-window NLSST predictors from the 11 and 12 micrometre bands",
                 PlanInput::Sensor}};
    if (id == "ndvi")
        return {{"retrieval::ndvi", "NDVI = (NIR - RED)/(NIR + RED) from surface reflectances",
                 PlanInput::Sensor}};
    if (id == "sea_ice_extent" || id == "sea_ice_concentration")
        return {{"retrieval::polarization_ratio", "PR(19) passive-microwave predictor",
                 PlanInput::Sensor},
                {"retrieval::gradient_ratio", "GR(37V/19V) passive-microwave predictor",
                 PlanInput::Sensor}};
    if (id == "precipitation")
        return {{"retrieval::gpi", "GPI cold-cloud precipitation index", PlanInput::Sensor},
                {"retrieval::gpi_with_water_vapour", "GPI adjusted for precipitable water",
                 PlanInput::Sensor}};
    if (id == "toa_energy_imbalance")
        return {{"retrieval::toa_energy_imbalance", "N = F + lambda dT from forcing and feedback",
                 PlanInput::Sensor}};
    if (id == "co2" || id == "ch4")
        return {{"harmonic::fit_harmonic", "polynomial + harmonic curve fit (Thoning et al.)",
                 PlanInput::Sensor},
                {"harmonic::harmonic_components", "amplitude and phase per harmonic",
                 PlanInput::Sensor},
                {"harmonic::harmonic_seasonal_amplitude", "total seasonal amplitude",
                 PlanInput::Sensor},
                {"harmonic::harmonic_trend_derivative", "growth rate from the fitted trend",
                 PlanInput::Sensor},
                {"retrieval::carbon_growth_rate_pgc_per_year",
                 "growth rate converted to PgC per year", PlanInput::Sensor}};
    return {};
}

// The chain every trend verdict needs, whatever the variable. Each step states
// what it is fed so it receives the series, the time axis or its own scalars.
std::vector<PlanStep> trend_chain() {
    return {
        {"series::interpolate_missing", "fill gaps before decomposing", PlanInput::Series},
        {"series::mean_spacing", "sampling interval used for per-decade rates", PlanInput::Time},
        {"stats::summarize", "level, spread and range of the record", PlanInput::Series},
        {"stats::autocorrelation", "lag-1 persistence the trend test must survive",
         PlanInput::Series},
        {"stats::stdev_value", "noise scale of the record", PlanInput::Series},
        {"stats::median_value", "robust centre, insensitive to outliers", PlanInput::Series},
        {"stats::quantile_value", "distribution of the record", PlanInput::Series},
        {"stats::ols_slope", "least-squares slope as the parametric counterweight",
         PlanInput::TimeSeries},
        {"series::phase_climatology", "calendar climatology of the seasonal cycle",
         PlanInput::Series},
        {"series::remove_seasonal_cycle", "deseasonalised series that is actually tested",
         PlanInput::Series},
        {"stl::stl_decompose", "trend + seasonal + remainder decomposition", PlanInput::Series},
        {"trend::sens_slope", "non-parametric slope with its rank-based interval",
         PlanInput::TimeSeries},
        {"trend::trend_free_prewhitening", "trend-free prewhitened series (Yue et al.)",
         PlanInput::TimeSeries},
        {"trend::trend_ensemble", "accuracy-weighted estimator ensemble", PlanInput::TimeSeries},
        {"power::variance_inflation_factor", "variance inflation from serial correlation",
         PlanInput::Derived},
        {"power::minimum_detectable_trend_power", "smallest trend this record could resolve",
         PlanInput::Derived},
        {"power::trend_power", "power of the observed slope", PlanInput::Derived},
        {"power::variable_is_feasible", "whether the variable can be monitored at this length",
         PlanInput::Derived},
        {"spatial::benjamini_hochberg", "false-discovery control over the tested statistics",
         PlanInput::Series},
        {"spatial::walker_threshold", "field-significance threshold for gridded variables",
         PlanInput::Derived},
        {"uncertainty::root_mean_square_error",
         "residual error against an independent reference record", PlanInput::Sensor},
        {"uncertainty::mean_bias", "bias against an independent reference record",
         PlanInput::Sensor},
    };
}

std::vector<std::string> variable_notes(const std::string& id, const std::string& method) {
    std::vector<std::string> notes;
    if (method == "harmonic")
        notes.push_back("well-mixed greenhouse gas: trend and seasonal cycle come from a "
                        "polynomial plus harmonic fit, not from an STL decomposition");
    else
        notes.push_back("seasonal cycle removed by STL decomposition before the trend test");
    if (method_steps(id).empty())
        notes.push_back("no retrieval formula is registered for this variable; it enters as a "
                        "measured series, so only the trend and significance chain applies");
    else
        notes.push_back("retrieval formulas are included because this variable is derived from "
                        "radiances or reflectances");
    notes.push_back("formulas that need sensor inputs or a reference record are listed as method "
                    "steps and are not executed on the series");
    return notes;
}

}  // namespace

std::size_t FormulaPlan::executable_count() const {
    std::size_t count = 0;
    for (const PlanStep& step : steps)
        if (step.executable()) ++count;
    return count;
}

FormulaPlan plan_for_variable(const std::string& variable_id, int seasonal_period) {
    FormulaPlan plan;
    const std::string id = variable_id.empty() ? std::string("global") : variable_id;
    plan.variable = id;
    plan.method = catalog_analysis_method(id) == AnalysisMethod::HarmonicFit ? "harmonic" : "stl";

    append(plan, method_steps(id));
    append(plan, trend_chain());
    if (seasonal_period >= 2)
        plan.notes.push_back("seasonality assumed at " + std::to_string(seasonal_period) +
                             " samples per cycle");
    for (const std::string& note : variable_notes(id, plan.method)) plan.notes.push_back(note);
    return plan;
}

}  // namespace climate
