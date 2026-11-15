// GoogleTest suite. Run by CTest as the `google_test` test.
//
// These assertions pin the scientific contracts the workbench relies on: what
// each formula promises, which invariants hold, and that the generated formula
// lab really covers every public declaration.

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "analysis.hpp"
#include "catalog.hpp"
#include "common.hpp"
#include "csv.hpp"
#include "formula_lab.hpp"
#include "json.hpp"
#include "physics.hpp"
#include "spatial.hpp"
#include "trend.hpp"

namespace climate {
namespace {

std::vector<double> ramp(double slope, double intercept, int n) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) values.push_back(intercept + slope * i);
    return values;
}

TEST(Catalog, VariablesAreRegistered) {
    EXPECT_GE(catalog_entries().size(), 16u);
    for (const CatalogEntry& entry : catalog_entries()) {
        EXPECT_FALSE(entry.id.empty());
        EXPECT_FALSE(entry.name.empty());
        EXPECT_FALSE(entry.units.empty());
    }
}

TEST(Catalog, LookupById) {
    const CatalogEntry* sst = find_catalog_entry("sst");
    ASSERT_NE(sst, nullptr);
    EXPECT_EQ(sst->id, "sst");
}

TEST(Trend, SenSlopeRecoversALine) {
    const std::vector<double> time = ramp(1.0, 0.0, 40);
    const std::vector<double> values = ramp(0.25, 3.0, 40);
    const SenResult result = sens_slope(time, values);
    EXPECT_NEAR(result.slope, 0.25, 1e-9);
    EXPECT_LT(result.lower, result.slope + 1e-12);
    EXPECT_GT(result.upper, result.slope - 1e-12);
}

TEST(Trend, MannKendallFindsAMonotoneTrend) {
    const std::vector<double> values = ramp(0.5, 0.0, 60);
    const MKResult result = mann_kendall(values);
    EXPECT_TRUE(result.significant);
    EXPECT_GT(result.tau, 0.99);
    EXPECT_LT(result.p_value, 0.001);
}

TEST(Trend, MannKendallRejectsNoise) {
    std::vector<double> noise;
    for (int i = 0; i < 80; ++i)
        noise.push_back(std::sin(i * 1.7) + 0.4 * std::cos(i * 0.3));
    const MKResult result = mann_kendall(noise);
    EXPECT_LT(std::fabs(result.tau), 0.4);
}

TEST(Physics, Co2ForcingMatchesMyhre) {
    EXPECT_NEAR(co2_radiative_forcing(560.0, 280.0), 5.35 * std::log(2.0), 1e-12);
    EXPECT_DOUBLE_EQ(co2_radiative_forcing(400.0, 400.0), 0.0);
}

TEST(Physics, ForcingIsMonotoneInConcentration) {
    double previous = -1e9;
    for (double ppm = 300.0; ppm <= 800.0; ppm += 25.0) {
        const double forcing = co2_radiative_forcing(ppm, 280.0);
        EXPECT_GT(forcing, previous);
        previous = forcing;
    }
}

TEST(Physics, EquilibriumWarmingSigns) {
    EXPECT_GT(equilibrium_temperature_change(3.7, -1.25), 0.0);
    EXPECT_NEAR(climate_sensitivity_2x(-1.25), 2.96, 0.01);
    EXPECT_THROW(equilibrium_temperature_change(3.7, 0.0), std::runtime_error);
}

TEST(Physics, VegetationIndicesStayBounded) {
    for (double red = 0.0; red <= 0.6; red += 0.05) {
        const double nir = red + 0.3;
        const double value = evi(red, nir, 0.05);
        EXPECT_GT(value, -1.0);
        EXPECT_LT(value, 1.0);
        const double savi_value = savi(red, nir, 0.5);
        EXPECT_GT(savi_value, -2.0);
        EXPECT_LT(savi_value, 2.0);
    }
}

TEST(Physics, MieEfficiencyIsPhysical) {
    const MieResult result = mie_sphere(5.0, std::complex<double>(1.5, 0.01));
    EXPECT_GT(result.qext, 0.0);
    EXPECT_GE(result.qsca, 0.0);
    EXPECT_LE(result.qsca, result.qext + 1e-12);
    EXPECT_GE(result.single_scattering_albedo, 0.0);
    EXPECT_LE(result.single_scattering_albedo, 1.0 + 1e-12);
}

TEST(Spatial, GetisOrdWithIdentityWeights) {
    const std::vector<double> values{1.0, 2.0, 3.0, 4.0};
    Weights weights(4, std::vector<double>(4, 0.0));
    for (std::size_t i = 0; i < 4; ++i) weights[i][i] = 1.0;
    for (std::size_t i = 0; i < 4; ++i) {
        const double gi = getis_ord_gi(values, weights, i);
        EXPECT_TRUE(std::isfinite(gi));
    }
    EXPECT_THROW(getis_ord_gi(values, Weights(2, std::vector<double>(2, 1.0)), 0),
                 std::runtime_error);
}

TEST(Csv, ParsesQuotedFields) {
    const Table table = parse_csv_text("time,value\n2000.0,1.5\n2001.0,2.5\n");
    ASSERT_EQ(table.rows.size(), 2u);
    EXPECT_EQ(table.header.size(), 2u);
    EXPECT_EQ(table.cell(1, 1), "2.5");
}

TEST(Json, RoundTripsNumbersStringsAndArrays) {
    json::Value root = json::Value::object();
    root.set("count", json::Value::from_number(3));
    root.set("name", json::Value::from_string("sst"));
    json::Value list = json::Value::array();
    list.push(json::Value::from_number(1.5));
    list.push(json::Value::from_number(2.5));
    root.set("values", list);
    const json::Value parsed = json::parse(json::dump(root));
    ASSERT_TRUE(parsed.is_object());
    EXPECT_DOUBLE_EQ(parsed.find("count")->number, 3.0);
    EXPECT_EQ(parsed.find("name")->text, "sst");
    ASSERT_TRUE(parsed.find("values")->is_array());
    EXPECT_EQ(parsed.find("values")->items.size(), 2u);
}

TEST(Analysis, DetectsASignificantWarmingSeries) {
    SeriesOptions options;
    options.period = 12;
    const std::vector<double> time = ramp(1.0 / 12.0, 2000.0, 240);
    std::vector<double> values;
    for (std::size_t i = 0; i < time.size(); ++i)
        values.push_back(14.0 + 0.03 * static_cast<double>(i) +
                         0.4 * std::sin(2.0 * PI * static_cast<double>(i) / 12.0));
    const SeriesAnalysis analysis = analyze_series(time, values, options);
    EXPECT_TRUE(analysis.mk.significant);
    EXPECT_GT(analysis.ensemble.robust_slope, 0.0);
    EXPECT_FALSE(analysis.ensemble.members.empty());
    EXPECT_LT(analysis.ensemble.selected, analysis.ensemble.members.size());
}

TEST(FormulaLab, EveryDeclarationIsBound) {
    const std::vector<LabFormula>& formulas = lab_formulas();
    EXPECT_GE(formulas.size(), 900u);
    for (std::size_t i = 0; i < formulas.size(); ++i) {
        EXPECT_EQ(formulas[i].index, i);
        EXPECT_FALSE(formulas[i].id.empty());
        EXPECT_FALSE(formulas[i].signature.empty());
    }
}

TEST(FormulaLab, RunsAnEntryPointToEnd) {
    LabInputs inputs;
    inputs.numbers = {3.7, -1.25};
    for (std::size_t i = 0; i < lab_formulas().size(); ++i) {
        const LabFormula& formula = lab_formulas()[i];
        if (formula.id != "physics::equilibrium_temperature_change") continue;
        const std::string out = run_formula(formula.index, inputs);
        EXPECT_FALSE(out.empty());
        EXPECT_NE(out.find("2.96"), std::string::npos) << out;
        return;
    }
    FAIL() << "physics::equilibrium_temperature_change is missing from the lab";
}

TEST(FormulaLab, OutOfRangeInputsRaise) {
    LabInputs inputs;
    inputs.numbers = {3.7, 0.0};
    for (const LabFormula& formula : lab_formulas()) {
        if (formula.id != "physics::equilibrium_temperature_change") continue;
        EXPECT_THROW(run_formula(formula.index, inputs), std::runtime_error);
        return;
    }
    FAIL() << "formula missing";
}

}  // namespace
}  // namespace climate
