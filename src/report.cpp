#include "report.hpp"

#include <sstream>

#include "json.hpp"

namespace climate {
namespace {

std::string format_number(double value, int precision = 6) {
    std::ostringstream stream;
    stream.setf(std::ios::fixed);
    stream.precision(precision);
    stream << value;
    return stream.str();
}

std::string format_p(double p) {
    std::ostringstream stream;
    if (p < 1e-12)
        stream << "< 1e-12";
    else if (p < 0.0001)
        stream << "< 0.0001";
    else
        stream << "= " << format_number(p, 4);
    return stream.str();
}

json::Value mk_to_json(const MKResult& mk) {
    json::Value value = json::Value::object();
    value.set("n", json::Value::from_number(static_cast<double>(mk.n)));
    value.set("s", json::Value::from_number(mk.s));
    value.set("variance", json::Value::from_number(mk.variance));
    value.set("variance_corrected", json::Value::from_number(mk.variance_corrected));
    value.set("correction_factor", json::Value::from_number(mk.correction_factor));
    value.set("z", json::Value::from_number(mk.z));
    value.set("p_value", json::Value::from_number(mk.p_value));
    value.set("tau", json::Value::from_number(mk.tau));
    value.set("critical_value", json::Value::from_number(mk.critical_value));
    value.set("significant", json::Value::from_bool(mk.significant));
    value.set("direction", json::Value::from_string(mk.direction));
    return value;
}

json::Value sen_to_json(const SenResult& sen) {
    json::Value value = json::Value::object();
    value.set("slope_per_time_unit", json::Value::from_number(sen.slope));
    value.set("intercept", json::Value::from_number(sen.intercept));
    value.set("lower", json::Value::from_number(sen.lower));
    value.set("upper", json::Value::from_number(sen.upper));
    value.set("alpha", json::Value::from_number(sen.alpha));
    value.set("pairs", json::Value::from_number(static_cast<double>(sen.pairs)));
    value.set("computed", json::Value::from_bool(sen.computed));
    return value;
}

json::Value strings_to_json(const std::vector<std::string>& values) {
    json::Value list = json::Value::array();
    for (const auto& value : values) list.push(json::Value::from_string(value));
    return list;
}

}  // namespace

std::string series_report_text(const SeriesAnalysis& analysis, const std::string& variable,
                               const std::string& units, const std::string& source) {
    std::ostringstream out;
    const std::string label = analysis.label.empty() ? variable : analysis.label;
    out << "Variable          : " << label << " (" << variable << ") [" << units << "]\n";
    out << "Source            : " << source << "\n";
    out << "Observations      : " << analysis.n;
    if (analysis.missing_filled > 0)
        out << " (" << analysis.missing_filled << " missing values interpolated)";
    out << "\n";
    out << "Period            : " << format_number(analysis.start_time, 3) << " to "
        << format_number(analysis.end_time, 3) << " (mean spacing "
        << format_number(analysis.spacing, 5) << ")\n";
    out << "Mean / SD         : " << format_number(analysis.summary.mean, 4) << " / "
        << format_number(analysis.summary.sd, 4) << "\n";
    out << "Median / IQR      : " << format_number(analysis.summary.median, 4) << " / ["
        << format_number(analysis.summary.q25, 4) << ", "
        << format_number(analysis.summary.q75, 4) << "]\n";
    out << "Minimum / maximum : " << format_number(analysis.summary.minimum, 4) << " / "
        << format_number(analysis.summary.maximum, 4) << "\n\n";

    out << "Decomposition\n";
    if (analysis.decomposed) {
        out << "  Method              : STL (period " << analysis.period
            << ", robustness iterations enabled)\n";
        out << "  Seasonal amplitude  : " << format_number(analysis.seasonal_amplitude, 5) << " "
            << units << "\n";
        out << "  Seasonal strength   : " << format_number(analysis.seasonal_strength, 3) << "\n";
        out << "  Trend strength      : " << format_number(analysis.trend_strength, 3) << "\n";
    } else {
        out << "  Method              : calendar-month climatology (period " << analysis.period
            << ")\n";
        out << "  Seasonal amplitude  : " << format_number(analysis.seasonal_amplitude, 5) << " "
            << units << "\n";
    }
    out << "  Standardised anomaly: latest " << format_number(analysis.latest_anomaly, 3)
        << " sigma, range [" << format_number(analysis.minimum_anomaly, 3) << ", "
        << format_number(analysis.maximum_anomaly, 3) << "]\n\n";

    out << "Trend test on " << (analysis.prewhiten ? "prewhitened " : "")
        << "deseasonalised series\n";
    out << "  Mann-Kendall S      : " << format_number(analysis.mk.s, 2) << ", Z = "
        << format_number(analysis.mk.z, 4) << ", p " << format_p(analysis.mk.p_value)
        << ", tau = " << format_number(analysis.mk.tau, 4) << "\n";
    out << "  Variance correction : x" << format_number(analysis.mk.correction_factor, 4)
        << (analysis.hamed_rao ? " (Hamed-Rao, serial correlation)" : " (none)") << "\n";
    out << "  Sen slope           : " << format_number(analysis.sen.slope, 6) << " " << units
        << "/time unit, 95% CI [" << format_number(analysis.sen.lower, 6) << ", "
        << format_number(analysis.sen.upper, 6) << "]\n";
    out << "                        " << format_number(analysis.sen_per_decade, 5) << " " << units
        << "/decade from " << analysis.sen.pairs << " pairwise slopes\n";
    out << "  Result              : " << analysis.mk.direction << " (";
    out << (analysis.mk.significant ? "significant" : "not significant")
        << " at alpha = " << format_number(analysis.alpha, 3) << ")\n\n";

    out << "Cross-check on the raw series\n";
    out << "  Mann-Kendall Z      : " << format_number(analysis.mk_raw.z, 4) << ", p "
        << format_p(analysis.mk_raw.p_value) << "\n";
    out << "  Sen slope           : " << format_number(analysis.sen_raw.slope, 6) << " " << units
        << "/time unit (" << format_number(analysis.sen_raw_per_decade, 5) << " " << units
        << "/decade)\n";
    if (analysis.period >= 2) {
        out << "\nSeasonal (within-block) Mann-Kendall\n";
        out << "  S = " << format_number(analysis.mk_seasonal.s, 2) << ", Z = "
            << format_number(analysis.mk_seasonal.z, 4) << ", p "
            << format_p(analysis.mk_seasonal.p_value) << ", " << analysis.mk_seasonal.direction
            << "\n";
    }

    if (!analysis.warnings.empty()) {
        out << "\nNotes\n";
        for (const auto& warning : analysis.warnings) out << "  - " << warning << "\n";
    }
    return out.str();
}

std::string series_report_markdown(const SeriesAnalysis& analysis, const std::string& variable,
                                   const std::string& units, const std::string& source) {
    std::ostringstream out;
    const std::string label = analysis.label.empty() ? variable : analysis.label;
    out << "# " << label << " (" << variable << ")\n\n";
    out << "| Field | Value |\n|---|---|\n";
    out << "| Units | " << units << " |\n";
    out << "| Source | `" << source << "` |\n";
    out << "| Observations | " << analysis.n << " |\n";
    out << "| Period | " << format_number(analysis.start_time, 3) << " to "
        << format_number(analysis.end_time, 3) << " |\n";
    out << "| Mean | " << format_number(analysis.summary.mean, 4) << " |\n";
    out << "| Standard deviation | " << format_number(analysis.summary.sd, 4) << " |\n";
    out << "| Seasonal amplitude | " << format_number(analysis.seasonal_amplitude, 5) << " |\n";
    out << "| Seasonal strength | " << format_number(analysis.seasonal_strength, 3) << " |\n";
    out << "| Trend strength | " << format_number(analysis.trend_strength, 3) << " |\n";
    out << "| Mann-Kendall Z | " << format_number(analysis.mk.z, 4) << " |\n";
    out << "| Mann-Kendall p | " << format_p(analysis.mk.p_value) << " |\n";
    out << "| Sen slope | " << format_number(analysis.sen.slope, 6) << " " << units
        << "/time unit |\n";
    out << "| Sen slope per decade | " << format_number(analysis.sen_per_decade, 5) << " " << units
        << "/decade |\n";
    out << "| Sen 95% CI | [" << format_number(analysis.sen.lower, 6) << ", "
        << format_number(analysis.sen.upper, 6) << "] |\n";
    out << "| Verdict | " << analysis.mk.direction << " |\n\n";
    if (!analysis.warnings.empty()) {
        out << "## Notes\n\n";
        for (const auto& warning : analysis.warnings) out << "- " << warning << "\n";
    }
    return out.str();
}

std::string series_report_json(const SeriesAnalysis& analysis, const std::string& variable,
                               const std::string& units, const std::string& source) {
    json::Value root = json::Value::object();
    root.set("variable", json::Value::from_string(variable));
    root.set("label", json::Value::from_string(analysis.label));
    root.set("units", json::Value::from_string(units));
    root.set("source", json::Value::from_string(source));
    root.set("samples", json::Value::from_number(static_cast<double>(analysis.n)));
    root.set("missing_filled", json::Value::from_number(static_cast<double>(analysis.missing_filled)));
    root.set("start_time", json::Value::from_number(analysis.start_time));
    root.set("end_time", json::Value::from_number(analysis.end_time));
    root.set("spacing", json::Value::from_number(analysis.spacing));
    root.set("period", json::Value::from_number(static_cast<double>(analysis.period)));
    root.set("alpha", json::Value::from_number(analysis.alpha));
    root.set("hamed_rao", json::Value::from_bool(analysis.hamed_rao));
    root.set("prewhitened", json::Value::from_bool(analysis.prewhiten));

    json::Value summary = json::Value::object();
    summary.set("mean", json::Value::from_number(analysis.summary.mean));
    summary.set("sd", json::Value::from_number(analysis.summary.sd));
    summary.set("median", json::Value::from_number(analysis.summary.median));
    summary.set("minimum", json::Value::from_number(analysis.summary.minimum));
    summary.set("maximum", json::Value::from_number(analysis.summary.maximum));
    summary.set("q25", json::Value::from_number(analysis.summary.q25));
    summary.set("q75", json::Value::from_number(analysis.summary.q75));
    root.set("summary", summary);

    json::Value decomposition = json::Value::object();
    decomposition.set("method", json::Value::from_string(analysis.decomposed ? "STL" : "climatology"));
    decomposition.set("seasonal_amplitude", json::Value::from_number(analysis.seasonal_amplitude));
    decomposition.set("seasonal_strength", json::Value::from_number(analysis.seasonal_strength));
    decomposition.set("trend_strength", json::Value::from_number(analysis.trend_strength));
    decomposition.set("latest_anomaly_sigma", json::Value::from_number(analysis.latest_anomaly));
    decomposition.set("minimum_anomaly_sigma", json::Value::from_number(analysis.minimum_anomaly));
    decomposition.set("maximum_anomaly_sigma", json::Value::from_number(analysis.maximum_anomaly));
    root.set("decomposition", decomposition);

    json::Value trend = json::Value::object();
    trend.set("series", json::Value::from_string("deseasonalised"));
    trend.set("mann_kendall", mk_to_json(analysis.mk));
    trend.set("sen", sen_to_json(analysis.sen));
    trend.set("sen_per_decade", json::Value::from_number(analysis.sen_per_decade));
    root.set("trend", trend);

    json::Value raw = json::Value::object();
    raw.set("series", json::Value::from_string("raw"));
    raw.set("mann_kendall", mk_to_json(analysis.mk_raw));
    raw.set("sen", sen_to_json(analysis.sen_raw));
    raw.set("sen_per_decade", json::Value::from_number(analysis.sen_raw_per_decade));
    root.set("raw_trend", raw);

    if (analysis.period >= 2) root.set("seasonal_mann_kendall", mk_to_json(analysis.mk_seasonal));
    root.set("warnings", strings_to_json(analysis.warnings));
    return json::dump(root, 2);
}

std::string co2_report_text(const CO2Analysis& analysis, const std::string& variable,
                            const std::string& units, const std::string& source) {
    std::ostringstream out;
    const std::string label = analysis.label.empty() ? variable : analysis.label;
    out << "Variable          : " << label << " (" << variable << ") [" << units << "]\n";
    out << "Source            : " << source << "\n";
    out << "Observations      : " << analysis.n;
    if (analysis.missing_filled > 0)
        out << " (" << analysis.missing_filled << " filled)";
    out << "\n";
    out << "Period            : " << format_number(analysis.start_time, 3) << " to "
        << format_number(analysis.end_time, 3) << "\n";
    out << "Model             : polynomial degree " << analysis.model.polynomial_order << " + "
        << analysis.model.harmonic_count << " harmonics\n";
    out << "Fit               : R^2 = " << format_number(analysis.model.r_squared, 6)
        << ", residual SD = " << format_number(analysis.model.residual_sd, 5)
        << ", used " << analysis.model.used_points << " points, rejected "
        << analysis.model.rejected_points << "\n";
    out << "Seasonal amplitude: " << format_number(analysis.seasonal_amplitude, 5) << " " << units
        << "\n\n";

    out << "Harmonic components\n";
    for (const auto& component : analysis.components) {
        out << "  n=" << component.number << ": amplitude "
            << format_number(component.amplitude, 5) << " " << units << ", phase "
            << format_number(component.phase_degrees, 3) << " deg\n";
    }

    out << "\nTrend and growth rate\n";
    out << "  Polynomial slope  : " << format_number(analysis.linear_coefficient, 6) << " " << units
        << "/year at record centre\n";
    out << "  Quadratic term    : " << format_number(analysis.quadratic_coefficient, 8) << " "
        << units << "/year^2\n";
    out << "  Growth at start   : " << format_number(analysis.growth_rate_start, 4) << " " << units
        << "/year\n";
    out << "  Growth at end     : " << format_number(analysis.growth_rate_end, 4) << " " << units
        << "/year (" << format_number(analysis.growth_rate_end_pgc_per_year, 4) << " PgC/year)\n";
    out << "  Mean growth       : " << format_number(analysis.growth_rate_mean, 4) << " " << units
        << "/year\n";
    out << "  Growth last 5 yr  : " << format_number(analysis.growth_rate_last_five_years, 4) << " "
        << units << "/year\n";
    out << "\nDeseasonalised trend test\n";
    out << "  Mann-Kendall Z    : " << format_number(analysis.mk.z, 4) << ", p "
        << format_p(analysis.mk.p_value) << ", tau = " << format_number(analysis.mk.tau, 4) << "\n";
    out << "  Sen slope         : " << format_number(analysis.sen.slope, 6) << " " << units
        << "/year [" << format_number(analysis.sen.lower, 6) << ", "
        << format_number(analysis.sen.upper, 6) << "]\n";
    out << "  Result            : " << analysis.mk.direction << "\n";

    if (!analysis.warnings.empty()) {
        out << "\nNotes\n";
        for (const auto& warning : analysis.warnings) out << "  - " << warning << "\n";
    }
    return out.str();
}

std::string co2_report_json(const CO2Analysis& analysis, const std::string& variable,
                            const std::string& units, const std::string& source) {
    json::Value root = json::Value::object();
    root.set("variable", json::Value::from_string(variable));
    root.set("label", json::Value::from_string(analysis.label));
    root.set("units", json::Value::from_string(units));
    root.set("source", json::Value::from_string(source));
    root.set("samples", json::Value::from_number(static_cast<double>(analysis.n)));
    root.set("missing_filled", json::Value::from_number(static_cast<double>(analysis.missing_filled)));
    root.set("start_time", json::Value::from_number(analysis.start_time));
    root.set("end_time", json::Value::from_number(analysis.end_time));

    json::Value model = json::Value::object();
    model.set("polynomial_order", json::Value::from_number(analysis.model.polynomial_order));
    model.set("harmonic_count", json::Value::from_number(analysis.model.harmonic_count));
    model.set("period", json::Value::from_number(analysis.model.period));
    model.set("r_squared", json::Value::from_number(analysis.model.r_squared));
    model.set("residual_sd", json::Value::from_number(analysis.model.residual_sd));
    model.set("used_points", json::Value::from_number(static_cast<double>(analysis.model.used_points)));
    model.set("rejected_points",
              json::Value::from_number(static_cast<double>(analysis.model.rejected_points)));
    root.set("model", model);

    json::Value components = json::Value::array();
    for (const auto& component : analysis.components) {
        json::Value item = json::Value::object();
        item.set("number", json::Value::from_number(component.number));
        item.set("amplitude", json::Value::from_number(component.amplitude));
        item.set("phase_degrees", json::Value::from_number(component.phase_degrees));
        components.push(item);
    }
    root.set("harmonics", components);
    root.set("seasonal_amplitude", json::Value::from_number(analysis.seasonal_amplitude));

    json::Value growth = json::Value::object();
    growth.set("linear_coefficient", json::Value::from_number(analysis.linear_coefficient));
    growth.set("quadratic_coefficient", json::Value::from_number(analysis.quadratic_coefficient));
    growth.set("start_per_year", json::Value::from_number(analysis.growth_rate_start));
    growth.set("end_per_year", json::Value::from_number(analysis.growth_rate_end));
    growth.set("mean_per_year", json::Value::from_number(analysis.growth_rate_mean));
    growth.set("last_five_years_per_year",
               json::Value::from_number(analysis.growth_rate_last_five_years));
    growth.set("end_pgc_per_year", json::Value::from_number(analysis.growth_rate_end_pgc_per_year));
    root.set("growth_rate", growth);

    json::Value trend = json::Value::object();
    trend.set("mann_kendall", mk_to_json(analysis.mk));
    trend.set("sen", sen_to_json(analysis.sen));
    root.set("trend", trend);
    root.set("warnings", strings_to_json(analysis.warnings));
    return json::dump(root, 2);
}

}  // namespace climate
