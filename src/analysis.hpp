#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "harmonic.hpp"
#include "stats.hpp"
#include "stl.hpp"
#include "trend.hpp"

namespace climate {

struct SeriesOptions {
    int period = 12;
    int seasonal_window = 13;
    int trend_window = 25;
    int low_pass_window = 13;
    int seasonal_degree = 1;
    int trend_degree = 1;
    int inner_iterations = 2;
    int outer_iterations = 5;
    double alpha = 0.05;
    double time_units_per_year = 1.0;  // 1.0 when time is expressed in years
    bool hamed_rao = true;
    bool prewhiten = false;
    bool test_deseasonalised = true;
};

struct SeriesAnalysis {
    std::string label;
    std::size_t n = 0;
    std::size_t missing_filled = 0;
    double start_time = 0.0;
    double end_time = 0.0;
    double spacing = 0.0;
    int period = 12;
    double alpha = 0.05;
    bool hamed_rao = true;
    bool prewhiten = false;

    Summary summary;

    bool decomposed = false;
    STLResult decomposition;
    std::vector<double> deseasonalised;
    std::vector<double> standardised_anomalies;
    double seasonal_amplitude = 0.0;
    double seasonal_strength = 0.0;
    double trend_strength = 0.0;

    MKResult mk;       // on the series used for trend testing
    SenResult sen;
    double sen_per_decade = 0.0;
    MKResult mk_raw;
    SenResult sen_raw;
    double sen_raw_per_decade = 0.0;
    MKResult mk_seasonal;
    TrendEnsemble ensemble;

    double latest_anomaly = 0.0;
    double minimum_anomaly = 0.0;
    double maximum_anomaly = 0.0;

    std::vector<std::string> warnings;
};

SeriesAnalysis analyze_series(const std::vector<double>& time, const std::vector<double>& values,
                              const SeriesOptions& options = {});

struct CO2Options {
    int polynomial_order = 2;
    int harmonic_count = 4;
    double period = 1.0;  // one year when time is fractional years
    int reject_iterations = 4;
    double reject_threshold = 3.0;
    double alpha = 0.05;
    double time_units_per_year = 1.0;
};

struct CO2Analysis {
    std::string label;
    std::size_t n = 0;
    std::size_t missing_filled = 0;
    double start_time = 0.0;
    double end_time = 0.0;

    HarmonicModel model;
    std::vector<double> fitted;
    std::vector<double> residual;
    std::vector<HarmonicComponent> components;
    double seasonal_amplitude = 0.0;
    double growth_rate_start = 0.0;
    double growth_rate_end = 0.0;
    double growth_rate_mean = 0.0;
    double growth_rate_last_five_years = 0.0;
    double growth_rate_end_pgc_per_year = 0.0;
    double linear_coefficient = 0.0;
    double quadratic_coefficient = 0.0;

    MKResult mk;
    SenResult sen;

    std::vector<std::string> warnings;
};

CO2Analysis analyze_co2(const std::vector<double>& time, const std::vector<double>& values,
                        const CO2Options& options = {});

}  // namespace climate
