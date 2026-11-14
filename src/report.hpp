#pragma once

#include <string>
#include <vector>

#include "analysis.hpp"

namespace climate {

std::string series_report_text(const SeriesAnalysis& analysis, const std::string& variable,
                               const std::string& units, const std::string& source);
std::string series_report_markdown(const SeriesAnalysis& analysis, const std::string& variable,
                                   const std::string& units, const std::string& source);
std::string series_report_json(const SeriesAnalysis& analysis, const std::string& variable,
                               const std::string& units, const std::string& source);

std::string co2_report_text(const CO2Analysis& analysis, const std::string& variable,
                            const std::string& units, const std::string& source);
std::string co2_report_json(const CO2Analysis& analysis, const std::string& variable,
                            const std::string& units, const std::string& source);

}  // namespace climate
