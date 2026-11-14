#include "intent.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <sstream>

#include "json.hpp"

namespace climate {
namespace {

std::string lowered(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (char c : text) out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    return out;
}

std::string collapsed(const std::string& text) {
    const std::string lower = lowered(text);
    std::string out;
    bool previous_space = true;
    for (char c : lower) {
        const bool space = std::isspace(static_cast<unsigned char>(c)) || c == '_' || c == '-';
        if (space) {
            if (!previous_space) out.push_back(' ');
        } else {
            out.push_back(c);
        }
        previous_space = space;
    }
    while (!out.empty() && out.back() == ' ') out.pop_back();
    return out;
}

bool contains_phrase(const std::string& haystack, const std::string& phrase) {
    if (phrase.empty()) return false;
    return haystack.find(phrase) != std::string::npos;
}

struct ProcessKeyword {
    const char* phrase;
    const char* process;
};

// Ordered by specificity: multi-word phrases are checked before generic ones.
const ProcessKeyword kProcessKeywords[] = {
    {"sea ice", "sea_ice"},        {"sea-ice", "sea_ice"},
    {"ice extent", "sea_ice"},     {"ice loss", "sea_ice"},
    {"ice cover", "sea_ice"},      {"cryosphere", "cryosphere"},
    {"glacier", "cryosphere"},     {"snow cover", "cryosphere"},
    {"permafrost", "cryosphere"},  {"ndvi", "vegetation"},
    {"vegetation", "vegetation"},  {"greenness", "vegetation"},
    {"greening", "vegetation"},    {"browning", "vegetation"},
    {"leaf area", "vegetation"},   {"drought", "drought"},
    {"precip", "precipitation"},   {"rainfall", "precipitation"},
    {"rain", "precipitation"},     {"rainy", "precipitation"},
    {"co2", "carbon"},             {"co₂", "carbon"},
    {"carbon dioxide", "carbon"},  {"carbon", "carbon"},
    {"methane", "carbon"},         {"ch4", "carbon"},
    {"greenhouse", "greenhouse"},  {"ocean heat", "ocean"},
    {"heat content", "ocean"},     {"ohc", "ocean"},
    {"energy imbalance", "energy_balance"},
    {"energy balance", "energy_balance"},
    {"radiative", "energy_balance"},
    {"albedo", "energy_balance"},  {"forcing", "energy_balance"},
    {"sea level", "ocean"},        {"ocean", "ocean"},
    {"sst", "temperature"},        {"sea surface", "temperature"},
    {"temperature", "temperature"},{"warming", "warming"},
    {"warm", "warming"},           {"hot", "warming"},
    {"changing", "climate_change"},  {"change", "climate_change"},
    {"shifting", "climate_change"},  {"trend", "climate_change"},
};

struct RegionKeyword {
    const char* phrase;
    const char* region;
};

const RegionKeyword kRegionKeywords[] = {
    {"north atlantic", "north_atlantic"}, {"atlantic", "north_atlantic"},
    {"northern hemisphere", "midlatitudes"}, {"mid-latitude", "midlatitudes"},
    {"mid latitude", "midlatitudes"},     {"midlatitudes", "midlatitudes"},
    {"europe", "midlatitudes"},           {"united states", "midlatitudes"},
    {"southern ocean", "antarctic"},      {"antarctic", "antarctic"},
    {"antarctica", "antarctic"},          {"arctic", "arctic"},
    {"north pole", "arctic"},             {"svalbard", "arctic"},
    {"greenland", "arctic"},              {"tropics", "tropics"},
    {"tropical", "tropics"},              {"equator", "tropics"},
    {"amazon", "amazon"},                 {"sahel", "sahel"},
    {"global", "global"},                 {"world", "global"},
    {"earth", "global"},                  {"planet", "global"},
};

bool is_year_token(const std::string& token) {
    if (token.size() != 4) return false;
    for (char c : token)
        if (!std::isdigit(static_cast<unsigned char>(c))) return false;
    const int year = std::atoi(token.c_str());
    return year >= 1850 && year <= 2200;
}

std::vector<std::string> tokenize(const std::string& text) {
    std::vector<std::string> tokens;
    std::string current;
    for (char c : text) {
        if (std::isalnum(static_cast<unsigned char>(c))) {
            current.push_back(c);
        } else if (!current.empty()) {
            tokens.push_back(current);
            current.clear();
        }
    }
    if (!current.empty()) tokens.push_back(current);
    return tokens;
}

int first_matching_process(const std::string& text, std::string& process) {
    for (const auto& keyword : kProcessKeywords) {
        if (contains_phrase(text, keyword.phrase)) {
            process = keyword.process;
            return 1;
        }
    }
    return 0;
}

int first_matching_region(const std::string& text, std::string& region) {
    for (const auto& keyword : kRegionKeywords) {
        if (contains_phrase(text, keyword.phrase)) {
            region = keyword.region;
            return 1;
        }
    }
    return 0;
}

void extract_period(const std::string& text, Intent& intent, int reference_year) {
    const std::vector<std::string> tokens = tokenize(text);
    std::vector<int> years;
    for (const auto& token : tokens) {
        if (is_year_token(token)) years.push_back(std::atoi(token.c_str()));
    }
    if (years.size() >= 2) {
        intent.start_year = std::min(years.front(), years.back());
        intent.end_year = std::max(years.front(), years.back());
        return;
    }
    if (years.size() == 1) {
        const std::string phrase = "since " + std::to_string(years[0]);
        if (contains_phrase(text, phrase) || contains_phrase(text, "after ") ||
            contains_phrase(text, "from ")) {
            intent.start_year = years[0];
            intent.end_year = reference_year;
        } else {
            intent.start_year = years[0];
            intent.end_year = reference_year;
        }
        return;
    }

    // "last N years", "past N decades", "recent N months".
    for (std::size_t i = 0; i + 1 < tokens.size(); ++i) {
        const std::string& token = tokens[i];
        const bool is_number = !token.empty() &&
                               std::all_of(token.begin(), token.end(), [](char c) {
                                   return std::isdigit(static_cast<unsigned char>(c));
                               });
        if (!is_number) continue;
        const int count = std::atoi(token.c_str());
        if (count <= 0) continue;
        const std::string unit = tokens[i + 1];
        int years = 0;
        if (unit == "year" || unit == "years" || unit == "yr" || unit == "yrs")
            years = count;
        else if (unit == "decade" || unit == "decades")
            years = 10 * count;
        else if (unit == "century" || unit == "centuries")
            years = 100 * count;
        if (years > 0) {
            intent.start_year = reference_year - years;
            intent.end_year = reference_year;
            return;
        }
    }

    if (contains_phrase(text, "last decade") || contains_phrase(text, "past decade")) {
        intent.start_year = reference_year - 10;
        intent.end_year = reference_year;
    } else if (contains_phrase(text, "since 1979") || contains_phrase(text, "satellite era")) {
        intent.start_year = 1979;
        intent.end_year = reference_year;
    }
}

}  // namespace

Intent parse_intent_text(const std::string& text, int reference_year) {
    Intent intent;
    intent.source_text = text;
    const std::string normalized = collapsed(text);
    if (normalized.empty()) {
        intent.error = "empty query";
        return intent;
    }

    std::string process;
    if (!first_matching_process(normalized, process)) {
        intent.error =
            "could not resolve a climate process; supported processes are: temperature, warming, "
            "sea_ice, cryosphere, vegetation, drought, precipitation, carbon, ocean, "
            "energy_balance, climate_change";
        return intent;
    }
    intent.process = process;

    std::string region;
    if (!first_matching_region(normalized, region)) {
        region = "global";
        intent.notes.push_back("region not stated; defaulting to global");
    }
    intent.region = region;

    extract_period(normalized, intent, reference_year);
    if (intent.start_year == 0)
        intent.notes.push_back("time period not stated; the full data record will be used");
    if (intent.start_year != 0 && intent.end_year != 0 && intent.start_year > intent.end_year)
        std::swap(intent.start_year, intent.end_year);

    intent.valid = true;
    return intent;
}

Intent parse_intent_json(const std::string& text, int reference_year) {
    Intent intent;
    intent.source_text = text;
    try {
        const json::Value root = json::parse(text);
        if (!root.is_object()) {
            intent.error = "intent JSON must be an object";
            return intent;
        }

        const auto read_string = [&](const char* key) -> std::string {
            const json::Value* value = root.find(key);
            if (value == nullptr || !value->is_string()) return std::string();
            return value->as_string(std::string());
        };

        const std::string raw_process = read_string("process");
        const std::string raw_region = read_string("region");
        const std::string raw_variable = read_string("variable");

        if (!raw_variable.empty()) {
            if (!valid_variable_id(raw_variable)) {
                intent.error = "unknown variable id '" + raw_variable + "'";
                return intent;
            }
            intent.variable = raw_variable;
        }

        if (!raw_process.empty()) {
            const std::string canonical = canonical_process(raw_process);
            if (canonical.empty()) {
                intent.error = "unknown process '" + raw_process + "'";
                return intent;
            }
            intent.process = canonical;
        }

        if (!raw_region.empty()) {
            const std::string canonical = canonical_region(raw_region);
            if (canonical.empty()) {
                intent.error = "unknown region '" + raw_region + "'";
                return intent;
            }
            intent.region = canonical;
        } else {
            intent.region = "global";
        }

        const json::Value* start = root.find("start_year");
        const json::Value* end = root.find("end_year");
        if (start != nullptr && start->is_number())
            intent.start_year = static_cast<int>(start->as_number(0.0));
        if (end != nullptr && end->is_number())
            intent.end_year = static_cast<int>(end->as_number(0.0));
        if (intent.start_year == 0 && intent.end_year == 0) {
            const std::string raw_period = read_string("time_period");
            if (!raw_period.empty()) extract_period(collapsed(raw_period), intent, reference_year);
        }
        if (intent.start_year != 0 && intent.end_year == 0) intent.end_year = reference_year;

        if (intent.process.empty() && intent.variable.empty()) {
            intent.error = "intent JSON must provide either process or variable";
            return intent;
        }
        intent.valid = true;
        return intent;
    } catch (const std::exception& exception) {
        intent.error = std::string("invalid intent JSON: ") + exception.what();
        return intent;
    }
}

std::string intent_to_json(const Intent& intent) {
    json::Value root = json::Value::object();
    root.set("process", json::Value::from_string(intent.process));
    root.set("region", json::Value::from_string(intent.region));
    root.set("start_year", json::Value::from_number(static_cast<double>(intent.start_year)));
    root.set("end_year", json::Value::from_number(static_cast<double>(intent.end_year)));
    if (!intent.variable.empty())
        root.set("variable", json::Value::from_string(intent.variable));
    else
        root.set("variable", json::Value());
    json::Value note_list = json::Value::array();
    for (const auto& note : intent.notes) note_list.push(json::Value::from_string(note));
    root.set("notes", note_list);
    root.set("valid", json::Value::from_bool(intent.valid));
    if (!intent.error.empty())
        root.set("error", json::Value::from_string(intent.error));
    else
        root.set("error", json::Value());
    return json::dump(root, 2);
}

Selection select_variables(const Intent& intent) {
    Selection selection;
    if (!intent.variable.empty()) {
        const CatalogEntry* entry = find_catalog_entry(intent.variable);
        if (entry == nullptr) {
            selection.reasons.push_back("variable '" + intent.variable + "' is not in the catalog");
            return selection;
        }
        selection.entries.push_back(*entry);
        selection.reasons.push_back("explicit variable '" + entry->id + "' selected from the catalog");
        return selection;
    }

    for (const auto& entry : catalog_entries()) {
        // climate_change is a broad intent: every variable that covers the region
        // is relevant, which is how the notes treat "is the Arctic changing?".
        const bool process_match =
            intent.process.empty() || intent.process == "climate_change" ||
            std::find(entry.processes.begin(), entry.processes.end(), intent.process) !=
                entry.processes.end();
        const bool region_match =
            intent.region.empty() ||
            std::find(entry.regions.begin(), entry.regions.end(), intent.region) !=
                entry.regions.end();
        if (process_match && region_match) {
            selection.entries.push_back(entry);
            selection.reasons.push_back(entry.id + ": indicates " + entry.processes.front() +
                                        " in the " + intent.region + " region");
        }
    }
    if (selection.entries.empty())
        selection.reasons.push_back("no catalog variable matches process '" + intent.process +
                                    "' in region '" + intent.region + "'");
    return selection;
}

bool valid_variable_id(const std::string& id) { return find_catalog_entry(id) != nullptr; }

}  // namespace climate
