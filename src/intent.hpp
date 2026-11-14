#pragma once

#include <string>
#include <vector>

#include "catalog.hpp"

namespace climate {

struct Intent {
    std::string process;    // canonical process name, or empty when unresolved
    std::string region;     // canonical region name, or empty when unresolved
    int start_year = 0;     // 0 means unspecified
    int end_year = 0;       // 0 means unspecified
    std::string variable;   // optional explicit catalog id
    std::string source_text;
    std::vector<std::string> notes;
    bool valid = false;
    std::string error;
};

struct Selection {
    std::vector<CatalogEntry> entries;
    std::vector<std::string> reasons;
};

// Deterministic intent parser used when no external language model provides the
// intermediate representation. It classifies process, region and time period
// from free text, mirroring the schema an LLM would emit.
Intent parse_intent_text(const std::string& text, int reference_year = 2025);

// Ingests the JSON intermediate representation produced by a language model
// ({"process": ..., "region": ..., "start_year": ..., "end_year": ...,
// "variable": ...}) and validates every field against the catalog whitelist.
Intent parse_intent_json(const std::string& text, int reference_year = 2025);

std::string intent_to_json(const Intent& intent);

// Deterministic compiler: resolves the intent into catalog variables. Unknown
// identifiers are rejected rather than silently ignored.
Selection select_variables(const Intent& intent);

bool valid_variable_id(const std::string& id);

}  // namespace climate
