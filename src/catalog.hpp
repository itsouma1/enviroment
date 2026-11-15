#pragma once

#include <string>
#include <vector>

namespace climate {

struct CatalogEntry {
    std::string id;
    std::string name;
    std::string units;
    std::vector<std::string> processes;
    std::vector<std::string> regions;
    std::string retrieval;        // how the variable is derived from radiances / reflectances
    std::string analysis_method;  // how the trend pipeline treats it
    std::string description;
};

// How the catalog says a variable is analysed. The engine chooses from this
// instead of applying one default method to every variable.
enum class AnalysisMethod { SeasonalStl, HarmonicFit };

AnalysisMethod catalog_analysis_method(const std::string& id);

const std::vector<CatalogEntry>& catalog_entries();

const CatalogEntry* find_catalog_entry(const std::string& id);

std::vector<std::string> catalog_processes();
std::vector<std::string> catalog_regions();

// Canonicalises a free-text process name, returning an empty string if unknown.
std::string canonical_process(const std::string& name);

// Canonicalises a free-text region name, returning an empty string if unknown.
std::string canonical_region(const std::string& name);

}  // namespace climate
