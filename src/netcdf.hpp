#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace climate {

// NetCDF type codes as stored in the classic (CDF-1 / CDF-2) header.
enum class NcType : int {
    Byte = 1,
    Char = 2,
    Short = 3,
    Int = 4,
    Float = 5,
    Double = 6,
};

struct NcDimension {
    std::string name;
    std::size_t length = 0;
    bool unlimited = false;
};

struct NcAttribute {
    std::string name;
    int type = 0;
    std::size_t count = 0;
    std::vector<double> numbers;
    std::string text;
};

struct NcVariable {
    std::string name;
    int type = 0;
    std::vector<int> dimension_ids;
    std::vector<std::size_t> shape;
    std::size_t element_count = 0;
    std::uint64_t data_offset = 0;
    std::uint32_t vsize = 0;
    bool record = false;
    std::vector<NcAttribute> attributes;
};

// Pure-C++ reader for the NetCDF-3 classic format (CDF-1 and CDF-2). NetCDF-4 /
// HDF5 files are detected and reported, and are converted by the Python
// data wrapper because HDF5 is a container format rather than a data layout.
class NetCDFFile {
public:
    static NetCDFFile open(const std::string& path);

    int version() const { return version_; }
    std::size_t record_count() const { return record_count_; }
    std::size_t record_stride() const { return record_stride_; }
    bool classic() const { return classic_; }

    const std::string& path() const { return path_; }
    const std::vector<NcDimension>& dimensions() const { return dimensions_; }
    const std::vector<NcVariable>& variables() const { return variables_; }
    const std::vector<NcAttribute>& global_attributes() const { return global_attributes_; }

    const NcVariable* find_variable(const std::string& name) const;

    // Returns the variable's values in C order (last dimension varies fastest),
    // converted to double without applying scale_factor or add_offset.
    std::vector<double> read_variable(const std::string& name) const;

    // Convenience: numeric attribute value, or fallback.
    double attribute_number(const NcVariable& variable, const std::string& name,
                            double fallback) const;

    std::string attribute_text(const NcVariable& variable, const std::string& name) const;

private:
    std::string path_;
    int version_ = 0;
    bool classic_ = true;
    std::size_t record_count_ = 0;
    std::size_t record_stride_ = 0;
    std::vector<NcDimension> dimensions_;
    std::vector<NcVariable> variables_;
    std::vector<NcAttribute> global_attributes_;
};

std::string netcdf_type_name(int type);

bool is_hdf5_file(const std::string& path);

// Converts time coordinate values with a "since" unit to fractional years.
double decimal_year_from_time(double value, const std::string& units);

struct NetCDFExportResult {
    bool ok = false;
    std::string error;
    std::string summary_json;
    std::string output_path;
};

// Exports a classic NetCDF file to the CSV schema used by this tool.
// mode "timeseries": time,value (spatially averaged with cosine-latitude weights)
// mode "grid": lat_min,lat_max,lon_min,lon_max,concentration (last time step)
NetCDFExportResult export_netcdf_to_csv(const std::string& path, const std::string& variable,
                                        const std::string& mode, const std::string& output);

}  // namespace climate
