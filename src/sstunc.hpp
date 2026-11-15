#ifndef CLIMATE_SSTUNC_HPP
#define CLIMATE_SSTUNC_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

struct SstUncertainty {
    double uncorrelated = 0.0;
    double simply_correlated = 0.0;
    double ensemble = 0.0;
    double total = 0.0;
};

double uncorrelated_variance(const std::vector<double>& weights,
                             const std::vector<double>& sigmas);

double annual_average_variance(const std::vector<double>& monthly_sigmas);

double covariance_variance(const std::vector<double>& weights, const Matrix& covariance);

double hemispheric_difference_variance(const std::vector<double>& north_weights,
                                       const std::vector<double>& south_weights,
                                       const Matrix& covariance);

double constant_covariance_variance(double sigma, double correlation,
                                    const std::vector<double>& weights);

double effective_sample_size(const std::vector<double>& weights, double correlation);

Matrix sample_correlated_errors(const Matrix& covariance, std::size_t samples, unsigned seed);

double ensemble_standard_deviation(const std::vector<double>& values);

double combine_uncertainty(double uncorrelated, double simply_correlated, double ensemble);

double uncorrelated_weight_sensitivity(double weight, double sigma);

double covariance_entry_sensitivity(double first_weight, double second_weight,
                                    bool diagonal);

double correlation_sensitivity(double sigma, const std::vector<double>& weights);

double total_component_sensitivity(double component, double total);

double covariance_condition_number(const Matrix& covariance);

}  // namespace climate

#endif
