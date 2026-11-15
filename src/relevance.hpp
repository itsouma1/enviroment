#ifndef CLIMATE_RELEVANCE_HPP
#define CLIMATE_RELEVANCE_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

double variance_explained(const std::vector<double>& response, const Matrix& predictors);

double gaussian_mutual_information(double variance_explained);

double relevance_score(double variance_explained, double mutual_information,
                       double penalty, double lambda_g, double lambda_i,
                       double lambda_penalty);

double effective_relevance_ratio(double lambda_g, double lambda_i, double lambda_penalty);

double marginal_relevance_gain(double variance_explained, double mutual_information,
                              double delta_variance, double delta_information,
                              double cost, double lambda_g, double lambda_i,
                              double lambda_penalty);

double aic_penalty(std::size_t variables);

double bic_penalty(std::size_t variables, std::size_t samples);

std::vector<double> variance_inflation_factors(const Matrix& predictors);

std::vector<std::size_t> greedy_relevance_selection(const std::vector<double>& response,
                                                    const Matrix& predictors, double kappa,
                                                    std::size_t penalty_form);

}  // namespace climate

#endif
