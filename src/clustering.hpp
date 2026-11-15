#ifndef CLIMATE_CLUSTERING_HPP
#define CLIMATE_CLUSTERING_HPP

#include <cstddef>
#include <vector>

#include "ml.hpp"

namespace climate {

// ---- Fuzzy C-means -------------------------------------------------------

struct FuzzyCMeansResult {
    Matrix membership;  // N x C
    Matrix centroids;   // C x D
    double objective = 0.0;
    int iterations = 0;
};

double fuzzy_objective(const Matrix& points, const Matrix& membership,
                       const Matrix& centroids, double fuzzifier);

FuzzyCMeansResult fuzzy_c_means(const Matrix& points, std::size_t clusters,
                                double fuzzifier, double tolerance, int max_iterations);

std::vector<double> fuzzy_membership_sensitivity(const std::vector<double>& point,
                                                 const Matrix& centroids,
                                                 const std::vector<double>& membership);

double fuzzy_centroid_trend(const std::vector<double>& trends,
                            const std::vector<double>& memberships, double fuzzifier);

// ---- Attribute-weighted DBSCAN ------------------------------------------

double attribute_distance(const std::vector<double>& first,
                          const std::vector<double>& second, double attribute_weight);

// Labels are 0..K-1 for clusters and -1 for noise.
std::vector<int> dbscan(const Matrix& points, double epsilon, std::size_t min_points);

double dbscan_critical_weight(double epsilon, double spatial_distance, double attribute_gap);

// Spatial coordinates in all but the last column, weighted attribute in the last.
std::vector<int> dbscan_attribute(const Matrix& points, double attribute_weight,
                                  double epsilon, std::size_t min_points);

double cluster_mean_slope(const std::vector<double>& slopes, const std::vector<int>& labels,
                          int cluster);

// ---- REOF partition sensitivities ---------------------------------------

double varimax_criterion(const Matrix& loadings);

std::vector<int> reof_partition(const Matrix& loadings, double threshold);

// Sensitivity of each RPC to the data: T x sites.
Matrix rpc_trend_sensitivity(const std::vector<double>& times, const Matrix& loadings,
                             const Matrix& rotation, std::size_t mode);

double trend_sign_flip_probability(double trend, double standard_error, double alpha);

// ---- Hierarchical cluster assignment ------------------------------------

std::vector<double> categorical_assignment_probabilities(
    const std::vector<double>& log_likelihoods);

}  // namespace climate

#endif
