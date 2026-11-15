#include "regional.hpp"

#include "stats.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double ols_slope(const std::vector<double>& values) {
    const double n = static_cast<double>(values.size());
    double mean_t = 0.5 * (n - 1.0), mean_y = 0.0;
    for (double value : values) mean_y += value;
    mean_y /= n;
    double covariance = 0.0, variance = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double dt = static_cast<double>(i) - mean_t;
        covariance += dt * (values[i] - mean_y);
        variance += dt * dt;
    }
    require(variance > 0.0, "ols_slope: zero variance");
    return covariance / variance;
}

double column_mean(const Matrix& data, std::size_t column) {
    double sum = 0.0;
    for (const auto& row : data) sum += row[column];
    return sum / static_cast<double>(data.size());
}

double column_sd(const Matrix& data, std::size_t column) {
    const double mean = column_mean(data, column);
    double sum = 0.0;
    for (const auto& row : data) {
        const double difference = row[column] - mean;
        sum += difference * difference;
    }
    return std::sqrt(sum / static_cast<double>(data.size() - 1));
}

// Eigen-decomposition of a small symmetric matrix by cyclic Jacobi rotations.
EigenSystem jacobi(const Matrix& input, int max_iterations) {
    const std::size_t n = input.size();
    Matrix a = input;
    Matrix v(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) v[i][i] = 1.0;
    for (int sweep = 0; sweep < max_iterations; ++sweep) {
        double off_diagonal = 0.0;
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = p + 1; q < n; ++q) off_diagonal += a[p][q] * a[p][q];
        if (off_diagonal < 1e-24) break;
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (std::fabs(a[p][q]) < 1e-18) continue;
                const double theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (std::size_t k = 0; k < n; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    EigenSystem result;
    result.values.resize(n);
    result.vectors = v;
    for (std::size_t i = 0; i < n; ++i) result.values[i] = a[i][i];
    // Sort eigenvectors by descending eigenvalue.
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](std::size_t i, std::size_t j) { return result.values[i] > result.values[j]; });
    std::vector<double> sorted_values(n);
    Matrix sorted_vectors(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        sorted_values[i] = result.values[order[i]];
        for (std::size_t row = 0; row < n; ++row)
            sorted_vectors[row][i] = result.vectors[row][order[i]];
    }
    result.values = sorted_values;
    result.vectors = sorted_vectors;
    return result;
}

// SVD of a small square matrix via the eigen-decomposition of B^T B.
bool small_svd(const Matrix& b, Matrix& u, std::vector<double>& singular_values,
               Matrix& w) {
    const std::size_t n = b.size();
    Matrix bt_b(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t k = 0; k < n; ++k) bt_b[i][j] += b[k][i] * b[k][j];
    const EigenSystem system = jacobi(bt_b, 100);
    w = system.vectors;
    singular_values.assign(n, 0.0);
    u.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t j = 0; j < n; ++j) {
        singular_values[j] = std::sqrt(std::max(0.0, system.values[j]));
        if (singular_values[j] < 1e-12) return false;
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k) sum += b[i][k] * w[k][j];
            u[i][j] = sum / singular_values[j];
        }
    return true;
}

}  // namespace

double sign_difference_statistic(std::size_t positive, std::size_t negative) {
    const std::size_t total = positive + negative;
    require(total > 0, "sign_difference_statistic: no trends");
    return (static_cast<double>(positive) - static_cast<double>(negative)) /
           static_cast<double>(total);
}

double sign_test_percentage(std::size_t positive, std::size_t total) {
    require(total > 0, "sign_test_percentage: no trends");
    return (static_cast<double>(positive) / static_cast<double>(total) - 1.0) * 100.0;
}

double sign_test_binomial_p_value(std::size_t observed, std::size_t trials) {
    require(observed <= trials && trials > 0, "sign_test_binomial_p_value: counts");
    const double n = static_cast<double>(trials);
    double tail = 0.0;
    for (std::size_t i = observed; i <= trials; ++i) {
        const double log_choose =
            std::lgamma(n + 1.0) - std::lgamma(static_cast<double>(i) + 1.0) -
            std::lgamma(n - static_cast<double>(i) + 1.0);
        tail += std::exp(log_choose - n * std::log(2.0));
    }
    return std::min(1.0, tail);
}

double sign_test_empirical_p_value(std::size_t observed,
                                   const std::vector<std::size_t>& ensemble_counts) {
    require(!ensemble_counts.empty(), "sign_test_empirical_p_value: empty ensemble");
    std::size_t exceedances = 0;
    for (std::size_t count : ensemble_counts)
        if (count >= observed) ++exceedances;
    return static_cast<double>(exceedances) / static_cast<double>(ensemble_counts.size());
}

double effective_sample_size(double samples, double spatial_autocorrelation) {
    require(samples > 0.0, "effective_sample_size: samples");
    require(spatial_autocorrelation > -1.0 && spatial_autocorrelation < 1.0,
            "effective_sample_size: autocorrelation");
    return samples * (1.0 - spatial_autocorrelation) / (1.0 + spatial_autocorrelation);
}

double sign_test_critical_value(double samples, double spatial_autocorrelation,
                                double alpha) {
    require(samples > 0.0, "sign_test_critical_value: samples");
    const double effective = effective_sample_size(samples, spatial_autocorrelation);
    return normal_quantile(1.0 - alpha) / std::sqrt(effective);
}

double sign_test_critical_value_derivative(double samples, double spatial_autocorrelation,
                                           double alpha) {
    require(samples > 0.0, "sign_test_critical_value_derivative: samples");
    require(spatial_autocorrelation >= 0.0 && spatial_autocorrelation < 1.0,
            "sign_test_critical_value_derivative: autocorrelation");
    const double z = normal_quantile(1.0 - alpha);
    return z / std::sqrt(samples) /
           (std::pow(1.0 - spatial_autocorrelation, 1.5) *
            std::sqrt(1.0 + spatial_autocorrelation));
}

double zonal_mean_at_latitude(const std::vector<double>& trends,
                              const std::vector<int>& latitude_index, int latitude) {
    require(trends.size() == latitude_index.size(), "zonal_mean_at_latitude: size");
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < trends.size(); ++i)
        if (latitude_index[i] == latitude) {
            sum += trends[i];
            ++count;
        }
    require(count > 0, "zonal_mean_at_latitude: empty latitude band");
    return sum / static_cast<double>(count);
}

std::vector<double> zonally_asymmetric_trend(const std::vector<double>& trends,
                                             const std::vector<int>& latitude_index) {
    require(trends.size() == latitude_index.size(), "zonally_asymmetric_trend: size");
    std::vector<double> asymmetric(trends.size());
    for (std::size_t i = 0; i < trends.size(); ++i)
        asymmetric[i] = trends[i] - zonal_mean_at_latitude(trends, latitude_index,
                                                           latitude_index[i]);
    return asymmetric;
}

double regional_mean(const std::vector<double>& field, const std::vector<int>& region_mask) {
    require(field.size() == region_mask.size(), "regional_mean: size");
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < field.size(); ++i)
        if (region_mask[i] != 0) {
            sum += field[i];
            ++count;
        }
    require(count > 0, "regional_mean: empty region");
    return sum / static_cast<double>(count);
}

double sign_agreement(const std::vector<double>& observed,
                      const std::vector<double>& estimated,
                      const std::vector<int>& land_mask) {
    require(observed.size() == estimated.size() && observed.size() == land_mask.size(),
            "sign_agreement: size");
    std::size_t agree = 0, count = 0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        if (land_mask[i] == 0) continue;
        const int first = (observed[i] > 0.0) - (observed[i] < 0.0);
        const int second = (estimated[i] > 0.0) - (estimated[i] < 0.0);
        if (first == second) ++agree;
        ++count;
    }
    require(count > 0, "sign_agreement: empty mask");
    return static_cast<double>(agree) / static_cast<double>(count);
}

double gev_location_score(double value, double location, double scale, double shape) {
    require(scale > 0.0, "gev_location_score: scale");
    const double z = (value - location) / scale;
    const double support = 1.0 + shape * z;
    require(support > 0.0, "gev_location_score: value outside support");
    if (std::fabs(shape) < 1e-12) {
        const double exponential = std::exp(-z);
        return (1.0 - exponential) / scale;
    }
    return ((1.0 + shape) / support - std::pow(support, -1.0 / shape - 1.0)) / scale;
}

double gev_trend_score(double value, double location, double scale, double shape,
                       double time) {
    return time * gev_location_score(value, location, scale, shape);
}

double shrinkage_log_prior(double site_trend, double cluster_mean,
                           double between_site_variance) {
    require(between_site_variance > 0.0, "shrinkage_log_prior: variance");
    const double difference = site_trend - cluster_mean;
    return -0.5 * difference * difference / between_site_variance;
}

double shrinkage_derivative(double site_trend, double cluster_mean,
                            double between_site_variance) {
    require(between_site_variance > 0.0, "shrinkage_derivative: variance");
    return -(site_trend - cluster_mean) / between_site_variance;
}

double pooled_posterior_variance(double at_site_variance, double between_site_variance) {
    require(at_site_variance > 0.0 && between_site_variance > 0.0,
            "pooled_posterior_variance: variances");
    return 1.0 / (1.0 / at_site_variance + 1.0 / between_site_variance);
}

double pooling_standard_error_ratio(double between_site_variance, double data_precision) {
    require(between_site_variance >= 0.0 && data_precision > 0.0,
            "pooling_standard_error_ratio: inputs");
    return 1.0 / std::sqrt(1.0 + between_site_variance * data_precision);
}

double pooled_cluster_mean_sensitivity(const std::vector<double>& precisions,
                                       std::size_t site) {
    require(site < precisions.size(), "pooled_cluster_mean_sensitivity: site");
    double total = 0.0;
    for (double precision : precisions) total += precision;
    require(total > 0.0, "pooled_cluster_mean_sensitivity: zero precision");
    return precisions[site] / total;
}

double pooled_cluster_mean(const std::vector<double>& site_trends,
                           const std::vector<double>& precisions) {
    require(!site_trends.empty() && site_trends.size() == precisions.size(),
            "pooled_cluster_mean: size");
    double weighted = 0.0, total = 0.0;
    for (std::size_t i = 0; i < site_trends.size(); ++i) {
        require(precisions[i] > 0.0, "pooled_cluster_mean: precision");
        weighted += precisions[i] * site_trends[i];
        total += precisions[i];
    }
    return weighted / total;
}

EigenSystem symmetric_eigen(const Matrix& matrix) {
    require(!matrix.empty() && matrix.size() == matrix[0].size(),
            "symmetric_eigen: square matrix required");
    for (std::size_t i = 0; i < matrix.size(); ++i)
        for (std::size_t j = 0; j < matrix.size(); ++j)
            require(std::fabs(matrix[i][j] - matrix[j][i]) < 1e-10,
                    "symmetric_eigen: matrix not symmetric");
    return jacobi(matrix, 100);
}

double north_test_standard_error(double eigenvalue, double effective_dof) {
    require(eigenvalue >= 0.0 && effective_dof > 0.0, "north_test_standard_error: inputs");
    return eigenvalue * std::sqrt(2.0 / effective_dof);
}

Matrix temporal_correlation_matrix(const Matrix& data) {
    require(data.size() > 1 && !data[0].empty(), "temporal_correlation_matrix: data");
    const std::size_t sites = data[0].size();
    Matrix correlation(sites, std::vector<double>(sites, 0.0));
    std::vector<double> means(sites), sds(sites);
    for (std::size_t s = 0; s < sites; ++s) {
        means[s] = column_mean(data, s);
        sds[s] = column_sd(data, s);
        require(sds[s] > 0.0, "temporal_correlation_matrix: zero variance");
    }
    for (std::size_t a = 0; a < sites; ++a)
        for (std::size_t b = 0; b < sites; ++b) {
            double sum = 0.0;
            for (const auto& row : data)
                sum += (row[a] - means[a]) * (row[b] - means[b]);
            correlation[a][b] =
                sum / (static_cast<double>(data.size() - 1) * sds[a] * sds[b]);
        }
    return correlation;
}

Matrix varimax_rotation(const Matrix& loadings, int max_iterations, double tolerance,
                        Matrix& rotation) {
    require(!loadings.empty() && !loadings[0].empty(), "varimax_rotation: loadings");
    const std::size_t sites = loadings.size();
    const std::size_t modes = loadings[0].size();
    rotation.assign(modes, std::vector<double>(modes, 0.0));
    for (std::size_t i = 0; i < modes; ++i) rotation[i][i] = 1.0;
    if (modes == 1) return loadings;
    Matrix current = loadings;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::vector<double> sums(modes, 0.0);
        for (std::size_t s = 0; s < sites; ++s)
            for (std::size_t k = 0; k < modes; ++k) sums[k] += current[s][k] * current[s][k];
        Matrix target(sites, std::vector<double>(modes, 0.0));
        for (std::size_t s = 0; s < sites; ++s)
            for (std::size_t k = 0; k < modes; ++k)
                target[s][k] = current[s][k] * current[s][k] * current[s][k] -
                               current[s][k] * sums[k] / static_cast<double>(sites);
        Matrix b(modes, std::vector<double>(modes, 0.0));
        for (std::size_t i = 0; i < modes; ++i)
            for (std::size_t j = 0; j < modes; ++j)
                for (std::size_t s = 0; s < sites; ++s)
                    b[i][j] += current[s][i] * target[s][j];
        Matrix u, w;
        std::vector<double> singular;
        if (!small_svd(b, u, singular, w)) break;
        Matrix candidate(modes, std::vector<double>(modes, 0.0));
        for (std::size_t i = 0; i < modes; ++i)
            for (std::size_t j = 0; j < modes; ++j)
                for (std::size_t k = 0; k < modes; ++k)
                    candidate[i][j] += u[i][k] * w[j][k];
        double change = 0.0;
        for (std::size_t i = 0; i < modes; ++i)
            for (std::size_t j = 0; j < modes; ++j)
                change = std::max(change, std::fabs(candidate[i][j] - rotation[i][j]));
        rotation = candidate;
        current.assign(sites, std::vector<double>(modes, 0.0));
        for (std::size_t s = 0; s < sites; ++s)
            for (std::size_t k = 0; k < modes; ++k)
                for (std::size_t j = 0; j < modes; ++j)
                    current[s][k] += loadings[s][j] * rotation[j][k];
        if (change < tolerance) break;
    }
    return current;
}

Matrix rotate_loadings(const Matrix& loadings, const Matrix& rotation) {
    require(!loadings.empty() && rotation.size() == loadings[0].size(),
            "rotate_loadings: size");
    const std::size_t sites = loadings.size();
    const std::size_t modes = rotation.size();
    Matrix rotated(sites, std::vector<double>(modes, 0.0));
    for (std::size_t s = 0; s < sites; ++s)
        for (std::size_t k = 0; k < modes; ++k)
            for (std::size_t j = 0; j < modes; ++j)
                rotated[s][k] += loadings[s][j] * rotation[j][k];
    return rotated;
}

std::vector<double> reof_rpc_trends(const Matrix& data, const Matrix& loadings,
                                    const Matrix& rotation) {
    require(!data.empty() && !data[0].empty() && data[0].size() == loadings.size(),
            "reof_rpc_trends: size");
    const std::size_t modes = rotation.size();
    std::vector<double> trends(modes, 0.0);
    for (std::size_t k = 0; k < modes; ++k) {
        std::vector<double> rpc(data.size(), 0.0);
        for (std::size_t t = 0; t < data.size(); ++t)
            for (std::size_t s = 0; s < loadings.size(); ++s)
                for (std::size_t j = 0; j < modes; ++j)
                    rpc[t] += data[t][s] * loadings[s][j] * rotation[j][k];
        trends[k] = ols_slope(rpc);
    }
    return trends;
}

double reof_back_projection(const std::vector<double>& trends, const Matrix& loadings,
                            const Matrix& rotation, std::size_t site) {
    require(site < loadings.size() && trends.size() == rotation.size(),
            "reof_back_projection: size");
    double value = 0.0;
    for (std::size_t j = 0; j < rotation.size(); ++j)
        for (std::size_t k = 0; k < rotation.size(); ++k)
            value += trends[k] * loadings[site][j] * rotation[j][k];
    return value;
}

double weatherhead_trend_standard_error(double noise_sd, double autocorrelation,
                                        double years) {
    require(noise_sd > 0.0 && years > 0.0, "weatherhead_trend_standard_error: inputs");
    require(std::fabs(autocorrelation) < 1.0,
            "weatherhead_trend_standard_error: autocorrelation");
    return noise_sd / std::sqrt(years) *
           std::sqrt((1.0 + autocorrelation) / (1.0 - autocorrelation));
}

double weatherhead_se_autocorrelation_derivative(double noise_sd, double autocorrelation,
                                                 double years) {
    const double se = weatherhead_trend_standard_error(noise_sd, autocorrelation, years);
    return se / (1.0 - autocorrelation * autocorrelation);
}

double weatherhead_se_noise_derivative(double noise_sd, double autocorrelation,
                                       double years) {
    return weatherhead_trend_standard_error(noise_sd, autocorrelation, years) / noise_sd;
}

double minimum_detectable_trend(double noise_sd, double autocorrelation, double years,
                                double alpha) {
    return normal_quantile(1.0 - 0.5 * alpha) *
           weatherhead_trend_standard_error(noise_sd, autocorrelation, years);
}

double detection_time(double noise_sd, double autocorrelation, double trend_magnitude,
                      double alpha) {
    require(noise_sd > 0.0 && trend_magnitude > 0.0,
            "detection_time: inputs");
    const double z = normal_quantile(1.0 - 0.5 * alpha);
    const double base = z * noise_sd / trend_magnitude *
                        std::sqrt((1.0 + autocorrelation) / (1.0 - autocorrelation));
    return std::pow(base, 2.0 / 3.0);
}

double detection_time_autocorrelation_derivative(double noise_sd, double autocorrelation,
                                                 double trend_magnitude, double alpha) {
    const double n = detection_time(noise_sd, autocorrelation, trend_magnitude, alpha);
    require(std::fabs(autocorrelation) < 1.0, "detection_time_derivative: autocorrelation");
    return 2.0 / 3.0 * n /
           ((1.0 - autocorrelation) * (1.0 + autocorrelation));
}

double detection_time_noise_derivative(double noise_sd, double autocorrelation,
                                       double trend_magnitude, double alpha) {
    const double n = detection_time(noise_sd, autocorrelation, trend_magnitude, alpha);
    return 2.0 / 3.0 * n / noise_sd;
}

double trend_pattern_correlation(const std::vector<double>& first,
                                 const std::vector<double>& second) {
    require(!first.empty() && first.size() == second.size(),
            "trend_pattern_correlation: size");
    const double n = static_cast<double>(first.size());
    double mean_first = 0.0, mean_second = 0.0;
    for (std::size_t i = 0; i < first.size(); ++i) {
        mean_first += first[i];
        mean_second += second[i];
    }
    mean_first /= n;
    mean_second /= n;
    double covariance = 0.0, variance_first = 0.0, variance_second = 0.0;
    for (std::size_t i = 0; i < first.size(); ++i) {
        const double a = first[i] - mean_first, b = second[i] - mean_second;
        covariance += a * b;
        variance_first += a * a;
        variance_second += b * b;
    }
    require(variance_first > 0.0 && variance_second > 0.0,
            "trend_pattern_correlation: zero variance");
    return covariance / std::sqrt(variance_first * variance_second);
}

double effective_spatial_dof(double points, double spatial_autocorrelation) {
    return effective_sample_size(points, spatial_autocorrelation);
}

double correlation_t_statistic(double correlation, double effective_dof) {
    require(std::fabs(correlation) < 1.0 && effective_dof > 2.0,
            "correlation_t_statistic: inputs");
    return correlation * std::sqrt((effective_dof - 2.0) / (1.0 - correlation * correlation));
}

double coherence_height(const std::vector<double>& pressures,
                        const std::vector<double>& correlations, double critical_value) {
    require(pressures.size() == correlations.size() && !pressures.empty(),
            "coherence_height: size");
    for (std::size_t i = 0; i < correlations.size(); ++i)
        if (std::fabs(correlations[i]) < critical_value) return pressures[i];
    return pressures.back();
}

double segmented_trend_difference(double second_period, double first_period) {
    return second_period - first_period;
}

std::size_t sign_test_max_count(std::size_t positive_sites, std::size_t negative_sites) {
    return std::max(positive_sites, negative_sites);
}

double sign_test_count_fraction(std::size_t positive_sites, std::size_t total_sites) {
    require(total_sites > 0, "sign_test_count_fraction: total");
    return static_cast<double>(positive_sites) / static_cast<double>(total_sites);
}

double empirical_critical_value(std::vector<double> null_counts, double alpha) {
    require(!null_counts.empty() && alpha > 0.0 && alpha < 1.0,
            "empirical_critical_value: inputs");
    std::sort(null_counts.begin(), null_counts.end());
    const std::size_t n = null_counts.size();
    std::size_t rank = static_cast<std::size_t>(std::ceil((1.0 - alpha) * static_cast<double>(n)));
    if (rank < 1) rank = 1;
    if (rank > n) rank = n;
    return null_counts[rank - 1];
}

double energy_budget_trend(double longwave_down, double shortwave_net, double latent_heat,
                           double sensible_heat, double storage) {
    return longwave_down + shortwave_net - latent_heat - sensible_heat + storage;
}

double budget_contribution(double term, double total) {
    require(std::fabs(total) > 0.0, "budget_contribution: zero total");
    return 100.0 * term / total;
}

double contribution_sensitivity(double term, double total, bool same_term) {
    require(std::fabs(total) > 0.0, "contribution_sensitivity: zero total");
    return 100.0 * ((same_term ? total : 0.0) - term) / (total * total);
}

double bayesian_posterior(double likelihood_inside, double likelihood_outside,
                          double prior_probability) {
    require(likelihood_inside >= 0.0 && likelihood_outside >= 0.0 &&
                prior_probability >= 0.0 && prior_probability <= 1.0,
            "bayesian_posterior: inputs");
    const double numerator = likelihood_inside * prior_probability;
    const double marginal = numerator + likelihood_outside * (1.0 - prior_probability);
    require(marginal > 0.0, "bayesian_posterior: zero marginal");
    return numerator / marginal;
}

double bayesian_log_odds(double likelihood_inside, double likelihood_outside,
                         double prior_probability) {
    require(likelihood_inside > 0.0 && likelihood_outside > 0.0 &&
                prior_probability > 0.0 && prior_probability < 1.0,
            "bayesian_log_odds: inputs");
    return std::log(prior_probability / (1.0 - prior_probability)) +
           std::log(likelihood_inside / likelihood_outside);
}

double posterior_prior_derivative(double likelihood_inside, double likelihood_outside,
                                  double prior_probability) {
    require(likelihood_inside >= 0.0 && likelihood_outside >= 0.0,
            "posterior_prior_derivative: likelihoods");
    const double marginal = likelihood_inside * prior_probability +
                            likelihood_outside * (1.0 - prior_probability);
    require(marginal > 0.0, "posterior_prior_derivative: zero marginal");
    const double posterior = likelihood_inside * prior_probability / marginal;
    return likelihood_inside / marginal -
           posterior * (likelihood_inside - likelihood_outside) / marginal;
}

double probability_weighted_centroid(const std::vector<double>& coordinates,
                                     const std::vector<double>& probabilities) {
    require(!coordinates.empty() && coordinates.size() == probabilities.size(),
            "probability_weighted_centroid: size");
    double weighted = 0.0, total = 0.0;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        weighted += coordinates[i] * probabilities[i];
        total += probabilities[i];
    }
    require(total > 0.0, "probability_weighted_centroid: zero probability");
    return weighted / total;
}

double centroid_probability_sensitivity(double coordinate, double centroid,
                                        double total_probability) {
    require(total_probability > 0.0, "centroid_probability_sensitivity: probability");
    return (coordinate - centroid) / total_probability;
}

double kolmogorov_smirnov_d_plus(const std::vector<double>& observed,
                                 const std::vector<double>& null_sample) {
    require(!observed.empty() && !null_sample.empty(), "kolmogorov_smirnov_d_plus: samples");
    std::vector<double> values = observed;
    values.insert(values.end(), null_sample.begin(), null_sample.end());
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    double maximum = 0.0;
    for (double value : values) {
        std::size_t observed_below = 0, null_below = 0;
        for (double sample : observed)
            if (sample <= value) ++observed_below;
        for (double sample : null_sample)
            if (sample <= value) ++null_below;
        const double null_cdf = static_cast<double>(null_below) /
                                static_cast<double>(null_sample.size());
        const double observed_cdf = static_cast<double>(observed_below) /
                                    static_cast<double>(observed.size());
        maximum = std::max(maximum, null_cdf - observed_cdf);
    }
    return maximum;
}

double optical_trend_precision(double prior_precision, const std::vector<double>& covariate,
                               double measurement_variance, double process_variance) {
    require(prior_precision >= 0.0 && measurement_variance > 0.0 &&
                process_variance >= 0.0,
            "optical_trend_precision: inputs");
    double precision = prior_precision;
    for (double value : covariate)
        precision += value * value / (measurement_variance + process_variance);
    return precision;
}

double elevation_trend_slope(const std::vector<double>& elevations,
                             const std::vector<double>& trends) {
    require(elevations.size() == trends.size() && elevations.size() > 1,
            "elevation_trend_slope: size");
    const double n = static_cast<double>(elevations.size());
    double mean_z = 0.0, mean_beta = 0.0;
    for (std::size_t i = 0; i < elevations.size(); ++i) {
        mean_z += elevations[i];
        mean_beta += trends[i];
    }
    mean_z /= n;
    mean_beta /= n;
    double covariance = 0.0, variance = 0.0;
    for (std::size_t i = 0; i < elevations.size(); ++i) {
        covariance += (elevations[i] - mean_z) * (trends[i] - mean_beta);
        variance += (elevations[i] - mean_z) * (elevations[i] - mean_z);
    }
    require(variance > 0.0, "elevation_trend_slope: zero variance");
    return covariance / variance;
}

double elevation_leverage(double elevation, const std::vector<double>& elevations) {
    require(elevations.size() > 1, "elevation_leverage: size");
    double mean = 0.0;
    for (double value : elevations) mean += value;
    mean /= static_cast<double>(elevations.size());
    double denominator = 0.0;
    for (double value : elevations) denominator += (value - mean) * (value - mean);
    require(denominator > 0.0, "elevation_leverage: zero variance");
    return (elevation - mean) / denominator;
}

double latitude_corrected_edw(double raw_edw, double latitude_sensitivity,
                              double latitude_elevation_slope) {
    return raw_edw - latitude_sensitivity * latitude_elevation_slope;
}

double lapse_rate_trend(double later_lapse_rate, double earlier_lapse_rate, double elapsed) {
    require(elapsed > 0.0, "lapse_rate_trend: elapsed");
    return (later_lapse_rate - earlier_lapse_rate) / elapsed;
}

}  // namespace climate
