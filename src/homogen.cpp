#include "homogen.hpp"

#include <cmath>
#include <map>
#include <stdexcept>

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double mean_of(const std::vector<double>& values, std::size_t begin, std::size_t end) {
    require(begin < end, "mean_of: empty segment");
    double total = 0.0;
    for (std::size_t i = begin; i < end; ++i) total += values[i];
    return total / static_cast<double>(end - begin);
}

double normal_cdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

}  // namespace

double absolute_mean_bias(const std::vector<double>& homogenised,
                          const std::vector<double>& truth) {
    require(!homogenised.empty() && homogenised.size() == truth.size(),
            "absolute_mean_bias: size mismatch");
    double total = 0.0;
    for (std::size_t i = 0; i < homogenised.size(); ++i)
        total += std::fabs(homogenised[i] - truth[i]);
    return total / static_cast<double>(homogenised.size());
}

double factor_of_exceedance(const std::vector<double>& homogenised,
                            const std::vector<double>& corrupted) {
    require(!homogenised.empty() && homogenised.size() == corrupted.size(),
            "factor_of_exceedance: size mismatch");
    int above = 0, below = 0;
    for (std::size_t i = 0; i < homogenised.size(); ++i) {
        if (homogenised[i] > corrupted[i]) ++above;
        if (homogenised[i] < corrupted[i]) ++below;
    }
    return 100.0 * static_cast<double>(above - below) /
           static_cast<double>(homogenised.size());
}

double relative_rmse(const std::vector<double>& homogenised,
                     const std::vector<double>& raw, const std::vector<double>& truth) {
    require(!homogenised.empty() && homogenised.size() == raw.size() &&
                homogenised.size() == truth.size(),
            "relative_rmse: size mismatch");
    double homogenised_error = 0.0, raw_error = 0.0;
    for (std::size_t i = 0; i < truth.size(); ++i) {
        homogenised_error += (homogenised[i] - truth[i]) * (homogenised[i] - truth[i]);
        raw_error += (raw[i] - truth[i]) * (raw[i] - truth[i]);
    }
    require(raw_error > 0.0, "relative_rmse: zero raw error");
    return std::sqrt(homogenised_error / raw_error);
}

double acmant_objective(const std::vector<double>& first, const std::vector<double>& second,
                        const std::vector<int>& break_indices, double weight) {
    require(first.size() == second.size(),
            "acmant_objective: size mismatch");
    require(!first.empty(), "acmant_objective: empty series");
    require(weight >= 0.0, "acmant_objective: negative weight");
    for (int index : break_indices)
        require(index >= 0 && static_cast<std::size_t>(index) < first.size(),
                "acmant_objective: break out of range");

    std::vector<std::size_t> boundaries{0};
    for (int index : break_indices)
        boundaries.push_back(static_cast<std::size_t>(index) + 1);
    boundaries.push_back(first.size());

    double objective = 0.0;
    for (std::size_t segment = 0; segment + 1 < boundaries.size(); ++segment) {
        const std::size_t begin = boundaries[segment];
        const std::size_t end = boundaries[segment + 1];
        if (begin >= end) continue;
        const double first_mean = mean_of(first, begin, end);
        const double second_mean = mean_of(second, begin, end);
        for (std::size_t i = begin; i < end; ++i) {
            const double first_deviation = first[i] - first_mean;
            const double second_deviation = second[i] - second_mean;
            objective += first_deviation * first_deviation +
                         weight * second_deviation * second_deviation;
        }
    }
    return objective;
}

void homogenisation_ensemble(const Matrix& members, std::vector<double>& mean,
                             std::vector<double>& variance) {
    require(!members.empty(), "homogenisation_ensemble: empty ensemble");
    const std::size_t length = members[0].size();
    for (const auto& member : members)
        require(member.size() == length, "homogenisation_ensemble: ragged member");
    mean.assign(length, 0.0);
    for (const auto& member : members)
        for (std::size_t t = 0; t < length; ++t) mean[t] += member[t];
    for (double& value : mean) value /= static_cast<double>(members.size());
    variance.assign(length, 0.0);
    if (members.size() < 2) return;
    for (const auto& member : members)
        for (std::size_t t = 0; t < length; ++t)
            variance[t] += (member[t] - mean[t]) * (member[t] - mean[t]);
    for (double& value : variance) value /= static_cast<double>(members.size() - 1);
}

double residual_cusum(const std::vector<double>& residuals, std::size_t candidate) {
    require(!residuals.empty() && candidate <= residuals.size(),
            "residual_cusum: candidate out of range");
    double partial = 0.0, total = 0.0;
    for (std::size_t i = 0; i < residuals.size(); ++i) {
        total += residuals[i];
        if (i < candidate) partial += residuals[i];
    }
    return partial - static_cast<double>(candidate) / static_cast<double>(residuals.size()) *
                         total;
}

double student_t_log_likelihood(const std::vector<double>& residuals, double degrees_of_freedom,
                                double scale) {
    require(degrees_of_freedom > 0.0 && scale > 0.0, "student_t_log_likelihood: parameters");
    const double constant = std::lgamma(0.5 * (degrees_of_freedom + 1.0)) -
                            std::lgamma(0.5 * degrees_of_freedom) -
                            0.5 * std::log(degrees_of_freedom * kPi * scale * scale);
    double total = 0.0;
    for (double residual : residuals) {
        const double standardised = residual / scale;
        total += constant - 0.5 * (degrees_of_freedom + 1.0) *
                                std::log(1.0 + standardised * standardised /
                                                   degrees_of_freedom);
    }
    return total;
}

double singh_maddala_density(double value, double a, double b, double q) {
    require(value > 0.0 && a > 0.0 && b > 0.0 && q > 0.0, "singh_maddala_density: inputs");
    const double ratio = value / b;
    return a * q / b * std::pow(ratio, a - 1.0) *
           std::pow(1.0 + std::pow(ratio, a), -(q + 1.0));
}

double sparse_jump_objective(const Matrix& observations, const std::vector<int>& states,
                             const Matrix& state_parameters, double jump_penalty) {
    require(!observations.empty() && observations.size() == states.size() &&
                states.size() == state_parameters.size(),
            "sparse_jump_objective: size mismatch");
    require(!observations[0].empty(), "sparse_jump_objective: empty observation");
    double objective = 0.0;
    for (std::size_t t = 0; t < observations.size(); ++t) {
        require(observations[t].size() == observations[0].size() &&
                    state_parameters[t].size() == observations[0].size(),
                "sparse_jump_objective: ragged");
        for (std::size_t j = 0; j < observations[t].size(); ++j) {
            const double difference = observations[t][j] - state_parameters[t][j];
            objective += difference * difference;
        }
        if (t > 0 && states[t] != states[t - 1]) objective += jump_penalty;
    }
    return objective;
}

double cautious_correction(double reported, double estimated_bias) {
    return reported - estimated_bias;
}

double surface_urban_heat_island(double baseline, double nonurban_minimum) {
    return baseline - nonurban_minimum;
}

double excursion_probability(double trend, double domain_trend, double standard_error) {
    require(standard_error > 0.0, "excursion_probability: standard error");
    return normal_cdf((trend - domain_trend) / standard_error);
}

double variational_lower_bound(double expected_log_likelihood, double kl_divergence) {
    return expected_log_likelihood - kl_divergence;
}

}  // namespace climate
