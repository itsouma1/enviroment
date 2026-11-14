#include "hmm.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace climate {
namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double emission(double observation, double mean, double variance) {
    require(variance > 0.0, "hmm: emission variance must be positive");
    const double difference = observation - mean;
    return std::exp(-0.5 * difference * difference / variance) /
           std::sqrt(2.0 * kPi * variance);
}

void validate(const GaussianHmm& model) {
    const std::size_t states = model.initial.size();
    require(states > 0, "hmm: no states");
    require(model.transition.size() == states && model.means.size() == states &&
                model.variances.size() == states,
            "hmm: state size mismatch");
    for (const auto& row : model.transition) require(row.size() == states, "hmm: ragged A");
    for (double variance : model.variances) require(variance > 0.0, "hmm: variance");
}

}  // namespace

HmmPosterior hmm_forward_backward(const GaussianHmm& model,
                                  const std::vector<double>& observations) {
    validate(model);
    require(!observations.empty(), "hmm_forward_backward: empty observations");
    const std::size_t states = model.initial.size();
    const std::size_t length = observations.size();

    Matrix alpha(length, std::vector<double>(states, 0.0));
    std::vector<double> scales(length, 1.0);

    double normaliser = 0.0;
    for (std::size_t state = 0; state < states; ++state) {
        alpha[0][state] = model.initial[state] * emission(observations[0], model.means[state],
                                                         model.variances[state]);
        normaliser += alpha[0][state];
    }
    require(normaliser > 0.0, "hmm_forward_backward: zero likelihood at t=0");
    for (double& value : alpha[0]) value /= normaliser;
    scales[0] = normaliser;

    for (std::size_t time = 1; time < length; ++time) {
        normaliser = 0.0;
        for (std::size_t state = 0; state < states; ++state) {
            double previous = 0.0;
            for (std::size_t from = 0; from < states; ++from)
                previous += alpha[time - 1][from] * model.transition[from][state];
            alpha[time][state] = previous * emission(observations[time], model.means[state],
                                                    model.variances[state]);
            normaliser += alpha[time][state];
        }
        require(normaliser > 0.0, "hmm_forward_backward: zero likelihood");
        for (double& value : alpha[time]) value /= normaliser;
        scales[time] = normaliser;
    }

    Matrix beta(length, std::vector<double>(states, 1.0));
    for (std::size_t time = length - 1; time-- > 0;) {
        for (std::size_t state = 0; state < states; ++state) {
            double sum = 0.0;
            for (std::size_t next = 0; next < states; ++next)
                sum += model.transition[state][next] *
                       emission(observations[time + 1], model.means[next],
                                model.variances[next]) *
                       beta[time + 1][next];
            beta[time][state] = sum / scales[time + 1];
        }
    }

    HmmPosterior posterior;
    posterior.log_likelihood = 0.0;
    for (double scale : scales) posterior.log_likelihood += std::log(scale);

    posterior.gamma.assign(length, std::vector<double>(states, 0.0));
    for (std::size_t time = 0; time < length; ++time) {
        double total = 0.0;
        for (std::size_t state = 0; state < states; ++state) {
            posterior.gamma[time][state] = alpha[time][state] * beta[time][state];
            total += posterior.gamma[time][state];
        }
        require(total > 0.0, "hmm_forward_backward: zero posterior");
        for (double& value : posterior.gamma[time]) value /= total;
    }

    posterior.xi.assign(length - 1, Matrix(states, std::vector<double>(states, 0.0)));
    for (std::size_t time = 0; time + 1 < length; ++time) {
        double total = 0.0;
        for (std::size_t from = 0; from < states; ++from) {
            for (std::size_t to = 0; to < states; ++to) {
                const double value =
                    alpha[time][from] * model.transition[from][to] *
                    emission(observations[time + 1], model.means[to], model.variances[to]) *
                    beta[time + 1][to] / scales[time + 1];
                posterior.xi[time][from][to] = value;
                total += value;
            }
        }
        require(total > 0.0, "hmm_forward_backward: zero transition posterior");
        for (auto& row : posterior.xi[time])
            for (double& value : row) value /= total;
    }
    return posterior;
}

GaussianHmm baum_welch(const GaussianHmm& model, const std::vector<double>& observations,
                       int max_iterations, double tolerance) {
    require(max_iterations > 0, "baum_welch: iterations must be positive");
    GaussianHmm current = model;
    const std::size_t states = current.initial.size();
    const std::size_t length = observations.size();

    double previous_likelihood = -std::numeric_limits<double>::infinity();
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        const HmmPosterior posterior = hmm_forward_backward(current, observations);
        if (std::fabs(posterior.log_likelihood - previous_likelihood) < tolerance) break;
        previous_likelihood = posterior.log_likelihood;

        for (std::size_t state = 0; state < states; ++state)
            current.initial[state] = posterior.gamma[0][state];

        for (std::size_t from = 0; from < states; ++from) {
            double denominator = 0.0;
            for (std::size_t time = 0; time + 1 < length; ++time)
                denominator += posterior.gamma[time][from];
            for (std::size_t to = 0; to < states; ++to) {
                double numerator = 0.0;
                for (std::size_t time = 0; time + 1 < length; ++time)
                    numerator += posterior.xi[time][from][to];
                current.transition[from][to] =
                    denominator > 0.0 ? numerator / denominator : 0.0;
            }
        }

        for (std::size_t state = 0; state < states; ++state) {
            double weight = 0.0, weighted_sum = 0.0;
            for (std::size_t time = 0; time < length; ++time) {
                weight += posterior.gamma[time][state];
                weighted_sum += posterior.gamma[time][state] * observations[time];
            }
            if (weight <= 0.0) continue;
            const double mean = weighted_sum / weight;
            double variance = 0.0;
            for (std::size_t time = 0; time < length; ++time) {
                const double difference = observations[time] - mean;
                variance += posterior.gamma[time][state] * difference * difference;
            }
            current.means[state] = mean;
            current.variances[state] = std::max(variance / weight, 1e-12);
        }
    }
    return current;
}

}  // namespace climate
