#include "copula.hpp"

#include <cmath>
#include <stdexcept>

namespace climate {
namespace {

double polynomial(double gamma0, double gamma1, double gamma2, double covariate) {
    return gamma0 + gamma1 * covariate + gamma2 * covariate * covariate;
}

double polynomial_derivative(double gamma1, double gamma2, double covariate) {
    return gamma1 + 2.0 * gamma2 * covariate;
}

}  // namespace

double copula_parameter(CopulaFamily family, double gamma0, double gamma1, double gamma2,
                        double covariate) {
    const double exponent = std::exp(polynomial(gamma0, gamma1, gamma2, covariate));
    switch (family) {
        case CopulaFamily::Gumbel:
        case CopulaFamily::Joe:
            return 1.0 + exponent;
        case CopulaFamily::Clayton:
            return exponent - 1.0;
        case CopulaFamily::HuelserReiss:
        case CopulaFamily::Galambos:
            return exponent;
    }
    throw std::runtime_error("copula_parameter: unknown family");
}

double copula_parameter_derivative(CopulaFamily family, double gamma0, double gamma1,
                                   double gamma2, double covariate) {
    const double exponent = std::exp(polynomial(gamma0, gamma1, gamma2, covariate));
    const double derivative = polynomial_derivative(gamma1, gamma2, covariate);
    switch (family) {
        case CopulaFamily::Gumbel:
        case CopulaFamily::Joe:
        case CopulaFamily::Clayton:
        case CopulaFamily::HuelserReiss:
        case CopulaFamily::Galambos:
            return exponent * derivative;
    }
    throw std::runtime_error("copula_parameter_derivative: unknown family");
}

}  // namespace climate
