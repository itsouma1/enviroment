#ifndef CLIMATE_COPULA_HPP
#define CLIMATE_COPULA_HPP

namespace climate {

enum class CopulaFamily { Gumbel, Joe, Clayton, HuelserReiss, Galambos };

// Time-varying dependence parameter gamma(nu) with the family transform that
// keeps it inside the valid range.
double copula_parameter(CopulaFamily family, double gamma0, double gamma1, double gamma2,
                        double covariate);

// d gamma / d nu.
double copula_parameter_derivative(CopulaFamily family, double gamma0, double gamma1,
                                   double gamma2, double covariate);

}  // namespace climate

#endif
