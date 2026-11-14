#pragma once

#include <vector>

namespace climate {

// Solves A x = b for a square system using Gauss-Jordan elimination with partial
// pivoting. Returns false if the matrix is singular to working precision.
bool solve_linear_system(std::vector<std::vector<double>> a, std::vector<double> b,
                         std::vector<double>& x);

// Least squares solution of the overdetermined system A x ~= b using Householder
// QR factorisation. A is given as m rows of n columns with m >= n. Returns false
// if the design matrix is rank deficient.
bool least_squares(const std::vector<std::vector<double>>& a, const std::vector<double>& b,
                   std::vector<double>& x);

// Dot product helper (kept here so callers do not re-implement it).
double dot_product(const std::vector<double>& a, const std::vector<double>& b);

}  // namespace climate
