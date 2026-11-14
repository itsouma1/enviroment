#include "linalg.hpp"

#include <cmath>
#include <cstddef>

namespace climate {

double dot_product(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = a.size() < b.size() ? a.size() : b.size();
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) sum += a[i] * b[i];
    return sum;
}

bool solve_linear_system(std::vector<std::vector<double>> a, std::vector<double> b,
                         std::vector<double>& x) {
    const std::size_t n = b.size();
    if (a.size() != n) return false;
    for (const auto& row : a)
        if (row.size() != n) return false;

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        double best = std::fabs(a[col][col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const double candidate = std::fabs(a[row][col]);
            if (candidate > best) {
                best = candidate;
                pivot = row;
            }
        }
        if (best < 1e-14) return false;
        if (pivot != col) {
            std::swap(a[pivot], a[col]);
            std::swap(b[pivot], b[col]);
        }
        const double diag = a[col][col];
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const double factor = a[row][col] / diag;
            if (factor == 0.0) continue;
            for (std::size_t k = col; k < n; ++k) a[row][k] -= factor * a[col][k];
            b[row] -= factor * b[col];
        }
    }

    x.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (std::fabs(a[i][i]) < 1e-14) return false;
        x[i] = b[i] / a[i][i];
    }
    return true;
}

bool least_squares(const std::vector<std::vector<double>>& a, const std::vector<double>& b,
                   std::vector<double>& x) {
    const std::size_t m = a.size();
    if (m == 0 || m != b.size()) return false;
    const std::size_t n = a[0].size();
    if (n == 0 || n > m) return false;
    for (const auto& row : a)
        if (row.size() != n) return false;

    std::vector<std::vector<double>> r(a);
    std::vector<double> rhs(b);

    for (std::size_t k = 0; k < n; ++k) {
        double norm = 0.0;
        for (std::size_t i = k; i < m; ++i) norm += r[i][k] * r[i][k];
        norm = std::sqrt(norm);
        if (norm < 1e-300) return false;

        const double alpha = r[k][k] > 0.0 ? -norm : norm;
        std::vector<double> v(m - k);
        for (std::size_t i = k; i < m; ++i) v[i - k] = r[i][k];
        v[0] -= alpha;

        double vnorm2 = 0.0;
        for (double value : v) vnorm2 += value * value;
        if (vnorm2 < 1e-300) continue;

        for (std::size_t j = k; j < n; ++j) {
            double projection = 0.0;
            for (std::size_t i = k; i < m; ++i) projection += v[i - k] * r[i][j];
            const double factor = 2.0 * projection / vnorm2;
            for (std::size_t i = k; i < m; ++i) r[i][j] -= factor * v[i - k];
        }
        double projection = 0.0;
        for (std::size_t i = k; i < m; ++i) projection += v[i - k] * rhs[i];
        const double factor = 2.0 * projection / vnorm2;
        for (std::size_t i = k; i < m; ++i) rhs[i] -= factor * v[i - k];
    }

    x.assign(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double sum = rhs[i];
        for (std::size_t j = i + 1; j < n; ++j) sum -= r[i][j] * x[j];
        if (std::fabs(r[i][i]) < 1e-14) return false;
        x[i] = sum / r[i][i];
    }
    return true;
}

}  // namespace climate
