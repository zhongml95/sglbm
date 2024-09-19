// hermite_basis.cpp
#include "hermite_basis.h"
#include <cmath>
#include <stdexcept>

namespace Polynomials {
namespace Hermite {

// Helper function to compute Hermite coefficients recursively
static std::vector<double> computeHermiteCoefficients(int n) {
    if (n == 0) {
        return {1.0};  // H_0(x) = 1
    } else if (n == 1) {
        return {0.0, 2.0};  // H_1(x) = 2x
    } else {
        // Recurrence relation: H_n(x) = 2x * H_{n-1}(x) - 2(n - 1) * H_{n-2}(x)
        std::vector<double> Hn_minus1 = computeHermiteCoefficients(n - 1);
        std::vector<double> Hn_minus2 = computeHermiteCoefficients(n - 2);
        std::vector<double> Hn(n + 1, 0.0);

        // Compute 2x * H_{n-1}
        for (int i = 0; i < Hn_minus1.size(); ++i) {
            Hn[i + 1] += 2.0 * Hn_minus1[i];
        }

        // Compute -2(n - 1) * H_{n-2}
        double c = -2.0 * (n - 1);
        for (int i = 0; i < Hn_minus2.size(); ++i) {
            Hn[i] += c * Hn_minus2[i];
        }

        // Normalize the coefficients by the largest power term
        for (int i = 0; i <= n; ++i) {
            Hn[i] /= Hn[n];
        }


        return Hn;
    }
}

// Compute coefficients for the Hermite polynomial of order n
std::vector<double> HermiteBasis::computeCoefficients(int n) const {
    return computeHermiteCoefficients(n);
}

// Construct the Jacobi matrix of size n for Hermite polynomials
std::vector<std::vector<double>> HermiteBasis::constructJacobiMatrix(int n) const {
    if (n <= 1) {
        throw std::invalid_argument("The order (n) must be greater than 1.");
    }

    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    // Fill the Jacobi matrix according to the recurrence relation for Hermite polynomials
    for (int i = 1; i < n; ++i) {
        J[i][i - 1] = std::sqrt(i / 2.0);
        J[i - 1][i] = J[i][i - 1];
    }

    return J;
}

// Evaluate the Hermite polynomial of order n at point x
double HermiteBasis::evaluatePolynomial(int n, double x) const {
    std::vector<double> coeffs = computeCoefficients(n);
    double result = 0.0;
    for (int i = 0; i <= n; ++i) {
        result += coeffs[i] * std::pow(x, i);
    }
    return result;
}

// Compute the derivative of the Hermite polynomial of order n at point x
double HermiteBasis::derivativePolynomial(int n, double x) const {
    std::vector<double> coeffs = computeCoefficients(n);
    double result = 0.0;
    for (int i = 1; i <= n; ++i) {
        result += i * coeffs[i] * std::pow(x, i - 1);
    }
    return result;
}

} // namespace Hermite
} // namespace Polynomials
