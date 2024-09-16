// legendre_basis.h
#ifndef LEGENDRE_BASIS_H
#define LEGENDRE_BASIS_H

#include <vector>

namespace Polynomials {
namespace Legendre {

// LegendreBasis class that provides methods for computing Legendre polynomials and related operations
class LegendreBasis {
public:
    // Compute coefficients for the Legendre polynomial of order n
    std::vector<double> computeCoefficients(int n) const;

    // Construct the Jacobi matrix of size n for Legendre polynomials
    std::vector<std::vector<double>> constructJacobiMatrix(int n) const;

    // Evaluate the Legendre polynomial of order n at point x
    double evaluatePolynomial(int n, double x) const;

    // Compute the derivative of the Legendre polynomial of order n at point x
    double derivativePolynomial(int n, double x) const;
};

} // namespace Legendre
} // namespace Polynomials

#endif // LEGENDRE_BASIS_H
