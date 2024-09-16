// hermite_basis.h
#ifndef HERMITE_BASIS_H
#define HERMITE_BASIS_H

#include <vector>

namespace Polynomials {
namespace Hermite {

// HermiteBasis class that provides methods for computing Hermite polynomials and related operations
class HermiteBasis {
public:
    // Compute coefficients for the Hermite polynomial of order n
    std::vector<double> computeCoefficients(int n) const;

    // Construct the Jacobi matrix of size n for Hermite polynomials
    std::vector<std::vector<double>> constructJacobiMatrix(int n) const;

    // Evaluate the Hermite polynomial of order n at point x
    double evaluatePolynomial(int n, double x) const;

    // Compute the derivative of the Hermite polynomial of order n at point x
    double derivativePolynomial(int n, double x) const;
};

} // namespace Hermite
} // namespace Polynomials

#endif // HERMITE_BASIS_H
