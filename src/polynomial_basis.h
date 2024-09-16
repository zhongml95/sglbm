// polynomial_basis.h
#ifndef POLYNOMIAL_BASIS_H
#define POLYNOMIAL_BASIS_H

#include <vector>

namespace Polynomials {

// Documentation of the required interface for a PolynomialBasis
// Any class used as a PolynomialBasis with the templated Quadrature class
// must implement the following methods:

/*
class PolynomialBasis {
public:
    // Compute polynomial coefficients of order n
    std::vector<double> computeCoefficients(int n) const;

    // Construct the Jacobi matrix of size n
    std::vector<std::vector<double>> constructJacobiMatrix(int n) const;

    // Evaluate the polynomial of order n at point x
    double evaluatePolynomial(int n, double x) const;

    // Compute the derivative of the polynomial of order n at point x
    double derivativePolynomial(int n, double x) const;
};
*/

// This file serves as a guideline for implementing PolynomialBasis classes.
// It is not meant to be instantiated directly.

} // namespace Polynomials

#endif // POLYNOMIAL_BASIS_H
