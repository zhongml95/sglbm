// polynomial.h
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "utils.h"
#include "quadrature_base.h"

namespace Polynomials {

// Abstract base class for all polynomial basis types
class PolynomialBasis {
public:
    virtual ~PolynomialBasis() = default;

    // Compute polynomial coefficients of order n
    virtual std::vector<double> computeCoefficients(int n) const = 0;

    // Construct the Jacobi matrix of size n
    virtual std::vector<std::vector<double>> constructJacobiMatrix(int n) const = 0;

    // Evaluate the polynomial of order n at point x
    virtual double evaluatePolynomial(int n, double x) const = 0;

    // Compute the derivative of the polynomial of order n at point x
    virtual double derivativePolynomial(int n, double x) const = 0;

    // Compute quadrature weight at point x for polynomial of order n
    virtual double computeQuadratureWeight(int n, double x) const = 0;

    // Dynamically create and return a quadrature object
    virtual std::shared_ptr<Quadrature::QuadratureBase> getQuadrature(int order, Quadrature::QuadratureMethod method) const = 0;

};

// LegendreBasis class that inherits from PolynomialBasis
class LegendreBasis : public PolynomialBasis {
public:
    // Implement all pure virtual functions
    std::vector<double> computeCoefficients(int n) const override;
    std::vector<std::vector<double>> constructJacobiMatrix(int n) const override;
    double evaluatePolynomial(int n, double x) const override;
    double derivativePolynomial(int n, double x) const override;
    double computeQuadratureWeight(int n, double x) const override;
    std::shared_ptr<Quadrature::QuadratureBase> getQuadrature(int order, Quadrature::QuadratureMethod method) const override;
};

// HermiteBasis class that inherits from PolynomialBasis
class HermiteBasis : public PolynomialBasis {
public:
    // Implement all pure virtual functions
    std::vector<double> computeCoefficients(int n) const override;
    std::vector<std::vector<double>> constructJacobiMatrix(int n) const override;
    double evaluatePolynomial(int n, double x) const override;
    double derivativePolynomial(int n, double x) const override;
    double computeQuadratureWeight(int n, double x) const override;
    std::shared_ptr<Quadrature::QuadratureBase> getQuadrature(int order, Quadrature::QuadratureMethod method) const override;
};

} // namespace Polynomials

#endif // POLYNOMIAL_H
