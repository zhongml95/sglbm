// polynomial.h
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "utils.h"
#include "quadrature_base.h"


// namespace olb {

// namespace uq {

namespace Polynomials {

// Abstract base class for all polynomial basis types
template <typename T>
class PolynomialBasis {
public:
    virtual ~PolynomialBasis() = default;

    // Compute polynomial coefficients of order n
    virtual std::vector<T> computeCoefficients(size_t n) const = 0;

    // Construct the Jacobi matrix of size n
    virtual std::vector<std::vector<T>> constructJacobiMatrix(size_t n) const = 0;

    // Evaluate the polynomial of order n at point x
    T evaluatePolynomial(size_t n, T x) const;

    // Compute the derivative of the polynomial of order n at point x
    T derivativePolynomial(size_t n, T x) const;

    // Dynamically create and return a quadrature object
    virtual std::shared_ptr<Quadrature::QuadratureBase<T>> getQuadrature(size_t nq, Quadrature::QuadratureMethod method) const = 0;

};

// Implement evaluatePolynomial using Horner's method for efficiency
template <typename T>
inline T PolynomialBasis<T>::evaluatePolynomial(size_t n, T x) const {
    std::vector<T> coeffs = computeCoefficients(n);

    // Evaluate polynomial using Horner's method
    T result = coeffs.back();
    for (size_t i = coeffs.size() - 1; i-- > 0;) {
        result = result * x + coeffs[i];
    }
    return result;
}

// Implement derivativePolynomial using Horner's method
template <typename T>
inline T PolynomialBasis<T>::derivativePolynomial(size_t n, T x) const {
    std::vector<T> coeffs = computeCoefficients(n);
    T result = 0.0;
    for (size_t i = coeffs.size() - 1; i > 0; --i) {
        result = result * x + i * coeffs[i];
    }
    return result;
}

// LegendreBasis class that inherits from PolynomialBasis
template <typename T>
class LegendreBasis : public PolynomialBasis<T> {
public:
    // Implement all pure virtual functions
    std::vector<T> computeCoefficients(size_t n) const override;
    std::vector<std::vector<T>> constructJacobiMatrix(size_t n) const override;
    std::shared_ptr<Quadrature::QuadratureBase<T>> getQuadrature(size_t nq, Quadrature::QuadratureMethod method) const override;
};

// HermiteBasis class that inherits from PolynomialBasis
template <typename T>
class HermiteBasis : public PolynomialBasis<T> {
public:
    // Implement all pure virtual functions
    std::vector<T> computeCoefficients(size_t n) const override;
    std::vector<std::vector<T>> constructJacobiMatrix(size_t n) const override;
    std::shared_ptr<Quadrature::QuadratureBase<T>> getQuadrature(size_t nq, Quadrature::QuadratureMethod method) const override;
};

} // namespace Polynomials

// } // namespace uq

// } // namespace olb

#endif // POLYNOMIAL_H
