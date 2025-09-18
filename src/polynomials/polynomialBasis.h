/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

// polynomial.h
#ifndef POLYNOMIAL_BASIS_H
#define POLYNOMIAL_BASIS_H


namespace olb {

namespace uq {

namespace polynomials {

  
enum class BasisKind {
  Unknown = 0,
  Legendre,
  Hermite,
  // Chebyshev
};

// Abstract base class for polynomial bases.
// Coefficient convention: coeffs[i] is the coefficient of x^i (ascending degree).
template <typename T>
class polynomialBasis {
public:
  explicit polynomialBasis(std::size_t order)
  : order_(order), coeffs_(order + 1) {}

  virtual ~polynomialBasis() = default;

  // --- Identification API ---
  virtual std::string name() const = 0;
  virtual BasisKind   kind() const noexcept { return BasisKind::Unknown; }

  // Return coefficients of the degree-n polynomial in ascending powers.
  virtual std::vector<T> computeCoefficients(std::size_t n) const = 0;

  // Horner evaluation using ascending coeffs.
  T evaluatePolynomial(std::size_t n, T x) const {
    const auto coeffs = coefficients(n);
    if (coeffs.empty()) return T(0);

    T acc = coeffs.back();
    // Post-decrement trick to include coeffs[0] at the end.
    for (std::size_t i = coeffs.size() - 1; i-- > 0; ) {
      acc = acc * x + coeffs[i];
    }
    return acc;
  }

  // Evaluate derivative at x (d/dx of sum c_i x^i is sum i c_i x^{i-1}).
T derivativePolynomial(std::size_t n, T x) const {
    const auto& coeffs = coefficients(n);  // <-- use precomputed coeffs_[n]
    if (coeffs.size() <= 1) return T(0);

    T acc = T(0);
    for (std::size_t i = coeffs.size() - 1; i > 0; --i) {
        acc = acc * x + static_cast<T>(i) * coeffs[i];
    }
    return acc;
}

  // Access precomputed coefficients for all polynomials up to order_.
  const std::vector<std::vector<T>>& coefficients() const { return coeffs_; }
  const std::vector<T>&              coefficients(std::size_t n) const {
    if (n > order_) {
      throw std::out_of_range("Requested polynomial order exceeds initialized order.");
    }
    return coeffs_[n];
  }

protected:
  // Call this from the DERIVED constructor to fill coeffs_ once.
  void buildAll() {
    coeffs_.resize(order_ + 1);
    for (std::size_t n = 0; n <= order_; ++n) {
      coeffs_[n] = computeCoefficients(n); // <-- reuse your function
    }
  }

  std::size_t                         order_{0};
  std::vector<std::vector<T>>         coeffs_; // coeffs_[n][i] is coeff of x^i

};


} // namespace Polynomials

} // namespace uq

} // namespace olb

#endif // POLYNOMIAL_BASIS_H
