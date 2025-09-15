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

#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

#include "quadrature/quadratureBase.h"
#include "matrixOperation.h"



namespace olb {

namespace uq {

namespace Quadrature {


template <typename T, typename PolynomialBasis>
class GaussQuadrature : public QuadratureBase<T> {
public:
  GaussQuadrature(std::size_t nq, QRMethod method = QRMethod::WilkinsonShiftQR)
      : nq(nq)
      , basis(std::make_shared<PolynomialBasis>())
      , method(method)
  {
    performQRDecomposition();
  }

  const std::vector<T>& getPoints() const override { return points; }

  const std::vector<T>& getWeights() const override { return weights; }

private:
  std::size_t                      nq;
  std::shared_ptr<PolynomialBasis> basis;
  QRMethod                         method;
  std::vector<T>                   points;
  std::vector<T>                   weights;

  std::vector<std::vector<T>> constructJacobiMatrix(std::size_t n) const
  {
    if (n <= 1) {
      throw std::invalid_argument("The order (n) must be greater than 1.");
    }

    std::vector<std::vector<T>> J(n, std::vector<T>(n, 0.0));
    if constexpr (std::is_same_v<PolynomialBasis, olb::uq::Polynomials::LegendreBasis<T>>) {
      // Fill the Jacobi matrix according to the recurrence relation for Legendre polynomials
      for (std::size_t i = 1; i < n; ++i) {
        T a         = i / std::sqrt(4.0 * i * i - 1.0);
        J[i][i - 1] = a;
        J[i - 1][i] = a;
      }
    }
    else if constexpr (std::is_same_v<PolynomialBasis, olb::uq::Polynomials::HermiteBasis<T>>) {
      // Fill the Jacobi matrix according to the recurrence relation for Hermite polynomials
      for (std::size_t i = 1; i < n; ++i) {
        // The recurrence relation coefficients for Hermite polynomials
        if (i > 0) {
          J[i][i - 1] = std::sqrt(i);
          J[i - 1][i] = J[i][i - 1];
        }
      }
    }
    else {
      throw std::invalid_argument("Unsupported basis in Quadrature::constructJacobiMatrix()");
    }

    return J;
  }

  void performQRDecomposition()
  {
    if (nq > 1) {
      // Step 1: Construct the Jacobi matrix using the basis
      auto J = constructJacobiMatrix(nq);

      // Step 2: Perform the appropriate QR decomposition
      if (method == QRMethod::HouseholderQR) {
        // std::cout << "Using HouseholderQR method." << std::endl;
        points = olb::uq::MatrixOperations::HouseholderQRDecomposition(J); // Perform Householder QR and get eigenvalues
      }
      else if (method == QRMethod::WilkinsonShiftQR) {
        // std::cerr << "Using Wilkinson's Shift QR method." << std::endl;
        points = olb::uq::MatrixOperations::WilkinsonShiftQRDecomposition(J); // Perform Wilkinson Shift QR
      }
      else {
        // std::cerr << "Warning: Unsupported method. Defaulting to Wilkinson's Shift QR." << std::endl;
        points = olb::uq::MatrixOperations::WilkinsonShiftQRDecomposition(J); // Default to Wilkinson Shift QR
      }
    }
    else if (nq == 1) {
      points.resize(1);
      points[0] = 0.0; // For nq = 1, the point is always 0
    }
    else {
      throw std::invalid_argument("Quadrature points number (nq) must be greater than or equal to 1.");
    }

    // Step 3: Compute quadrature weights
    computeWeights();
  }

  void computeWeights()
  {
    weights.resize(points.size());

    // If the “PolynomialBasis” is actually LegendreBasis<T>:
    if constexpr (std::is_same_v<PolynomialBasis, olb::uq::Polynomials::LegendreBasis<T>>) {
      for (std::size_t i = 0; i < points.size(); ++i) {
        T x        = points[i];
        T Pn_prime = basis->derivativePolynomial(nq, x);
        weights[i] = 2.0 / ((1.0 - x * x) * Pn_prime * Pn_prime);
      }
    }
    // Else if it's HermiteBasis<T>:
    else if constexpr (std::is_same_v<PolynomialBasis, olb::uq::Polynomials::HermiteBasis<T>>) {
      for (std::size_t i = 0; i < points.size(); ++i) {
        T x         = points[i];
        T Hn_minus1 = basis->evaluatePolynomial(nq - 1, x);
        weights[i]  = std::pow(2, nq - 1) * std::tgamma(nq + 1) * std::sqrt(M_PI) / std::pow(nq * Hn_minus1, 2);
      }
    }
    else {
      throw std::invalid_argument("Unsupported basis in Quadrature::computeWeights()");
    }

    // (optionally normalize)
    T sum = std::accumulate(weights.begin(), weights.end(), T(0));
    for (auto& w : weights) {
      w /= sum;
    }
  }

};

} // namespace Quadrature

} // namespace uq
} // namespace olb


#endif // GAUSSQUADRATURE_H