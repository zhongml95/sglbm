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

#ifndef GAUSS_HERMITE_QUADRATURE_H
#define GAUSS_HERMITE_QUADRATURE_H

#include "quadrature/quadratureBase.h"
#include "matrixOperation.h"
#include "polynomials/hermiteBasis.h"

namespace olb { namespace uq { namespace Quadrature {

template <typename T>
class GaussHermiteQuadrature : public QuadratureBase<T> {
public:
  GaussHermiteQuadrature(std::size_t nq,
                         QRMethod method = QRMethod::WilkinsonShiftQR)
  : nq_(nq), basis_(std::make_shared<olb::uq::polynomials::HermiteBasis<T>>(nq)),
    method_(method)
  {
    performQRDecomposition();
  }

  const std::vector<T>& getPoints()  const override { return points_;  }
  const std::vector<T>& getWeights() const override { return weights_; }

private:
  std::size_t                                         nq_;
  std::shared_ptr<olb::uq::polynomials::HermiteBasis<T>> basis_;
  QRMethod                                            method_;
  std::vector<T>                                      points_;
  std::vector<T>                                      weights_;

  std::vector<std::vector<T>> constructJacobiMatrix(std::size_t n) const {
    if (n <= 1) throw std::invalid_argument("n must be > 1");
    // Symmetric tridiagonal for probabilists’ Hermite: beta_i = sqrt(i)
    std::vector<std::vector<T>> J(n, std::vector<T>(n, T(0)));
    for (std::size_t i = 1; i < n; ++i) {
      T a = std::sqrt(static_cast<T>(i));
      J[i][i-1] = a;
      J[i-1][i] = a;
    }
    return J;
  }

  void performQRDecomposition() {
    if (nq_ == 0) throw std::invalid_argument("nq must be >= 1");

    if (nq_ == 1) {
      points_.assign(1, T(0));
    } else {
      auto J = constructJacobiMatrix(nq_);
      if (method_ == QRMethod::HouseholderQR) {
        points_ = olb::uq::MatrixOperations::HouseholderQRDecomposition(J);
      } else {
        points_ = olb::uq::MatrixOperations::WilkinsonShiftQRDecomposition(J);
      }
    }
    computeWeights();
  }

  void computeWeights() {
    weights_.resize(points_.size());
    // Your formula (matches probabilists’ convention in your basis):
    // w_i = 2^{n-1} n! sqrt(pi) / [ n * He_{n-1}(x_i) ]^2
    for (std::size_t i = 0; i < points_.size(); ++i) {
      T x          = points_[i];
      T Hn_minus_1 = basis_->evaluatePolynomial(nq_ - 1, x);
      // NOTE: std::tgamma(n+1) = n!
      T numer = std::pow(static_cast<T>(2), static_cast<T>(nq_ - 1))
              * std::tgamma(static_cast<T>(nq_) + 1)
              * std::sqrt(M_PI);
      T denom = static_cast<T>(nq_) * Hn_minus_1;
      weights_[i] = numer / (denom * denom);
    }

    // Normalize weights to sum to 1 (integral of exp(-x^2) over R)
    const T intervalLength = static_cast<T>(2);
    for (auto& w : weights_) {
      w /= intervalLength;
    }
  }
};

}}} // ns

#endif // GAUSS_HERMITE_QUADRATURE_H