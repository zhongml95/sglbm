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

#ifndef GAUSS_LEGENDRE_QUADRATURE_H
#define GAUSS_LEGENDRE_QUADRATURE_H

#include "quadrature/quadratureBase.h"
#include "matrixOperation.h"
#include "polynomials/legendreBasis.h"

namespace olb { namespace uq { namespace Quadrature {

template <typename T>
class GaussLegendreQuadrature : public QuadratureBase<T> {
public:
  GaussLegendreQuadrature(std::size_t nq,
                          QRMethod method = QRMethod::WilkinsonShiftQR)
  : nq_(nq),
    basis_(std::make_shared<olb::uq::polynomials::LegendreBasis<T>>(nq)),
    method_(method)
  {
    if (nq_ < 1) {
      throw std::invalid_argument("GaussLegendreQuadrature: nq must be >= 1");
    }
    performQRDecomposition_();
  }

  const std::vector<T>& getPoints()  const override { return points_;  }
  const std::vector<T>& getWeights() const override { return weights_; }

private:
  std::size_t                                        nq_;
  std::shared_ptr<olb::uq::polynomials::LegendreBasis<T>> basis_;
  QRMethod                                           method_;
  std::vector<T>                                     points_;
  std::vector<T>                                     weights_;

  // Jacobi matrix for Legendre polynomials (symmetric tridiagonal)
  std::vector<std::vector<T>> constructJacobiMatrix_(std::size_t n) const {
    if (n <= 1) {
      throw std::invalid_argument("constructJacobiMatrix_: n must be > 1");
    }
    std::vector<std::vector<T>> J(n, std::vector<T>(n, T(0)));
    for (std::size_t i = 1; i < n; ++i) {
      // a_i = i / sqrt(4 i^2 - 1)
      const T a = static_cast<T>(i) /
                  std::sqrt(static_cast<T>(4 * i * i - 1));
      J[i][i - 1] = a;
      J[i - 1][i] = a;
    }
    return J;
  }

  void performQRDecomposition_() {
    // Eigenvalues of Jacobi matrix = Gauss–Legendre nodes
    if (nq_ == 1) {
      points_.assign(1, T(0)); // single node at 0
    } else {
      auto J = constructJacobiMatrix_(nq_);
      if (method_ == QRMethod::HouseholderQR) {
        points_ = olb::uq::MatrixOperations::HouseholderQRDecomposition(J);
      } else {
        points_ = olb::uq::MatrixOperations::WilkinsonShiftQRDecomposition(J);
      }
    }
    computeWeights_();
  }

  void computeWeights_() {
    // w_i = 2 / [(1 - x_i^2) * (P'_n(x_i))^2]
    weights_.resize(points_.size());
    for (std::size_t i = 0; i < points_.size(); ++i) {
      const T x  = points_[i];
      const T dP = basis_->derivativePolynomial(nq_, x); // cached in basis
      weights_[i] = static_cast<T>(2) /
                    ((static_cast<T>(1) - x * x) * dP * dP);
    }
    // Do NOT renormalize; formula already yields correct normalization (sum=2).
  }
};

}}} // ns

#endif // GAUSS_LEGENDRE_QUADRATURE_H