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

// quadrature.h
#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "filesIO.h"
#include "quadratureBase.h"
#include "polynomial.h"
#include "matrixOperation.h"
#include "genzKeisterLibrary.h"
#include "fft.h"

namespace olb {

namespace uq {

namespace Quadrature {
template <typename>
struct always_false : std::false_type {};


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
      static_assert(always_false<PolynomialBasis>::value, "Unsupported basis in Quadrature::constructJacobiMatrix()");
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
      static_assert(always_false<PolynomialBasis>::value, "Unsupported basis in Quadrature::computeWeights()");
    }

    // (optionally normalize)
    T sum = std::accumulate(weights.begin(), weights.end(), T(0));
    for (auto& w : weights) {
      w /= sum;
    }
  }

};

template <typename T>
class GenzKeister16 : public QuadratureBase<T> {
public:
  GenzKeister16(std::size_t order) {
    std::size_t level = getLevelFromOrder(order);

    std::vector<T> rawPoints, rawWeights;
    GenzKeisterLibrary<T>::getPointsAndWeights(level, rawPoints, rawWeights);

    points = std::move(rawPoints);
    weights = std::move(rawWeights);
  }

  const std::vector<T>& getPoints() const override { return points; }
  const std::vector<T>& getWeights() const override { return weights; }

private:
  std::vector<T> points;
  std::vector<T> weights;

  static std::size_t getLevelFromOrder(std::size_t order) {
    static const std::vector<std::size_t> RULES_16 = {1, 3, 7, 9, 17, 19, 31};
    if (order >= RULES_16.size()) {
      throw std::out_of_range("Genz-Keister order too large. Max supported index is " +
                              std::to_string(RULES_16.size() - 1));
    }
    return RULES_16[order];
  }
};

template <typename T>
class GenzKeister18 : public QuadratureBase<T> {
public:
  GenzKeister18(std::size_t order) {
    std::size_t level = getLevelFromOrder(order);

    std::vector<T> rawPoints, rawWeights;
    GenzKeisterLibrary<T>::getPointsAndWeights(level, rawPoints, rawWeights);

    points = std::move(rawPoints);
    weights = std::move(rawWeights);
  }

  const std::vector<T>& getPoints() const override { return points; }
  const std::vector<T>& getWeights() const override { return weights; }

private:
  std::vector<T> points;
  std::vector<T> weights;

  static std::size_t getLevelFromOrder(std::size_t order) {
    static const std::vector<std::size_t> RULES_18 = {1, 3, 9, 19, 37};
    if (order >= RULES_18.size()) {
      throw std::out_of_range("Genz-Keister order too large. Max supported index is " +
                              std::to_string(RULES_18.size() - 1));
    }
    return RULES_18[order];
  }
};

template <typename T>
class GenzKeister22 : public QuadratureBase<T> {
public:
  GenzKeister22(std::size_t order) {
    std::size_t level = getLevelFromOrder(order);

    std::vector<T> rawPoints, rawWeights;
    GenzKeisterLibrary<T>::getPointsAndWeights(level, rawPoints, rawWeights);

    points = std::move(rawPoints);
    weights = std::move(rawWeights);
  }

  const std::vector<T>& getPoints() const override { return points; }
  const std::vector<T>& getWeights() const override { return weights; }

private:
  std::vector<T> points;
  std::vector<T> weights;

  static std::size_t getLevelFromOrder(std::size_t order) {
    static const std::vector<std::size_t> RULES_22 = {1, 3, 9, 19, 41};
    if (order >= RULES_22.size()) {
      throw std::out_of_range("Genz-Keister order too large. Max supported index is " +
                              std::to_string(RULES_22.size() - 1));
    }
    return RULES_22[order];
  }
};

template <typename T>
class GenzKeister24 : public QuadratureBase<T> {
public:
  GenzKeister24(std::size_t order) {
    std::size_t level = getLevelFromOrder(order);

    std::vector<T> rawPoints, rawWeights;
    GenzKeisterLibrary<T>::getPointsAndWeights(level, rawPoints, rawWeights);

    points = std::move(rawPoints);
    weights = std::move(rawWeights);
  }

  const std::vector<T>& getPoints() const override { return points; }
  const std::vector<T>& getWeights() const override { return weights; }

private:
  std::vector<T> points;
  std::vector<T> weights;

  static std::size_t getLevelFromOrder(std::size_t order) {
    static const std::vector<std::size_t> RULES_24 = {1, 3, 9, 19, 43};
    if (order >= RULES_24.size()) {
      throw std::out_of_range("Genz-Keister order too large. Max supported index is " +
                              std::to_string(RULES_24.size() - 1));
    }
    return RULES_24[order];
  }
};


template <typename T>
class ClenshawCurtisQuadrature : public QuadratureBase<T> {
public:
    ClenshawCurtisQuadrature(std::size_t order) {
        if (order == 0) {
            points = {T(0.0)};
            weights = {T(1.0)};
            return;
        } else if (order == 1) {
            points = {T(-1.0), T(1.0)};
            weights = {T(0.5), T(0.5)};
            return;
        }

        std::vector<T> theta(order + 1);
        for (std::size_t k = 0; k <= order; ++k) {
            theta[k] = (order - k) * M_PI / order;
        }
        points.resize(order + 1);
        for (std::size_t k = 0; k <= order; ++k) {
            points[k] = std::cos(theta[k]);
        }


        // Match Python: steps = np.arange(1, order, 2)
        std::vector<T> steps;
        for (std::size_t k = 1; k < order; k += 2) {
            steps.push_back(static_cast<T>(k));
        }

        std::size_t length = steps.size();
        std::size_t remains = order - length;

        std::vector<T> beta;

        // Step 1: 2.0 / (k * (k - 2)) for each step
        for (std::size_t i = 0; i < steps.size(); ++i) {
            T k = steps[i];
            beta.push_back(2.0 / (k * (k - 2.0)));
        }

        // Step 2: Append [1.0 / steps.back()]
        beta.push_back(1.0 / steps.back());

        // Step 3: Append remains zeros
        for (std::size_t i = 0; i < remains; ++i) {
            beta.push_back(0.0);
        }

        std::size_t N = beta.size();
        std::vector<T> new_beta(N - 1);

        for (std::size_t i = 0; i < N - 1; ++i) {
            new_beta[i] = -beta[i] - beta[N - 1 - i];
        }

        beta = new_beta;

        std::vector<T> gamma(order, -1.0);
        gamma[length] += static_cast<T>(order);
        gamma[remains] += static_cast<T>(order);
        T denom = static_cast<T>(order * order - 1 + (order % 2));
        for (auto& g : gamma) {
            g /= denom;
        }

        assert(beta.size() == gamma.size());

        std::vector<T> combined(order);
        for (std::size_t i = 0; i < order; ++i) {
          combined[i] = beta[i] + gamma[i];
        }

      std::vector<std::complex<T>> combined_complex(combined.size(), 0.0);
      for (std::size_t i = 0; i < combined.size(); ++i) {
          combined_complex[i] = std::complex<T>(combined[i], 0.0);
      }

      olb::uq::fft::ifft_dispatch<T>(combined_complex);

      std::vector<T> half_weights;
      for (std::size_t i = 0; i < (combined_complex.size()/2 + 1); ++i) {
          half_weights.push_back(std::real(combined_complex[i]));
      }

      // Step 2: Mirror and divide by 2
      weights = half_weights;

      int start = static_cast<int>(half_weights.size()) - 2 + (order % 2);
      for (int i = start; i >= 0; --i) {
          weights.push_back(half_weights[i]);
      }

      for (std::size_t i = 0; i < weights.size(); ++i) {
        weights[i] /= 2.0;
      }
    }

    const std::vector<T>& getPoints() const override { return points; }
    const std::vector<T>& getWeights() const override { return weights; }

private:
    std::vector<T> points;
    std::vector<T> weights;
};


// Then define makeQuadrature()
template <typename T>
std::unique_ptr<QuadratureBase<T>> makeQuadrature(
    const std::shared_ptr<olb::uq::Polynomials::PolynomialBasis<T>>& basis,
    std::size_t order,
    QuadratureMethod method)
{
  using namespace olb::uq::Polynomials;

  if (auto legendre = std::dynamic_pointer_cast<LegendreBasis<T>>(basis)) {
    switch (method) {
      case QuadratureMethod::GaussQuadrature:
        return std::make_unique<GaussQuadrature<T, LegendreBasis<T>>>(order);
      case QuadratureMethod::ClenshawCurtis:
        return std::make_unique<ClenshawCurtisQuadrature<T>>(order);
      default:
        throw std::runtime_error("Unsupported quadrature method for Legendre basis.");
    }
  }
  else if (auto hermite = std::dynamic_pointer_cast<HermiteBasis<T>>(basis)) {
    switch (method) {
      case QuadratureMethod::GaussQuadrature:
        return std::make_unique<GaussQuadrature<T, HermiteBasis<T>>>(order);

      case QuadratureMethod::GenzKeister16:
        if (order > 6) {
          throw std::runtime_error("Genz-Keister-16 index too large: " + std::to_string(order));
        }
        return std::make_unique<GenzKeister16<T>>(order);

      case QuadratureMethod::GenzKeister18:
        if (order > 4) {
          throw std::runtime_error("Genz-Keister-18 index too large: " + std::to_string(order));
        }
        return std::make_unique<GenzKeister18<T>>(order);

      case QuadratureMethod::GenzKeister22:
        if (order > 4) {
          throw std::runtime_error("Genz-Keister-22 index too large: " + std::to_string(order));
        }
        return std::make_unique<GenzKeister22<T>>(order);

      case QuadratureMethod::GenzKeister24:
        if (order > 4) {
          throw std::runtime_error("Genz-Keister-24 index too large: " + std::to_string(order));
        }
        return std::make_unique<GenzKeister24<T>>(order);

      default:
        throw std::runtime_error("Unknown quadrature method for Hermite basis.");
    }
  }
  else {
    throw std::runtime_error("Unsupported basis in makeQuadrature");
  }
}




} // namespace Quadrature

} // namespace uq

} // namespace olb

#endif // QUADRATURE_H
