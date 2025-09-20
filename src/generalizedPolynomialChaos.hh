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

// generalized_polynomial_chaos.hh
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_HH
#define GENERALIZED_POLYNOMIAL_CHAOS_HH

#include "generalizedPolynomialChaos.h"
#include "matrixOperation.h"
#include "quadrature/quadrature.h"
#include "stochasticCollocationGrid.h"

#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <functional>

// Include the polynomial basis and quadrature headers
// #include "quadrature.h"

// Namespace aliases

namespace olb {

namespace uq {

// Constructor
template <typename T>
GeneralizedPolynomialChaos<T>::GeneralizedPolynomialChaos( std::size_t _order, std::size_t _nq,
                                                           const olb::uq::StochasticCollocationGrid<T>& _grid,
                                                           std::vector<std::shared_ptr<olb::uq::polynomials::polynomialBasis<T>>> _polynomialBases)
    : order(_order)
    , nq(_nq)
    , polynomialBases(_polynomialBases)
    , grid(_grid)
{
  totalNq = grid.totalNq;
  points  = grid.points;
  weights = grid.weightsMultiplied;

  randomNumberDimension = _polynomialBases.size();

  // Calculate multi-indices
  calculateMultiIndices(randomNumberDimension, order, inds);
  No = inds.size();

  phiRan.resize(totalNq * No, 0.0);
  phiRan_T.resize(totalNq * No, 0.0);
  t2Product.resize(No, 0.0);
  t2Product_inv.resize(No, 0.0);
  t3Product.resize(No * No * No, 0.0);

  // Evaluate polynomials at quadrature points
  evaluatePhiRan();

  // Compute tensors
  computeTensors();
}


// Evaluate n_order polynomial at point k
template <typename T>
T GeneralizedPolynomialChaos<T>::evaluate(std::size_t n_order, std::size_t k)
{
  T result = 1.0;
  for (std::size_t phi_i = 0; phi_i < randomNumberDimension; ++phi_i) {
    result *= phi1D[phi_i][inds[n_order][phi_i]][k];
  }
  return result;
}

// Evaluate phiRan matrix
template <typename T>
void GeneralizedPolynomialChaos<T>::evaluatePhiRan()
{
  phi1D.resize(randomNumberDimension);
  for (std::size_t d = 0; d < randomNumberDimension; ++d) {
    auto maxDeg = order; // total degree -> per-dim degree max is 'order'
    phi1D[d].assign(maxDeg+1, std::vector<T>(totalNq));
    for (std::size_t k = 0; k < totalNq; ++k) {
      for (std::size_t p = 0; p <= maxDeg; ++p) {
        phi1D[d][p][k] = polynomialBases[d]->evaluatePolynomial(p, points[k][d]);
      }
    }
  }
  for (std::size_t k = 0; k < totalNq; ++k) {
    for (std::size_t i = 0; i < No; ++i) {
      phiRan[k * No + i] = evaluate(i, k);
      phiRan_T[i * totalNq + k] = phiRan[k * No + i];
    }
  }
}

// Helper functions
template <typename T>
void GeneralizedPolynomialChaos<T>::calculateMultiIndices(std::size_t d, std::size_t n,
                                                          std::vector<std::vector<std::size_t>>& indices)
{

  std::vector<std::size_t> index(d, 0);

  std::function<void(std::size_t, std::size_t, std::size_t)> recursiveFunction = [&](std::size_t pos, std::size_t sum,
                                                                                     std::size_t maxOrder) {
    if (pos == d - 1) {
      index[pos] = maxOrder - sum;
      indices.push_back(index);
      return;
    }

    for (std::size_t i = 0; i <= maxOrder - sum; ++i) {
      index[pos] = i;
      recursiveFunction(pos + 1, sum + i, maxOrder);
    }
  };

  for (std::size_t order = 0; order <= n; ++order) {
    recursiveFunction(0, 0, order);
  }
}

// Compute tensors (t2Product and t3Product) with sparse t3
template <typename T>
void GeneralizedPolynomialChaos<T>::computeTensors()
{
  // small threshold to prune numerical near-zeros in triple products
  const T eps = static_cast<T>(1e-14);

    // t2 diagonal inner products
    std::fill(t2Product.begin(),     t2Product.end(),     T(0));
    std::fill(t2Product_inv.begin(), t2Product_inv.end(), T(0));

    for (std::size_t i = 0; i < No; ++i) {
      T acc = T(0);
      for (std::size_t m = 0; m < totalNq; ++m) {
        const T v = phiRan[m * No + i];
        acc += v * v * weights[m];
      }
      t2Product[i] = acc;
    }
    for (std::size_t i = 0; i < No; ++i) {
      const T denom = t2Product[i];
      t2Product_inv[i] = (std::abs(denom) > eps) ? T(1) / denom : T(0);
    }

    // Build sparse t3
    t3RowPtr.assign(No * No + 1, 0);
    t3KIndex.clear();
    t3Product.clear();

    for (std::size_t i = 0; i < No; ++i) {
      for (std::size_t j = 0; j < No; ++j) {
        const std::size_t row = i * No + j;
        t3RowPtr[row] = t3KIndex.size();

        for (std::size_t k = 0; k < No; ++k) {
          T sum = T(0);
          for (std::size_t m = 0; m < totalNq; ++m) {
            const std::size_t base = m * No;
            sum += phiRan[base + i] * phiRan[base + j] * phiRan[base + k] * weights[m];
          }
          if (std::abs(sum) > eps) {
            t3KIndex.push_back(k);
            t3Product.push_back(sum); // keep "t3Product" as values
          }
        }
      }
    }
    t3RowPtr[No * No] = t3KIndex.size();
}


// Transformation functions
template <typename T>
void GeneralizedPolynomialChaos<T>::chaosToRandom(const std::vector<T>& chaosCoefficients,
                                                  std::vector<T>&       randomVariables)
{
  randomVariables.resize(totalNq, 0.0);

  for (std::size_t k = 0; k < totalNq; ++k) {
    auto startIt       = phiRan.begin() + k * No;
    randomVariables[k] = std::inner_product(chaosCoefficients.begin(), chaosCoefficients.end(), startIt, 0.0);
  }
}

template <typename T>
void GeneralizedPolynomialChaos<T>::randomToChaos(const std::vector<T>& randomVariables,
                                                  std::vector<T>&       chaosCoefficients)
{
  chaosCoefficients.resize(No, 0.0);
  std::vector<T> weightedRandomVariables(totalNq);

  // Compute weighted random variables
  for (std::size_t k = 0; k < totalNq; ++k) {
    weightedRandomVariables[k] = weights[k] * randomVariables[k];
  }

  // Compute chaos coefficients
  for (std::size_t i = 0; i < No; ++i) {
    auto startIt = phiRan_T.begin() + i * totalNq;
    chaosCoefficients[i] =
        std::inner_product(weightedRandomVariables.begin(), weightedRandomVariables.end(), startIt, 0.0);
    chaosCoefficients[i] *= t2Product_inv[i];
  }
}

// Chaos operations
template <typename T>
void GeneralizedPolynomialChaos<T>::chaosProduct(const std::vector<T>& chaos1, const std::vector<T>& chaos2,
                                                 std::vector<T>& product)
{
  product.assign(No, T(0));

  for (std::size_t i = 0; i < No; ++i) {
    T coeffSum = T(0);

    // For each j, walk the sparse row (i, j, :)
    for (std::size_t j = 0; j < No; ++j) {
      const T c1j = chaos1[j];
      if (c1j == T(0)) continue;

      const std::size_t row     = i * No + j;           // CSR row id for (i,j)
      const std::size_t begin   = t3RowPtr[row];
      const std::size_t end     = t3RowPtr[row + 1];

      // Accumulate c1[j] * c2[k] * G_{i,j,k} over the row's nonzeros
      for (std::size_t p = begin; p < end; ++p) {
        const std::size_t k = t3KIndex[p];
        coeffSum += c1j * chaos2[k] * t3Product[p];     // t3Product holds sparse values
      }
    }

    product[i] = coeffSum * t2Product_inv[i];
  }
}

template <typename T>
void GeneralizedPolynomialChaos<T>::chaosSum(const std::vector<T>& chaos1, const std::vector<T>& chaos2,
                                             std::vector<T>& sum)
{
  sum.resize(No);
  for (std::size_t i = 0; i < No; ++i) {
    sum[i] = chaos1[i] + chaos2[i];
  }
}

// Statistical moments
template <typename T>
T GeneralizedPolynomialChaos<T>::mean(const std::vector<T>& chaosCoefficients)
{
  return chaosCoefficients[0];
}

template <typename T>
T GeneralizedPolynomialChaos<T>::std(const std::vector<T>& chaosCoefficients)
{
  T variance = 0.0;
  for (std::size_t i = 1; i < No; ++i) {
    variance += t2Product[i] * chaosCoefficients[i] * chaosCoefficients[i];
  }
  return std::sqrt(variance);
}

template <typename T>
void GeneralizedPolynomialChaos<T>::convert2affinePCE(const Distribution<T>& distribution, std::vector<T>& chaos)
{
  switch (distribution.type) {
  case DistributionType::Uniform: {
    T a1     = 0.5 * (distribution.param1 + distribution.param2);
    T a2     = 0.5 * (distribution.param2 - distribution.param1);
    chaos[0] = a1;
    chaos[1] = a2;
    break;
  }
  case DistributionType::Normal: {
    chaos[0] = distribution.param1;
    chaos[1] = distribution.param2;
    break;
  }
  // Add cases for other distributions
  default:
    throw std::runtime_error("Unsupported distribution type for convert2affinePCE.");
  }
}

// Getters
template <typename T>
std::size_t GeneralizedPolynomialChaos<T>::getPolynomialsOrder() const
{
  return No;
}

template <typename T>
std::size_t GeneralizedPolynomialChaos<T>::getQuadraturePointsNumber() const
{
  return totalNq;
}


template <typename T>
void GeneralizedPolynomialChaos<T>::getPointsAndWeights(std::vector<std::vector<T>>& _points,
                                                        std::vector<T>& _weights)
{
  _points  = this->points;
  _weights = this->weights;
}


template <typename T>
std::vector<std::vector<T>> GeneralizedPolynomialChaos<T>::getStochasticCollocationSample()
{
  return points;
}

template <typename T>
std::vector<T> GeneralizedPolynomialChaos<T>::getWeightsMultiplied() const
{
  return weights;
}

template <typename T>
void GeneralizedPolynomialChaos<T>::getTensors(std::vector<T>& t2Product, std::vector<T>& t2Product_inv,
                                               std::vector<T>& t3Product)
{
  t2Product     = this->t2Product;
  t2Product_inv = this->t2Product_inv;
  t3Product     = this->t3Product;
}

template <typename T>
std::shared_ptr<polynomials::polynomialBasis<T>>
GeneralizedPolynomialChaos<T>::getPolynomialBasis(std::size_t dimension) const
{
  if (dimension >= randomNumberDimension) {
    throw std::out_of_range("Dimension is out of bounds");
  }

  // Cast the void pointer back to the correct polynomial basis type
  return polynomialBases[dimension];
}

template <typename T>
std::vector<std::vector<std::size_t>> GeneralizedPolynomialChaos<T>::getMultiIndices() const
{
  return inds;
}

template <typename T>
void GeneralizedPolynomialChaos<T>::getPhiRan(std::vector<T>& phiRan)
{
  phiRan = this->phiRan;
}

template <typename T>
void GeneralizedPolynomialChaos<T>::getCoefficients(std::vector<std::vector<std::vector<T>>>& polynomialCoeffs)
{
  polynomialCoeffs = this->coefficients;
}

} // namespace uq

} // namespace olb

#endif // GENERALIZED_POLYNOMIAL_CHAOS_HH
