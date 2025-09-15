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

#ifndef SMOLYAK_H
#define SMOLYAK_H

#include <vector>
#include <memory>
#include <numeric>
#include <map>
#include <cmath>
#include <iostream>
#include <functional>
#include <algorithm>

#include "matrixOperation.h"
#include "polynomial.h"
#include "quadratureBase.h"
#include "quadrature.h"  // for get1DRule to call polynomialBasis->getQuadrature()

namespace olb {
namespace uq {
namespace smolyak {


// helper: fast integer binomial coefficient  C(n,k)
static constexpr std::size_t binom(std::size_t n, std::size_t k)
{
  if (k > n)              return 0;
  if (k == 0 || k == n)   return 1;
  if (k > n - k)          k = n - k;          // use symmetry
  std::size_t res = 1;
  for (std::size_t i = 1; i <= k; ++i)
  {
    res = res * (n - k + i) / i;
  }
  return res;
}

template <typename T>
void get1DRule(
    int level_i,
    int dim,
    const std::vector<std::shared_ptr<olb::uq::Polynomials::PolynomialBasis<T>>>& polynomialBases,
    olb::uq::Quadrature::QuadratureMethod quadratureMethod,
    std::vector<T>& points,
    std::vector<T>& weights)
{
  std::size_t m = static_cast<std::size_t>(level_i);
    if (quadratureMethod == olb::uq::Quadrature::QuadratureMethod::ClenshawCurtis) {
      // Chaospy’s growth=True behaviour:  order' = 2^ℓ   (ℓ>0)
      m = (level_i == 0) ? 0u : (1u << level_i);
    }
    else if ((quadratureMethod == olb::uq::Quadrature::QuadratureMethod::GaussQuadrature)) {
      // classic n-point Gauss rule needs order = n-1
      m = static_cast<std::size_t>(level_i + 1);
    }
  auto quadrature = olb::uq::Quadrature::makeQuadrature<T>(
      polynomialBases[dim], m, quadratureMethod
  );
  if (!quadrature) {
    throw std::runtime_error("Failed to create quadrature for level " + std::to_string(level_i) +
                             " and dimension " + std::to_string(dim));
  }
  points  = quadrature->getPoints();
  weights = quadrature->getWeights();
}

// ---------------------------------------------------------------------
//  Smolyak combination coefficient
//      c(i) = (-1)^(L-|i|1) * C(d-1, L-|i|1)
//      with the convention  C(n,k)=0  for k<0  or  k>n
// ---------------------------------------------------------------------
template<typename T>
T smolyakCoeff(std::size_t L,
               std::size_t d,
               const std::vector<int>& idx)
{
    const int k = std::accumulate(idx.begin(), idx.end(), 0); // |i|₁
    const int q = static_cast<int>(L) - k;                    // L-|i|₁

    /* outside the range of the binomial coefficient -> zero                */
    if (q < 0 || q > static_cast<int>(d) - 1) return T{0};

    /* (-1)^(q)  */
    const T sign = (q & 1) ? T(-1) : T(1);

    /* binomial(d-1, q)  – any integer implementation is fine;
       std::size_t is enough for the small numbers occurring here          */
    auto binom = [](int n, int r) -> std::size_t {
        if (r < 0 || r > n) return 0;
        std::size_t res = 1;
        for (int i = 1; i <= r; ++i)
            res = res * (n - r + i) / i;
        return res;
    };

    return sign * static_cast<T>( binom(static_cast<int>(d) - 1, q) );
}



// Generate sparse grid
template <typename T>
void generateSparseGrid(
    std::size_t level,  // Smolyak level L
    std::size_t randomNumberDimension,
    const std::vector<std::shared_ptr<olb::uq::Polynomials::PolynomialBasis<T>>>& polynomialBases,
    olb::uq::Quadrature::QuadratureMethod quadratureMethod,
    std::vector<std::vector<T>>& outputPoints,  // Transposed output
    std::vector<T>& outputWeightsMultiplied)
{
  // Implement sparse grid generation logic here
  // -----------------------------------------------------------------
  // 1. collect raw nodes and weights (duplicates possible)
  // -----------------------------------------------------------------
  T tol = 1e-14;
  std::vector<std::vector<T>> rawPts;
  std::vector<T>              rawWts;

  std::vector<int> multiIndex(randomNumberDimension, 0);

  // recursion over all multi-indices 0 ≤ i_k ≤ level
  std::function<void(std::size_t)> recurse;
  recurse = [&](std::size_t axis) {
    if (axis == randomNumberDimension) {
      if (std::accumulate(multiIndex.begin(), multiIndex.end(), std::size_t{0}) > level)
        return;

      T coeff = smolyakCoeff<T>(level, randomNumberDimension, multiIndex);
      if (coeff == T{}) return;

      // ---- tensor product of the 1-D rules --------------------
      std::vector<std::vector<T>> axesPts(randomNumberDimension);
      std::vector<std::vector<T>> axesWts(randomNumberDimension);
      for (std::size_t k = 0; k < randomNumberDimension; ++k) {
        get1DRule<T>(
          multiIndex[k],
          static_cast<int>(k),
          polynomialBases,
          quadratureMethod,
          axesPts[k],
          axesWts[k]
        );
      }

      std::vector<T> point(randomNumberDimension, T{});
      std::function<void(std::size_t, T)> tensorRec;
      tensorRec = [&](std::size_t ax, T wAcc) {
        if (ax == randomNumberDimension) {
          rawPts.push_back(point);
          rawWts.push_back(wAcc * coeff);
          return;
        }
        for (std::size_t j = 0; j < axesPts[ax].size(); ++j) {
          point[ax] = axesPts[ax][j];
          tensorRec(ax + 1, wAcc * axesWts[ax][j]);
        }
      };
      tensorRec(0, T{1});
      return;
    }

    for (std::size_t k = 0; k <= level; ++k) {
      multiIndex[axis] = k;
      recurse(axis + 1);
    }
  };
  recurse(0);

  // -----------------------------------------------------------------
  // 2. merge duplicate nodes coming from nested rules
  // -----------------------------------------------------------------
constexpr double SCALE = 1e10;          // keep 10 decimal places

using IntKey = std::vector<long long>;  // integerised coordinate

struct KeyLess
{
    bool operator()(const IntKey& a, const IntKey& b) const
    {
        return std::lexicographical_compare(a.begin(), a.end(),
                                            b.begin(), b.end());
    }
};

std::map<IntKey, T, KeyLess> merged;

auto makeKey = [](const std::vector<T>& v)
{
    IntKey key(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
        key[i] = static_cast<long long>(std::llround(v[i] * SCALE));
    return key;
};

// accumulate
for (std::size_t n = 0; n < rawPts.size(); ++n)
    merged[ makeKey(rawPts[n]) ] += rawWts[n];

  // -----------------------------------------------------------------
  // 3. filter out (numerically) zero-weight nodes
  // -----------------------------------------------------------------
  std::vector<std::vector<T>> uniqPts;
  std::vector<T>              uniqWts;

  for (const auto& [ikey, wt] : merged)
  {
      if (std::abs(wt) > tol) {
          std::vector<T> pt(ikey.size());
          for (std::size_t i = 0; i < ikey.size(); ++i)
              pt[i] = static_cast<T>(ikey[i]) / SCALE;

          uniqPts.push_back(std::move(pt));
          uniqWts.push_back(wt);
      }
  }
  outputPoints            = MatrixOperations::transposeMatrix(std::move(uniqPts));
  outputWeightsMultiplied = std::move(uniqWts);

}

} // namespace smolyak
} // namespace uq
} // namespace olb

#endif // SMOLYAK_H
