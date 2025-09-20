// uq/stochasticCollocationGrid.hh
#pragma once
#include "stochasticCollocationGrid.h"
#include <cmath>
#include <cassert>

namespace olb { namespace uq {

namespace detail {

// odometer over nq^dim, fastest varying dim = 0
inline std::vector<std::vector<std::size_t>>
buildIndexList(std::size_t dim, std::size_t nqPerDim) {
  const std::size_t total = static_cast<std::size_t>(std::pow(nqPerDim, dim));
  std::vector<std::vector<std::size_t>> idxList;
  idxList.reserve(total);
  std::vector<std::size_t> idx(dim, 0);

  for (std::size_t k = 0; k < total; ++k) {
    idxList.push_back(idx);
    for (std::size_t j = 0; j < dim; ++j) {
      if (++idx[j] < nqPerDim) break;
      idx[j] = 0;
    }
  }
  return idxList;
}

} // namespace detail

template <typename T>
StochasticCollocationGrid<T> createTensorCollocationGrid(
  std::size_t dim,
  std::size_t nq,
  const std::vector<olb::uq::Quadrature::QuadratureMethod>& quadratureMethods)
{
  if (dim == 0) {
    throw std::invalid_argument("createTensorCollocationGrid: dim must be >= 1");
  }
  if (nq == 0) {
    throw std::invalid_argument("createTensorCollocationGrid: nq must be >= 1");
  }

  StochasticCollocationGrid<T> g;
  g.isSparse = false;

  // 1) Per-dimension 1D rules (nodes/weights length = nq)
  g.points.resize(dim);
  g.weights.resize(dim);

  for (std::size_t i = 0; i < dim; ++i) {
    auto q = Quadrature::makeQuadrature<T>(nq, quadratureMethods[i]);
    g.points[i]  = q->getPoints();   // length nq
    g.weights[i] = q->getWeights();  // length nq
  }

  // 2) Sizes
  g.nqPerDim = g.points[0].size();       // == nq
  // integer pow to avoid fp round-trip
  auto ipow = [](std::size_t a, std::size_t b) {
    std::size_t r = 1;
    while (b--) r *= a;
    return r;
  };
  g.totalNq = ipow(g.nqPerDim, dim);

  // 3) Build index list (tensor lexicographic order)
  //    Convention: dimension 0 is least significant / fastest varying.
  g.indexList.assign(g.totalNq, std::vector<std::size_t>(dim, 0));
  for (std::size_t k = 0; k < g.totalNq; ++k) {
    std::size_t idx = k;
    // fill from dim-1 down to 0 (last dimension fastest)
    for (std::size_t i = dim; i-- > 0;) {
      g.indexList[k][i] = idx % g.nqPerDim;
      idx /= g.nqPerDim;
    }
  }

  // 4) Precompute combined weights
  g.pointsTensor.assign(g.totalNq, std::vector<T>(dim, T(0)));
  g.weightsMultiplied.assign(g.totalNq, T(1));
  T wSum = 0;
  for (std::size_t k = 0; k < g.totalNq; ++k) {
    T w = T(1);
    for (std::size_t d = 0; d < dim; ++d) {
      const std::size_t j = g.indexList[k][d];
      w *= g.weights[d][j];
      g.pointsTensor[k][d] = g.points[d][j];    // assign coordinate for sample k, dim d
    }
    g.weightsMultiplied[k] = w;
    wSum += w;
  }

  for (auto& wt : g.weightsMultiplied) {
    wt /= wSum; // normalize
  }
  g.points = g.pointsTensor; // assign full tensor points
  return g;
}

template <typename T>
StochasticCollocationGrid<T> createSparseCollocationGrid( std::size_t dim,
  std::size_t level,
  const std::vector<olb::uq::Quadrature::QuadratureMethod>& quadratureMethods)
{
  StochasticCollocationGrid<T> g;
  g.isSparse = true;
  g.nqPerDim = 0;              // not meaningful for sparse
  // g.points.resize(dim);        // points[d].size() will be totalNq after generation

  // Fill dimension-major flattened points and combined weights
  Quadrature::generateSparseGrid<T>(
    level, dim, quadratureMethods, g.points, g.weightsMultiplied
  );

  // totalNq is the number of flattened samples
  assert(!g.points.empty());
  g.totalNq = g.points.size();
  // indexList/weights (per-dim) are not used for sparse
  return g;
}

template <typename T>
StochasticCollocationGrid<T> createCollocationGrid(
  std::size_t nqOrLevel,
  const std::vector<olb::uq::Quadrature::QuadratureMethod>& quadratureMethods,
  bool sparse)
{
  const std::size_t dim = quadratureMethods.size();
  return sparse
    ? createSparseCollocationGrid<T>(dim, nqOrLevel, quadratureMethods)
    : createTensorCollocationGrid<T>(dim, nqOrLevel, quadratureMethods);
}

template <typename T>
StochasticCollocationGrid<T> createCollocationGrid(
  std::size_t nqOrLevel,
  std::size_t polyOrder,
  const std::vector<std::string>& rules,
  bool sparse)
{
  std::vector<olb::uq::Quadrature::QuadratureMethod> methods;
  for (const auto& rule : rules) {
    methods.push_back(Quadrature::toQuadratureMethod(rule, nqOrLevel, polyOrder));
  }
  return createCollocationGrid<T>(nqOrLevel, methods, sparse);
}

}} // namespace olb::uq
