// uq/stochasticCollocationGrid.h
#pragma once
#include <vector>
#include <memory>
#include <cstddef>
#include "quadrature/quadrature.h"
#include "polynomials/polynomial.h"

namespace olb { namespace uq {

template <typename T>
struct StochasticCollocationGrid {
  // Always dimension-major:
  // - Tensor grid:   points[i] has size nqPerDim
  // - Sparse grid:   points[i] has size totalNq   (one coord per sample along dim i)
  std::vector<std::vector<T>> points;
  std::vector<std::vector<T>> pointsTensor; // Tensor product of points (tensor only)

  // Per-dimension weights (tensor only). For sparse, leave empty.
  std::vector<std::vector<T>> weights;

  // One combined weight per sample (sparse always; tensor optional precompute).
  std::vector<T>              weightsMultiplied;

  // Tensor-only: map from flat sample index -> tuple of per-dim indices.
  // For sparse, usually empty.
  std::vector<std::vector<std::size_t>> indexList;

  std::size_t nqPerDim = 0;  // tensor only (0 for sparse)
  std::size_t totalNq  = 0;  // number of samples (both tensor & sparse)
  bool isSparse = false;
};

// Build a tensor (dense) grid: matches your old dense branch.
template <typename T>
StochasticCollocationGrid<T> createTensorCollocationGrid(
  const std::vector<std::shared_ptr<polynomials::polynomialBasis<T>>>& bases,
  std::size_t nq,
  Quadrature::QuadratureMethod method);

// Build a Smolyak sparse grid: points[i] is length totalNq; weightsMultiplied size totalNq.
template <typename T>
StochasticCollocationGrid<T> createSparseCollocationGrid(
  const std::vector<std::shared_ptr<polynomials::polynomialBasis<T>>>& bases,
  std::size_t level,
  Quadrature::QuadratureMethod method);

// Convenience: choose tensor vs sparse.
template <typename T>
StochasticCollocationGrid<T> createCollocationGrid(
  const std::vector<std::shared_ptr<polynomials::polynomialBasis<T>>>& bases,
  std::size_t nqOrLevel,
  Quadrature::QuadratureMethod method,
  bool sparse);

}} // namespace olb::uq
