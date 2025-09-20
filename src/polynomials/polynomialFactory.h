#pragma once
#ifndef POLYNOMIAL_FACTORY_H
#define POLYNOMIAL_FACTORY_H

#include <memory>
#include <vector>

#include "../distribution.h"
#include "polynomialBasis.h"
#include "hermiteBasis.h"
#include "legendreBasis.h"

namespace olb {
namespace uq {
namespace polynomials {

// Single-dist convenience
template <typename T>
std::shared_ptr<polynomialBasis<T>>
createPolynomialBasis(olb::uq::DistributionType dist, std::size_t order)
{
  switch (dist) {
    case olb::uq::DistributionType::Uniform:
      return std::make_shared<LegendreBasis<T>>(order);
    case olb::uq::DistributionType::Normal:
      return std::make_shared<HermiteBasis<T>>(order);
    default:
      // fallback
      return std::make_shared<LegendreBasis<T>>(order);
  }
}

// Vector version
template <typename T>
std::vector<std::shared_ptr<polynomialBasis<T>>>
createPolynomialBases(const std::vector<olb::uq::DistributionType>& dists,
                      std::size_t order)
{
  std::vector<std::shared_ptr<polynomialBasis<T>>> bases;
  bases.reserve(dists.size());

  for (auto dist : dists) {
    switch (dist) {
      case olb::uq::DistributionType::Uniform:
        bases.push_back(std::make_shared<LegendreBasis<T>>(order));
        break;
      case olb::uq::DistributionType::Normal:
        bases.push_back(std::make_shared<HermiteBasis<T>>(order));
        break;
      default:
        bases.push_back(std::make_shared<LegendreBasis<T>>(order)); // fallback
        break;
    }
  }
  return bases;
}

template <typename T>
std::vector<std::shared_ptr<polynomialBasis<T>>>
createPolynomialBases(const std::vector<Distribution<T>>& dists,
                      std::size_t order)
{
  std::vector<std::shared_ptr<polynomialBasis<T>>> bases;
  bases.reserve(dists.size());

  for (const auto& dist : dists) {
    switch (dist.getType()) {   // assume Distribution<T> has getType()
      case olb::uq::DistributionType::Uniform:
        bases.push_back(std::make_shared<LegendreBasis<T>>(order));
        break;

      case olb::uq::DistributionType::Normal:
        bases.push_back(std::make_shared<HermiteBasis<T>>(order));
        break;

      // You can extend here for Beta, Exponential, etc.
      default:
        bases.push_back(std::make_shared<LegendreBasis<T>>(order)); // fallback
        break;
    }
  }

  return bases;
}


} // namespace polynomials
} // namespace uq
} // namespace olb

#endif // POLYNOMIAL_FACTORY_H
