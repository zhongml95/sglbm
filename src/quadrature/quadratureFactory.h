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

// quadratureFactory.h
#ifndef QUADRATURE_FACTORY_H
#define QUADRATURE_FACTORY_H

#include "polynomial.h"
#include "quadratureBase.h"

#include "quadratureRules/gaussQuadrature.h"
#include "quadratureRules/clenshawCurtisQuadrature.h"
#include "quadratureRules/genzKeisterQuadrature.h"

namespace olb {

namespace uq {

namespace Quadrature {


template <typename>
struct always_false : std::false_type {};

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

#endif // QUADRATURE_FACTORY_H
