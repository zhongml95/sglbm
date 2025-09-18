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

#include "polynomials/polynomial.h"
#include "quadratureBase.h"

#include "quadratureRules/gaussHermiteQuadrature.h"
#include "quadratureRules/gaussLegendreQuadrature.h"
#include "quadratureRules/clenshawCurtisQuadrature.h"
#include "quadratureRules/genzKeisterQuadrature.h"

namespace olb {

namespace uq {

namespace Quadrature {

template <typename T>
std::unique_ptr<QuadratureBase<T>>
makeQuadrature(std::size_t order,
               QuadratureMethod method,
               QRMethod qrMethod = QRMethod::WilkinsonShiftQR)
{
  switch (method) {
    case QuadratureMethod::GaussLegendre:
      return std::make_unique<GaussLegendreQuadrature<T>>(order, qrMethod);

    case QuadratureMethod::GaussHermite:
      return std::make_unique<GaussHermiteQuadrature<T>>(order, qrMethod);

    case QuadratureMethod::ClenshawCurtis:
      return std::make_unique<ClenshawCurtisQuadrature<T>>(order);

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
  }

  // Should be unreachable if all enum cases handled
  throw std::runtime_error("Unknown QuadratureMethod in makeQuadrature");
}


} // namespace Quadrature

} // namespace uq

} // namespace olb

#endif // QUADRATURE_FACTORY_H
