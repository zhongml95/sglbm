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

#ifndef GENZ_KEISTER_QUADRATURE_H
#define GENZ_KEISTER_QUADRATURE_H


#include "quadrature/quadratureRules/genzKeisterLibrary.h"
#include "quadrature/quadratureRules/fft.h"
#include "quadrature/quadratureBase.h"

namespace olb {

namespace uq {

namespace Quadrature {

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


} // namespace Quadrature

} // namespace uq

} // namespace olb

#endif // GENZ_KEISTER_QUADRATURE_H