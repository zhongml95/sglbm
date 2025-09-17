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

// quadrature_base.h
#ifndef QUADRATURE_BASE_H
#define QUADRATURE_BASE_H

#include "filesIO.h"

namespace olb {

namespace uq {

namespace Quadrature {

enum class QRMethod { HouseholderQR, WilkinsonShiftQR };
enum class QuadratureMethod { GaussQuadrature, GenzKeister16, GenzKeister18, GenzKeister22, GenzKeister24, ClenshawCurtis };

template <typename T>
class QuadratureBase {
public:
  virtual ~QuadratureBase()                        = default;
  virtual const std::vector<T>& getPoints() const  = 0;
  virtual const std::vector<T>& getWeights() const = 0;
};

// Helper: map string to QuadratureMethod
inline olb::uq::Quadrature::QuadratureMethod
toQuadratureMethod(const std::string& rule, unsigned nq, unsigned order) {
    using QM = olb::uq::Quadrature::QuadratureMethod;
    std::string s = rule;
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (s == "gauss" || s == "gaussquadrature" || s == "legendre")
        return QM::GaussQuadrature;
    if (s == "clenshaw" || s == "clenshawcurtis" || s == "cc")
        return QM::ClenshawCurtis;
    if (s == "gk16" || s == "genzkeister16")
        return QM::GenzKeister16;
    if (s == "gk18" || s == "genzkeister18")
        return QM::GenzKeister18;
    if (s == "gk22" || s == "genzkeister22")
        return QM::GenzKeister22;
    if (s == "gk24" || s == "genzkeister24")
        return QM::GenzKeister24;
    if (s == "gk" || s == "genzkeister") {
        // choose a level based on nq or order
        if ((nq ? nq : order) <= 3) return QM::GenzKeister16;
        if ((nq ? nq : order) <= 5) return QM::GenzKeister18;
        if ((nq ? nq : order) <= 7) return QM::GenzKeister22;
        return QM::GenzKeister24;
    }

    // default
    return QM::GaussQuadrature;
}


} // namespace Quadrature

} // namespace uq

} // namespace olb

#endif // QUADRATURE_BASE_H
