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

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

namespace olb {

namespace uq {
  
// Enumeration to specify the distribution type
enum class DistributionType {
  Uniform,
  Normal,
  // MultivariateNormal
};

  
// --- helpers: string <-> enum ---
inline DistributionType distributionTypeFromString(std::string s) {
    // case-insensitive, with a few aliases
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    if (s == "uniform") return DistributionType::Uniform;
    if (s == "normal" || s == "gaussian" || s == "gauss") return DistributionType::Normal;
    throw std::invalid_argument("Unknown distribution type: " + s);
}

inline std::string distributionTypeToString(DistributionType t) {
    switch (t) {
        case DistributionType::Uniform: return "Uniform";
        case DistributionType::Normal:  return "Normal";
        default:                        return "Unknown";
    }
}

inline std::vector<DistributionType>
distributionTypesFromStrings(const std::vector<std::string>& strs) {
    std::vector<DistributionType> types;
    types.reserve(strs.size());
    for (auto& s : strs) {
        types.push_back(distributionTypeFromString(s));
    }
    return types;
}



// Struct to hold distribution information
template <typename T>
struct Distribution {
  DistributionType type;

  // Parameters for univariate distributions
  T param1 = 0.0; // For Uniform: lower bound, for Normal: mean
  T param2 = 1.0; // For Uniform: upper bound, for Normal: standard deviation

  // Constructor for univariate distributions
  Distribution(DistributionType type, T param1 = 0.0, T param2 = 1.0)
      : type(type)
      , param1(param1)
      , param2(param2)
  {}

  DistributionType getType() const noexcept {
    return type;
  }
};

// Factory functions for cleaner syntax
template <typename T>
inline Distribution<T> uniform(T min, T max)
{
  return Distribution(DistributionType::Uniform, min, max);
}

template <typename T>
inline Distribution<T> normal(T mean, T stddev)
{
  return Distribution(DistributionType::Normal, mean, stddev);
}

// Joint distribution as std::vector<Distribution>
template <typename T>
inline std::vector<Distribution<T>> joint(const std::vector<Distribution<T>>& dists)
{
  return dists;
}

// Affine transformation function
template <typename T>
T affine(T x, const Distribution<T>& dist)
{
  switch (dist.type) {
  case DistributionType::Uniform:
    return 0.5 * (dist.param2 - dist.param1) * x + 0.5 * (dist.param1 + dist.param2);
  case DistributionType::Normal:
    return dist.param1 + dist.param2 * x;
  default:
    throw std::runtime_error("Unsupported distribution type for affine transformation.");
  }
}


// --- single factories ---
template <typename T>
inline Distribution<T> createDistribution(DistributionType type, T p1, T p2) {
    switch (type) {
        case DistributionType::Uniform:
            // p1=min, p2=max
            if (!(p1 <= p2)) {
                throw std::invalid_argument("Uniform: require min <= max.");
            }
            return uniform<T>(p1, p2);

        case DistributionType::Normal:
            // p1=mean, p2=stddev
            if (!(p2 > T(0))) {
                throw std::invalid_argument("Normal: stddev must be > 0.");
            }
            return normal<T>(p1, p2);
    }
    throw std::invalid_argument("Unsupported DistributionType.");
}

template <typename T>
inline Distribution<T> createDistribution(const std::string& typeStr, T p1, T p2) {
    return createDistribution<T>(distributionTypeFromString(typeStr), p1, p2);
}

// --- batch factory ---
template <typename T>
inline std::vector<Distribution<T>>
createDistributions(const std::vector<std::string>& types,
                    const std::vector<T>& p1,
                    const std::vector<T>& p2)
{
    if (types.size() != p1.size() || types.size() != p2.size()) {
        throw std::invalid_argument("createDistributions: size mismatch among types, p1, p2.");
    }
    std::vector<Distribution<T>> out;
    out.reserve(types.size());
    for (size_t i = 0; i < types.size(); ++i) {
        out.push_back(createDistribution<T>(types[i], p1[i], p2[i]));
    }
    return out;
}


} // namespace uq

} // namespace olb

#endif // DISTRIBUTION_H
