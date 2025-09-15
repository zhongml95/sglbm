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
#ifndef CLENSHAW_CURTIS_QUADRATURE_H
#define CLENSHAW_CURTIS_QUADRATURE_H

#include <stdexcept> // For std::runtime_error
#include <memory>    // For std::unique_ptr
#include "quadrature/quadratureBase.h"
#include "quadrature/quadratureRules/fft.h"


namespace olb {

namespace uq {

namespace Quadrature {

template <typename T>
class ClenshawCurtisQuadrature : public QuadratureBase<T> {
public:
    ClenshawCurtisQuadrature(std::size_t order) {
        if (order == 0) {
            points = {T(0.0)};
            weights = {T(1.0)};
            return;
        } else if (order == 1) {
            points = {T(-1.0), T(1.0)};
            weights = {T(0.5), T(0.5)};
            return;
        }

        std::vector<T> theta(order + 1);
        for (std::size_t k = 0; k <= order; ++k) {
            theta[k] = (order - k) * M_PI / order;
        }
        points.resize(order + 1);
        for (std::size_t k = 0; k <= order; ++k) {
            points[k] = std::cos(theta[k]);
        }


        // Match Python: steps = np.arange(1, order, 2)
        std::vector<T> steps;
        for (std::size_t k = 1; k < order; k += 2) {
            steps.push_back(static_cast<T>(k));
        }

        std::size_t length = steps.size();
        std::size_t remains = order - length;

        std::vector<T> beta;

        // Step 1: 2.0 / (k * (k - 2)) for each step
        for (std::size_t i = 0; i < steps.size(); ++i) {
            T k = steps[i];
            beta.push_back(2.0 / (k * (k - 2.0)));
        }

        // Step 2: Append [1.0 / steps.back()]
        beta.push_back(1.0 / steps.back());

        // Step 3: Append remains zeros
        for (std::size_t i = 0; i < remains; ++i) {
            beta.push_back(0.0);
        }

        std::size_t N = beta.size();
        std::vector<T> new_beta(N - 1);

        for (std::size_t i = 0; i < N - 1; ++i) {
            new_beta[i] = -beta[i] - beta[N - 1 - i];
        }

        beta = new_beta;

        std::vector<T> gamma(order, -1.0);
        gamma[length] += static_cast<T>(order);
        gamma[remains] += static_cast<T>(order);
        T denom = static_cast<T>(order * order - 1 + (order % 2));
        for (auto& g : gamma) {
            g /= denom;
        }

        assert(beta.size() == gamma.size());

        std::vector<T> combined(order);
        for (std::size_t i = 0; i < order; ++i) {
          combined[i] = beta[i] + gamma[i];
        }

      std::vector<std::complex<T>> combined_complex(combined.size(), 0.0);
      for (std::size_t i = 0; i < combined.size(); ++i) {
          combined_complex[i] = std::complex<T>(combined[i], 0.0);
      }

      olb::uq::fft::ifft_dispatch<T>(combined_complex);

      std::vector<T> half_weights;
      for (std::size_t i = 0; i < (combined_complex.size()/2 + 1); ++i) {
          half_weights.push_back(std::real(combined_complex[i]));
      }

      // Step 2: Mirror and divide by 2
      weights = half_weights;

      int start = static_cast<int>(half_weights.size()) - 2 + (order % 2);
      for (int i = start; i >= 0; --i) {
          weights.push_back(half_weights[i]);
      }

      for (std::size_t i = 0; i < weights.size(); ++i) {
        weights[i] /= 2.0;
      }
    }

    const std::vector<T>& getPoints() const override { return points; }
    const std::vector<T>& getWeights() const override { return weights; }

private:
    std::vector<T> points;
    std::vector<T> weights;
};



} // namespace Quadrature

} // namespace uq
} // namespace olb


#endif // CLENSHAW_CURTIS_QUADRATURE_H