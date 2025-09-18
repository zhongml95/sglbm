#pragma once
#include <vector>
#include <cstddef>
#include "polynomialBasis.h"

namespace olb { namespace uq { namespace polynomials {

// Legendre P_n on [-1,1], L2-orthogonal with weight 1.
template <typename T>
class LegendreBasis : public polynomialBasis<T> {
public:
  explicit LegendreBasis(std::size_t order)
  : polynomialBasis<T>(order) {
    // compute and SAVE the whole coefficients matrix once:
    this->buildAll();
  }

  std::string name() const override { return "Legendre"; }
  BasisKind   kind() const noexcept override { return BasisKind::Legendre; }

  std::vector<T> computeCoefficients(std::size_t n) const override {
    if (n == 0) return { T(1) };
    if (n == 1) return { T(0), T(1) };

    // Recurrence: n P_n = (2n-1) x P_{n-1} - (n-1) P_{n-2}
    std::vector<T> Pnm2{ T(1) };          // P_0
    std::vector<T> Pnm1{ T(0), T(1) };    // P_1
    std::vector<T> Pn;

    for (std::size_t k = 2; k <= n; ++k) {
      Pn.assign(k + 1, T(0));
      const T a = T(2*k - 1) / T(k);
      const T b = T(k - 1)   / T(k);

      // a * x * P_{k-1} (shift by 1)
      for (std::size_t i = 0; i < Pnm1.size(); ++i)
        Pn[i + 1] += a * Pnm1[i];

      // - b * P_{k-2}
      for (std::size_t i = 0; i < Pnm2.size(); ++i)
        Pn[i] -= b * Pnm2[i];

      Pnm2.swap(Pnm1);
      Pnm1.swap(Pn);
    }
    return Pnm1;
  }
};

}}} // ns
