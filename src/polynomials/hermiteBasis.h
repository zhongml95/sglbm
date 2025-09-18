#pragma once
#include <vector>
#include <cstddef>
#include "polynomialBasis.h"


namespace olb { namespace uq { namespace polynomials {

// Probabilist's Hermite He_n (H0=1, H1=x, He_{n} = x He_{n-1} - (n-1) He_{n-2})
template <typename T>
class HermiteBasis : public polynomialBasis<T> {
public:
  explicit HermiteBasis(std::size_t order)
  : polynomialBasis<T>(order) {
    // compute and SAVE the whole coefficients matrix once:
    this->buildAll();
  }

  std::string name() const override { return "Hermite"; }
  BasisKind   kind() const noexcept override { return BasisKind::Hermite; }

  std::vector<T> computeCoefficients(std::size_t n) const override {
    if (n == 0) return { T(1) };
    if (n == 1) return { T(0), T(1) };

    std::vector<T> Hnm2{ T(1) };         // He_0
    std::vector<T> Hnm1{ T(0), T(1) };   // He_1
    std::vector<T> Hn;

    for (std::size_t k = 2; k <= n; ++k) {
      Hn.assign(k + 1, T(0));
      // x * He_{k-1}
      for (std::size_t i = 0; i < Hnm1.size(); ++i)
        Hn[i + 1] += Hnm1[i];
      // - (k-1) * He_{k-2}
      for (std::size_t i = 0; i < Hnm2.size(); ++i)
        Hn[i] -= T(k - 1) * Hnm2[i];

      Hnm2.swap(Hnm1);
      Hnm1.swap(Hn);
    }
    return Hnm1;
  }
};

}}} // ns
