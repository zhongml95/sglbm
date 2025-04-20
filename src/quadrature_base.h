// quadrature_base.h
#ifndef QUADRATURE_BASE_H
#define QUADRATURE_BASE_H

#include "utils.h"


// namespace olb {

// namespace uq {

namespace Quadrature {

enum class QuadratureMethod {
    HouseholderQR,
    WilkinsonShiftQR
};

template<typename T>
class QuadratureBase {
public:
    virtual ~QuadratureBase() = default;
    virtual const std::vector<T>& getPoints() const = 0;
    virtual const std::vector<T>& getWeights() const = 0;
};

} // namespace Quadrature

// } // namespace uq

// } // namespace olb

#endif // QUADRATURE_BASE_H
