// quadrature_base.h
#ifndef QUADRATURE_BASE_H
#define QUADRATURE_BASE_H

#include "utils.h"

namespace Quadrature {

enum class QuadratureMethod {
    HouseholderQR,
    GSL
};

class QuadratureBase {
public:
    virtual ~QuadratureBase() = default;
    virtual const std::vector<double>& getPoints() const = 0;
    virtual const std::vector<double>& getWeights() const = 0;
};

} // namespace Quadrature

#endif // QUADRATURE_BASE_H
