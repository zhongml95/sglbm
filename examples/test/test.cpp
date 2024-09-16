#include <cmath>
// #include "../../src/sglbm.h"
// #include "../../src/postprocessing_sglbm.h"
#include "../../src/generalized_polynomial_chaos.h"

#define USE_GSL

int main() {
    int order = 3;
    int nq = 7;
    std::vector<double> parameters1(1, -1.0);
    std::vector<double> parameters2(1,  1.0);
    std::vector<int> parameterType(1, 0);
    Quadrature::QuadratureMethod points_weights_method = Quadrature::QuadratureMethod::HouseholderQR;

    GeneralizedPolynomialChaos ops( order, nq, parameters1, parameters2, parameterType, points_weights_method );

}