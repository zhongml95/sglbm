#include <cmath>
#include "../../src/sglbm.h"
#include "../../src/postprocessing_sglbm.h"

int main() {
    int order = 3;
    int nq = 7;
    std::vector<double> parameters1(1, -1.0);
    std::vector<double> parameters2(1,  1.0);
    std::vector<int> parameterType(1, 0);
    int points_weights_method = 1;

    polynomials ops( nq, order, parameters1, parameters2, parameterType, points_weights_method );

}