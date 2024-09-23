#include <cmath>
#include <iostream>
#include "../../src/generalized_polynomial_chaos.h"


int main() {
    // Set parameters
    int order = 4; // Polynomial order
    int nq = 8;    // Number of quadrature points

    // Choose the quadrature method
    Quadrature::QuadratureMethod quadratureMethod = Quadrature::QuadratureMethod::GSL;

    // Initialize polynomial bases
    std::vector<std::shared_ptr<Polynomials::PolynomialBasis>> polynomialBases;
    polynomialBases.push_back(std::make_shared<Polynomials::LegendreBasis>());

    // Create GeneralizedPolynomialChaos object
    GeneralizedPolynomialChaos ops(order, nq, polynomialBases, quadratureMethod);

    // // Get points and weights
    std::vector<std::vector<double>> points, weights;
    ops.getPointsAndWeights(points, weights);

    // // Output the points and weights
    std::cout << "Quadrature Points and Weights:\n";
    for (size_t dim = 0; dim < points.size(); ++dim) {
        std::cout << "Dimension " << dim + 1 << ":\n";
        std::cout << "Points: ";
        for (const auto& p : points[dim]) {
            std::cout << p << " ";
        }
        std::cout << "\nWeights: ";
        for (const auto& w : weights[dim]) {
            std::cout << w << " ";
        }
        std::cout << "\n\n";
    }

    std::cout << "weightsMultiplied: ";
    for (const auto& w : ops.getWeightsMultiplied()) {
        std::cout << w << " ";
    }
    std::cout << "\n\n";

    // Output the t2Product tensor
    std::vector<double> t2Product, t2Product_inv, t3Product;
    ops.getTensors(t2Product, t2Product_inv, t3Product);
    std::cout << "t2Product tensor:\n";
    for (size_t i = 0; i < t2Product.size(); ++i) {
        std::cout << t2Product[i] << " ";
    }
    std::cout << "\n\n";

    // Output phiRand matrix
    std::cout << "phiRan matrix:\n";
    std::vector<double> phiRan;
    ops.getPhiRan(phiRan);
    
    for (size_t i = 0; i < ops.getQuadraturePointsNumber(); ++i) {
        for (size_t j = 0; j < ops.getPolynomialsOrder(); ++j) {
            std::cout << phiRan[i * ops.getPolynomialsOrder() + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n";

    // Output the coefficients of polynomials
    std::vector<std::vector<std::vector<double>>> polynomialCoeffs;
    ops.getCoefficients(polynomialCoeffs);

    std::cout << "Polynomial coefficients:\n";
    for (size_t i = 0; i < polynomialCoeffs.size(); ++i) {
        std::cout << "Dimension " << i + 1 << ":\n";
        for (size_t j = 0; j < polynomialCoeffs[i].size(); ++j) {
            std::cout << "Order " << j << ": ";
            for (size_t k = 0; k < polynomialCoeffs[i][j].size(); ++k) {
                std::cout << polynomialCoeffs[i][j][k] << " ";
            }
            std::cout << "\n";
        }
    }
    std::cout << "\n";

    // Check the inds
    std::vector<std::vector<int>> inds = ops.getMultiIndices();
    std::cout << "Multi-indices:\n";
    for (size_t i = 0; i < inds.size(); ++i) {
        for (size_t j = 0; j < inds[i].size(); ++j) {
            std::cout << inds[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    return 0;
}
