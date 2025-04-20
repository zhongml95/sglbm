// generalized_polynomial_chaos.hh
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_HH
#define GENERALIZED_POLYNOMIAL_CHAOS_HH

#include "generalized_polynomial_chaos.h"

#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <functional>
#include "matrix_operation.h"

// Include the polynomial basis and quadrature headers
// #include "quadrature.h"

// Namespace aliases
template<typename T>
using LegendreBasis = Polynomials::LegendreBasis<T>;
template<typename T>
using HermiteBasis  = Polynomials::HermiteBasis<T>;
// using Quadrature = Quadrature::Quadrature;

// namespace olb {

// namespace uq {

// Constructor
template<typename T>
GeneralizedPolynomialChaos<T>::GeneralizedPolynomialChaos(size_t _order,
                                                       size_t _nq,
                                                       const std::vector<Distribution<T>>& _distributions,
                                                       Quadrature::QuadratureMethod _quadratureMethod)
    : pointsWeightsMethod(0),
      No(0),
      nq(_nq),
      totalNq(0),
      order(_order),
      randomNumberDimension(0),
      quadratureMethod(_quadratureMethod),
      distributions(_distributions) // Initialize distributions
{
    // Map distributions to polynomial bases
    polynomialBases.clear();
    for (const auto& dist : distributions) {
        polynomialBases.push_back(createPolynomialBasis(dist));
    }

    // Set randomNumberDimension to the size of distributions
    randomNumberDimension = distributions.size();

    // Calculate multi-indices
    calculateMultiIndices(randomNumberDimension, order, inds);
    No = inds.size();

    // Initialize polynomial coefficients, quadratures, and matrices
    initializeQuadratures();
    initializeMatrices();

    // Evaluate polynomials at quadrature points
    evaluatePhiRan();

    // Compute tensors
    computeTensors();
}


// Initialize quadratures
template<typename T>
void GeneralizedPolynomialChaos<T>::initializeQuadratures() {
    points.resize(randomNumberDimension);
    weights.resize(randomNumberDimension);

    totalNq = std::pow(nq, randomNumberDimension);

    for (size_t i = 0; i < randomNumberDimension; ++i) {
        auto quadrature = polynomialBases[i]->getQuadrature(nq, quadratureMethod);
        points[i] = quadrature->getPoints();
        weights[i] = quadrature->getWeights();
    }
}



// Initialize matrices
template<typename T>
void GeneralizedPolynomialChaos<T>::initializeMatrices() {
    phiRan.resize(totalNq * No, 0.0);
    phiRan_T.resize(totalNq * No, 0.0);
    t2Product.resize(No, 0.0);
    t2Product_inv.resize(No, 0.0);
    t3Product.resize(No * No * No, 0.0);

    // Generate pointsWeightsIndexList
    pointsWeightsIndexList.resize(totalNq, std::vector<size_t>(randomNumberDimension));
    for (size_t i = 0; i < totalNq; ++i) {
        pointsWeightsIndexList[i] = findIndex(i, randomNumberDimension, nq);
    }


    // Compute weightsMultiplied and pointsTensor
    weightsMultiplied.resize(totalNq, 1.0);
    pointsTensor.resize(totalNq, std::vector<T>(randomNumberDimension));
    for (size_t k = 0; k < totalNq; ++k) {
        pointsTensor[k].resize(randomNumberDimension);
        for (size_t dim = 0; dim < randomNumberDimension; ++dim) {
            size_t idx = pointsWeightsIndexList[k][dim];
            weightsMultiplied[k] *= weights[dim][idx];
            pointsTensor[k][dim] = points[dim][idx];
        }
    }
}

// Initialize polynomial coefficients
template<typename T>
void GeneralizedPolynomialChaos<T>::initializePolynomialCoefficients() {
    // std::cout << "Initializing polynomial coefficients..." << std::endl;
    coefficients.resize(randomNumberDimension);
    for (size_t phi_i = 0; phi_i < randomNumberDimension; ++phi_i) {
        auto basis = std::static_pointer_cast<LegendreBasis>(polynomialBases[phi_i]);

        coefficients[phi_i].resize(No);
        for (size_t i = 0; i < No; ++i) {
            coefficients[phi_i][i] = basis->computeCoefficients(i);
        }
    }

}

// Evaluate n_order polynomial at point k
template<typename T>
T GeneralizedPolynomialChaos<T>::evaluate(size_t n_order, size_t k) {
    T result = 1.0;
    for (size_t i = 0; i < randomNumberDimension; ++i) {
        result *= evaluate(inds[n_order][i], k, i);
    }
    return result;
}

// Evaluate n_order polynomial at point k and dimension phi_i
template<typename T>
T GeneralizedPolynomialChaos<T>::evaluate(size_t n_order, size_t k, size_t phi_i) {
    T x = points[phi_i][pointsWeightsIndexList[k][phi_i]];
    return evaluate(n_order, x, phi_i);
}

// Evaluate n_order polynomial at given multi-index
template<typename T>
T GeneralizedPolynomialChaos<T>::evaluate(size_t n_order, const std::vector<size_t>& idx) {
    T result = 1.0;
    for (size_t i = 0; i < randomNumberDimension; ++i) {
        result *= evaluate(inds[n_order][i], points[i][idx[i]], i);
    }
    return result;
}

// Evaluate polynomial basis at given order, point x, and dimension phi_i
template<typename T>
T GeneralizedPolynomialChaos<T>::evaluate(size_t n_order, T x, size_t phi_i) {
    if (phi_i < 0 || phi_i >= polynomialBases.size()) {
        throw std::out_of_range("Invalid dimension index phi_i.");
    }
    return polynomialBases[phi_i]->evaluatePolynomial(n_order, x);
}

// Evaluate the polynomial at kth point up to order_max
template<typename T>
T GeneralizedPolynomialChaos<T>::evaluate_polynomial(size_t order_max, size_t k) {
    T sum = 0.0;
    for (size_t i = 0; i <= order_max; ++i) {
        sum += evaluate(i, k);
    }
    return sum;
}

// Evaluate phiRan matrix
template<typename T>
void GeneralizedPolynomialChaos<T>::evaluatePhiRan() {
    // std::cout << "Evaluating phiRan matrix..." << std::endl;
    for (size_t k = 0; k < totalNq; ++k) {
        for (size_t i = 0; i < No; ++i) {
            phiRan[k * No + i] = evaluate(i, pointsWeightsIndexList[k]);
            // std::cout << phiRan[k * No + i] << " ";
            phiRan_T[i * totalNq + k] = phiRan[k * No + i];
        }
        // std::cout << std::endl;
    }
}

// Helper functions
template<typename T>
void GeneralizedPolynomialChaos<T>::calculateMultiIndices(size_t d, size_t n, std::vector<std::vector<size_t>>& indices) {

    std::vector<size_t> index(d, 0);

    std::function<void(size_t, size_t, size_t)> recursiveFunction = [&](size_t pos, size_t sum, size_t maxOrder) {
        if (pos == d - 1) {
            index[pos] = maxOrder - sum;
            indices.push_back(index);
            return;
        }

        for (size_t i = 0; i <= maxOrder - sum; ++i) {
            index[pos] = i;
            recursiveFunction(pos + 1, sum + i, maxOrder);
        }
    };

    for (size_t order = 0; order <= n; ++order) {
        recursiveFunction(0, 0, order);
    }
}

template<typename T>
std::vector<size_t> GeneralizedPolynomialChaos<T>::findIndex(size_t idx, size_t dimension, size_t nq) {
    if (dimension == 1) {
        return {idx};
    }

    // General case for dimension > 1
    std::vector<size_t> index(dimension);
    for (size_t i = dimension; i-- > 0;) {  // Loop from dimension-1 to 0
        index[i] = idx % nq;
        idx /= nq;
    }
    return index;
}

// Compute tensors (t2Product and t3Product)
template<typename T>
void GeneralizedPolynomialChaos<T>::computeTensors() {
    const std::string directoryT2Product = "./t2Product/";
    if (!directoryExists(directoryT2Product)) {
        createDirectory(directoryT2Product);
    }
    const std::string directoryT3Product = "./t3Product/";
    if (!directoryExists(directoryT3Product)) {
        createDirectory(directoryT3Product);
    }
    const std::string t2ProductFile = directoryT2Product + "dims_" + std::to_string(randomNumberDimension) + "_order_" + std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";
    const std::string t3ProductFile = directoryT3Product + "dims_" + std::to_string(randomNumberDimension) + "_order_" + std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";

        for (size_t i = 0; i < No; ++i) {
            for (size_t m = 0; m < totalNq; ++m) {
                t2Product[i] += phiRan[m * No + i] * phiRan[m * No + i] * weightsMultiplied[m];
            }
        }

    for (size_t i = 0; i < No; ++i) {
        t2Product_inv[i] = 1.0 / t2Product[i];
    }


        for (size_t i = 0; i < No; ++i) {
            for (size_t j = 0; j < No; ++j) {
                for (size_t k = 0; k < No; ++k) {
                    T sum = 0.0;
                    for (size_t m = 0; m < totalNq; ++m) {
                        sum += phiRan[m * No + i] * phiRan[m * No + j] * phiRan[m * No + k] * weightsMultiplied[m];
                    }
                    t3Product[i * No * No + j * No + k] = sum;
                }
            }
        }
}

// Transformation functions
template<typename T>
void GeneralizedPolynomialChaos<T>::chaosToRandom(const std::vector<T>& chaosCoefficients, std::vector<T>& randomVariables) {
    randomVariables.resize(totalNq, 0.0);

    for (size_t k = 0; k < totalNq; ++k) {
        auto startIt = phiRan.begin() + k * No;
        randomVariables[k] = std::inner_product(chaosCoefficients.begin(), chaosCoefficients.end(), startIt, 0.0);
    }
}

template<typename T>
void GeneralizedPolynomialChaos<T>::randomToChaos(const std::vector<T>& randomVariables, std::vector<T>& chaosCoefficients) {
    chaosCoefficients.resize(No, 0.0);
    std::vector<T> weightedRandomVariables(totalNq);

    // Compute weighted random variables
    for (size_t k = 0; k < totalNq; ++k) {
        weightedRandomVariables[k] = weightsMultiplied[k] * randomVariables[k];
    }

    // Compute chaos coefficients
    for (size_t i = 0; i < No; ++i) {
        auto startIt = phiRan_T.begin() + i * totalNq;
        chaosCoefficients[i] = std::inner_product(weightedRandomVariables.begin(), weightedRandomVariables.end(), startIt, 0.0);
        chaosCoefficients[i] *= t2Product_inv[i];
    }
}

// Chaos operations
template<typename T>
void GeneralizedPolynomialChaos<T>::chaosProduct(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& product) {
    product.resize(No, 0.0);
    std::vector<T> precomputedProductsFlat(No * No);

    for (size_t j = 0; j < No; ++j) {
        for (size_t k = 0; k < No; ++k) {
            precomputedProductsFlat[j * No + k] = chaos1[j] * chaos2[k];
        }
    }

    for (size_t i = 0; i < No; ++i) {
        T sum = 0.0;
        for (size_t j = 0; j < No; ++j) {
            for (size_t k = 0; k < No; ++k) {
                size_t flatIndex = i * No * No + j * No + k;
                sum += precomputedProductsFlat[j * No + k] * t3Product[flatIndex];
            }
        }
        product[i] = sum * t2Product_inv[i];
    }
}

template<typename T>
void GeneralizedPolynomialChaos<T>::chaosDivide(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& division) {
    std::vector<T> random1(totalNq, 0.0);
    std::vector<T> random2(totalNq, 0.0);
    chaosToRandom(chaos1, random1);
    chaosToRandom(chaos2, random2);
    randomToChaos(MatrixOperations::vectorDivision(random1, random2), division);
}

template<typename T>
void GeneralizedPolynomialChaos<T>::chaosSum(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& sum) {
    sum.resize(No);
    for (size_t i = 0; i < No; ++i) {
        sum[i] = chaos1[i] + chaos2[i];
    }
}

// Statistical moments
template<typename T>
T GeneralizedPolynomialChaos<T>::mean(const std::vector<T>& chaosCoefficients) {
    return chaosCoefficients[0];
}

template<typename T>
T GeneralizedPolynomialChaos<T>::std(const std::vector<T>& chaosCoefficients) {
    T variance = 0.0;
    for (size_t i = 1; i < No; ++i) {
        variance += t2Product[i] * chaosCoefficients[i] * chaosCoefficients[i];
    }
    return std::sqrt(variance);
}

template<typename T>
void GeneralizedPolynomialChaos<T>::convert2affinePCE(const Distribution<T>& distribution, std::vector<T>&chaos)
{
    switch (distribution.type) {
        case DistributionType::Uniform: {
            T a1 = 0.5 * (distribution.param1 + distribution.param2);
            T a2 = 0.5 * (distribution.param2 - distribution.param1);
            chaos[0] = a1;
            chaos[1] = a2;
            break;
        }
        case DistributionType::Normal: {
            chaos[0] = distribution.param1;
            chaos[1] = distribution.param2;
            break;
        }
        // Add cases for other distributions
        default:
            throw std::runtime_error("Unsupported distribution type for GPC.");
    }
}

// Getters
template<typename T>
size_t GeneralizedPolynomialChaos<T>::getPolynomialsOrder() const {
    return No;
}

template<typename T>
size_t GeneralizedPolynomialChaos<T>::getQuadraturePointsNumber() const {
    return totalNq;
}


template<typename T>
void GeneralizedPolynomialChaos<T>::getPointsAndWeights(std::vector<std::vector<T>>& points, std::vector<std::vector<T>>& weights) {
    points = this->points;
    weights = this->weights;
}

template<typename T>
std::vector<std::vector<T>> GeneralizedPolynomialChaos<T>::getStochasticCollocationSample() {
    std::vector<std::vector<T>> samples(totalNq, std::vector<T>(randomNumberDimension));

    for (size_t j = 0; j < randomNumberDimension; ++j) {
        for (size_t i = 0; i < totalNq; ++i) {
            samples[i][j] = affine(points[j][pointsWeightsIndexList[i][j]], distributions[j]);
        }
    }

    return samples;

}

template<typename T>
std::vector<T> GeneralizedPolynomialChaos<T>::getWeightsMultiplied() const {
    return weightsMultiplied;
}

template<typename T>
void GeneralizedPolynomialChaos<T>::getTensors(std::vector<T>& t2Product, std::vector<T>& t2Product_inv, std::vector<T>& t3Product) {
    t2Product = this->t2Product;
    t2Product_inv = this->t2Product_inv;
    t3Product = this->t3Product;
}

template<typename T>
std::shared_ptr<Polynomials::PolynomialBasis<T>> GeneralizedPolynomialChaos<T>::getPolynomialBasis(size_t dimension) const {
    if (dimension < 0 || dimension >= randomNumberDimension) {
        throw std::out_of_range("Dimension is out of bounds");
    }

    // Cast the void pointer back to the correct polynomial basis type
    return polynomialBases[dimension];
}

template<typename T>
std::vector<std::vector<size_t>> GeneralizedPolynomialChaos<T>::getMultiIndices() const {
    return inds;
}

template<typename T>
void GeneralizedPolynomialChaos<T>::getPhiRan(std::vector<T>& phiRan) {
    phiRan = this->phiRan;
}

template<typename T>
void GeneralizedPolynomialChaos<T>::getCoefficients(std::vector<std::vector<std::vector<T>>>& polynomialCoeffs) {
    polynomialCoeffs = this->coefficients;
}


// } // namespace uq

// } // namespace olb

#endif // GENERALIZED_POLYNOMIAL_CHAOS_HH