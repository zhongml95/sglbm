// generalized_polynomial_chaos.h
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_H
#define GENERALIZED_POLYNOMIAL_CHAOS_H

#include <vector>
#include <memory>
#include <cmath>
#include <string>

#include "utils.h"

#include "distribution.h"

// Include the polynomial basis and quadrature headers
#include "polynomial.h"

// #include "quadrature.h"

// namespace olb {

// namespace uq {
    
template<typename T>
class GeneralizedPolynomialChaos {
public:
    // Constructor
    GeneralizedPolynomialChaos(size_t order,
                               size_t nq,
                               const std::vector<Distribution<T>>& distributions,
                               Quadrature::QuadratureMethod quadratureMethod);

    // Evaluation functions
    T evaluate(size_t n_order, size_t k);
    T evaluate(size_t n_order, size_t k, size_t phi_i);
    T evaluate(size_t n_order, const std::vector<size_t>& idx);
    T evaluate(size_t n_order, T x, size_t phi_i);
    T evaluate_polynomial(size_t order_max, size_t k);

    // Compute phiRan matrix
    void evaluatePhiRan();

    // Compute tensors
    void computeTensors();

    // Transformation functions
    void chaosToRandom(const std::vector<T>& chaosCoefficients, std::vector<T>& randomVariables);
    void randomToChaos(const std::vector<T>& randomVariables, std::vector<T>& chaosCoefficients);

    // Chaos operations
    void chaosProduct(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& product);
    void chaosDivide(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& division);
    void chaosSum(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& sum);

    // Statistical moments
    T mean(const std::vector<T>& chaosCoefficients);
    T std(const std::vector<T>& chaosCoefficients);

    void convert2affinePCE(const Distribution<T>& distribution, std::vector<T>&chaos);

    // Getters
    size_t getPolynomialsOrder() const;
    size_t getQuadraturePointsNumber() const;
    void getPointsAndWeights(std::vector<std::vector<T>>& points, std::vector<std::vector<T>>& weights);
    std::vector<std::vector<T>> getStochasticCollocationSample();
    void getTensors(std::vector<T>& t2Product, std::vector<T>& t2Product_inv, std::vector<T>& t3Product);
    std::vector<T> getWeightsMultiplied() const;
    // Template function to get the polynomial basis at a specific dimension (i)
    std::shared_ptr<Polynomials::PolynomialBasis<T>> getPolynomialBasis(size_t i) const;
    std::vector<std::vector<size_t>> getMultiIndices() const;
    const std::vector<Distribution<T>>& getDistributions() const {
        return distributions;
    }

    void getPhiRan(std::vector<T>& phiRan);
    void getCoefficients(std::vector<std::vector<std::vector<T>>>& polynomialCoeffs);


    // void get_polynomial_coefficients(std::vector<std::vector<T>>& polynomialCoeffs) {
    //     polynomialCoeffs = this->polynomialCoeffs;
    // }


private:
    size_t pointsWeightsMethod;
    size_t No; // Number of polynomials
    size_t nq; // Number of quadrature points per dimension
    size_t totalNq; // Total number of quadrature points
    size_t order;
    size_t randomNumberDimension;
    std::vector<std::vector<size_t>> inds; // Multi-indices
    std::vector<std::vector<T>> points;  // Points for each dimension
    std::vector<std::vector<T>> weights; // Weights for each dimension
    std::vector<std::vector<T>> pointsTensor; // Tensor product of points
    std::vector<T> weightsMultiplied;    // Combined weights
    std::vector<std::vector<size_t>> pointsWeightsIndexList;
    std::vector<std::vector<std::vector<T>>> coefficients; // Coefficients of polynomials

    Quadrature::QuadratureMethod quadratureMethod;

    std::vector<T> phiRan;   // Evaluated polynomials at quadrature points
    std::vector<T> phiRan_T; // Transpose of phiRan
    std::vector<T> t2Product;
    std::vector<T> t2Product_inv;
    std::vector<T> t3Product;

    // Distributions for each dimension
    std::vector<Distribution<T>> distributions;

    // Polynomial bases for each dimension
    std::vector<std::shared_ptr<Polynomials::PolynomialBasis<T>>> polynomialBases;

    // Initialization functions
    void initializeQuadratures();
    void initializeMatrices();

    void initializePolynomialCoefficients();

    // Helper functions
    std::vector<size_t> findIndex(size_t idx, size_t dimension, size_t nq);
    void calculateMultiIndices(size_t d, size_t n, std::vector<std::vector<size_t>>& indices);

    std::shared_ptr<Polynomials::PolynomialBasis<T>> createPolynomialBasis(const Distribution<T>& dist) {
        switch (dist.type) {
            case DistributionType::Uniform:
                return std::make_shared<Polynomials::LegendreBasis<T>>();
            case DistributionType::Normal:
                return std::make_shared<Polynomials::HermiteBasis<T>>();
            // Add cases for other distributions
            default:
                throw std::runtime_error("Unsupported distribution type for GPC.");
        }
    }

};

// } // namespace uq

// } // namespace olb

#endif // GENERALIZED_POLYNOMIAL_CHAOS_H
