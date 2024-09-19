// generalized_polynomial_chaos.h
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_H
#define GENERALIZED_POLYNOMIAL_CHAOS_H

#include <vector>
#include <memory>
#include <cmath>
#include <string>

#include "utils.h"

// Include the polynomial basis and quadrature headers
#include "legendre_basis.h"
#include "hermite_basis.h"
#include "quadrature.h"

class GeneralizedPolynomialChaos {
public:
    // Constructor
    GeneralizedPolynomialChaos(int order,
                               int nq,
                               const std::vector<double>& parameter1,
                               const std::vector<double>& parameter2,
                               const std::vector<int>& polynomialTypes,
                               Quadrature::QuadratureMethod quadratureMethod);

    // Evaluation functions
    double evaluate(int n_order, int k);
    double evaluate(int n_order, int k, int phi_i);
    double evaluate(int n_order, const std::vector<int>& idx);
    double evaluate(int n_order, double x, int phi_i);
    double evaluate_polynomial(int order_max, int k);

    // Compute phiRan matrix
    void evaluatePhiRan();

    // Compute tensors
    void computeTensors();

    // Transformation functions
    void chaosToRandom(const std::vector<double>& chaosCoefficients, std::vector<double>& randomVariables);
    void randomToChaos(const std::vector<double>& randomVariables, std::vector<double>& chaosCoefficients);

    // Chaos operations
    void chaosProduct(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& product);
    void chaosSum(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& sum);

    // Statistical moments
    double mean(const std::vector<double>& chaosCoefficients);
    double std(const std::vector<double>& chaosCoefficients);

    void convert2affinePCE(double para1, double para2, int polynomialType, std::vector<double>&domain);

    // Getters
    int getPolynomialsOrder() const;
    int getQuadraturePointsNumber() const;
    double getParameter1(int i) const;
    double getParameter2(int i) const;
    void getPointsAndWeights(std::vector<std::vector<double>>& points, std::vector<std::vector<double>>& weights);
    void getTensors(std::vector<double>& t2Product, std::vector<double>& t2Product_inv, std::vector<double>& t3Product);
    std::vector<double> getWeightsMultiplied() const;
    // Template function to get the polynomial basis at a specific dimension (i)
    template <typename PolynomialBasis>
    std::shared_ptr<PolynomialBasis> getPolynomialBasis(int i) const;
    std::vector<std::vector<int>> getMultiIndices() const;

    void getPhiRan(std::vector<double>& phiRan);
    void getCoefficients(std::vector<std::vector<std::vector<double>>>& polynomialCoeffs);


    // void get_polynomial_coefficients(std::vector<std::vector<double>>& polynomialCoeffs) {
    //     polynomialCoeffs = this->polynomialCoeffs;
    // }


private:
    size_t pointsWeightsMethod;
    int No; // Number of polynomials
    int nq; // Number of quadrature points per dimension
    int totalNq; // Total number of quadrature points
    int order;
    int randomNumberDimension;
    std::vector<int> polynomialTypes;
    std::vector<double> parameter1;
    std::vector<double> parameter2;
    std::vector<std::vector<int>> inds; // Multi-indices
    std::vector<std::vector<double>> points;  // Points for each dimension
    std::vector<std::vector<double>> weights; // Weights for each dimension
    std::vector<double> weightsMultiplied;    // Combined weights
    std::vector<std::vector<int>> pointsWeightsIndexList;
    std::vector<std::vector<std::vector<double>>> coefficients; // Coefficients of polynomials

    Quadrature::QuadratureMethod quadratureMethod;

    std::vector<double> phiRan;   // Evaluated polynomials at quadrature points
    std::vector<double> phiRan_T; // Transpose of phiRan
    std::vector<double> t2Product;
    std::vector<double> t2Product_inv;
    std::vector<double> t3Product;

    // Polynomial bases and quadratures for each dimension
    std::vector<std::shared_ptr<void>> polynomialBases;
    std::vector<std::shared_ptr<void>> quadratures;

    // Initialization functions
    void initializePolynomialBases();
    void initializeQuadratures();
    void initializeMatrices();

    void initializePolynomialCoefficients();

    // Helper functions
    std::vector<int> findIndex(int idx, int dimension, int nq);
    void calculateMultiIndices(int d, int n, std::vector<std::vector<int>>& indices);

    // File I/O for tensor storage
    // bool fileExists(const std::string& filename);
    // void createDirectory(const std::string& directory);
    // void readVector1D(const std::string& filename, std::vector<double>& data);
    // void saveVector1D(const std::string& filename, const std::vector<double>& data);
};

#endif // GENERALIZED_POLYNOMIAL_CHAOS_H
