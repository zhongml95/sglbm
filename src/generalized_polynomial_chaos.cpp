// generalized_polynomial_chaos.cpp

#include "generalized_polynomial_chaos.h"
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

// Include the polynomial basis and quadrature headers
// #include "quadrature.h"

// Namespace aliases
using LegendreBasis = Polynomials::LegendreBasis;
using HermiteBasis = Polynomials::HermiteBasis;
// using Quadrature = Quadrature::Quadrature;

// Helper functions for file I/O and directory management
// namespace {
//     bool directoryExists(const std::string& directory) {
//         struct stat info;
//         if (stat(directory.c_str(), &info) != 0)
//             return false;
//         else if (info.st_mode & S_IFDIR)
//             return true;
//         else
//             return false;
//     }

//     void createDirectory(const std::string& directory) {
//         mkdir(directory.c_str(), 0777);
//     }

//     bool fileExists(const std::string& filename) {
//         std::ifstream infile(filename);
//         return infile.good();
//     }

//     void readVector1D(const std::string& filename, std::vector<double>& data) {
//         std::ifstream infile(filename, std::ios::binary);
//         if (!infile) {
//             throw std::runtime_error("Failed to open file for reading: " + filename);
//         }
//         size_t size = 0;
//         infile.read(reinterpret_cast<char*>(&size), sizeof(size));
//         data.resize(size);
//         infile.read(reinterpret_cast<char*>(data.data()), size * sizeof(double));
//         infile.close();
//     }

//     void saveVector1D(const std::string& filename, const std::vector<double>& data) {
//         std::ofstream outfile(filename, std::ios::binary);
//         if (!outfile) {
//             throw std::runtime_error("Failed to open file for writing: " + filename);
//         }
//         size_t size = data.size();
//         outfile.write(reinterpret_cast<const char*>(&size), sizeof(size));
//         outfile.write(reinterpret_cast<const char*>(data.data()), size * sizeof(double));
//         outfile.close();
//     }
// }

// Constructor
GeneralizedPolynomialChaos::GeneralizedPolynomialChaos(int _order,
                                                       int _nq,
                                                       const std::vector<std::shared_ptr<Polynomials::PolynomialBasis>>& _polynomialBases,
                                                       Quadrature::QuadratureMethod _quadratureMethod)
    : order(_order),
      nq(_nq),
      polynomialBases(_polynomialBases),
      quadratureMethod(_quadratureMethod) {

    randomNumberDimension = static_cast<int>(polynomialBases.size());

    // Calculate multi-indices
    calculateMultiIndices(randomNumberDimension, order, inds);
    No = static_cast<int>(inds.size());

    // Initialize polynomial coefficients, quadratures, and matrices
    initializePolynomialCoefficients();
    initializeQuadratures();
    initializeMatrices();

    // Evaluate polynomials at quadrature points
    evaluatePhiRan();

    // Compute tensors
    computeTensors();
}



// Initialize quadratures
void GeneralizedPolynomialChaos::initializeQuadratures() {
    points.resize(randomNumberDimension);
    weights.resize(randomNumberDimension);

    totalNq = static_cast<int>(std::pow(nq, randomNumberDimension));

    for (int i = 0; i < randomNumberDimension; ++i) {
        auto quadrature = polynomialBases[i]->getQuadrature(nq, quadratureMethod);
        points[i] = quadrature->getPoints();
        weights[i] = quadrature->getWeights();
    }
}



// Initialize matrices
void GeneralizedPolynomialChaos::initializeMatrices() {
    // std::cout << "Initializing matrices..." << std::endl;
    // Resize vectors
    phiRan.resize(totalNq * No, 0.0);
    phiRan_T.resize(totalNq * No, 0.0);
    t2Product.resize(No, 0.0);
    t2Product_inv.resize(No, 0.0);
    t3Product.resize(No * No * No, 0.0);


    // Generate pointsWeightsIndexList
    pointsWeightsIndexList.resize(totalNq, std::vector<int>(randomNumberDimension));
    for (int i = 0; i < totalNq; ++i) {
        pointsWeightsIndexList[i] = findIndex(i, randomNumberDimension, nq);
    }

    // Compute weightsMultiplied
    weightsMultiplied.resize(totalNq, 1.0);
    for (int k = 0; k < totalNq; ++k) {
        for (int dim = 0; dim < randomNumberDimension; ++dim) {
            int idx = pointsWeightsIndexList[k][dim];
            weightsMultiplied[k] *= weights[dim][idx];
        }
    }
}

// Initialize polynomial coefficients
void GeneralizedPolynomialChaos::initializePolynomialCoefficients() {
    coefficients.resize(randomNumberDimension);
    for (int phi_i = 0; phi_i < randomNumberDimension; ++phi_i) {
        auto basis = std::static_pointer_cast<LegendreBasis>(polynomialBases[phi_i]);

        coefficients[phi_i].resize(No);
        for (int i = 0; i < No; ++i) {
            coefficients[phi_i][i] = basis->computeCoefficients(i);
        }
    }

}

// Evaluate n_order polynomial at point k
double GeneralizedPolynomialChaos::evaluate(int n_order, int k) {
    double result = 1.0;
    for (int i = 0; i < randomNumberDimension; ++i) {
        result *= evaluate(inds[n_order][i], k, i);
    }
    return result;
}

// Evaluate n_order polynomial at point k and dimension phi_i
double GeneralizedPolynomialChaos::evaluate(int n_order, int k, int phi_i) {
    double x = points[phi_i][pointsWeightsIndexList[k][phi_i]];
    return evaluate(n_order, x, phi_i);
}

// Evaluate n_order polynomial at given multi-index
double GeneralizedPolynomialChaos::evaluate(int n_order, const std::vector<int>& idx) {
    double result = 1.0;
    for (int i = 0; i < randomNumberDimension; ++i) {
        result *= evaluate(inds[n_order][i], points[i][idx[i]], i);
    }
    return result;
}

// Evaluate polynomial basis at given order, point x, and dimension phi_i
double GeneralizedPolynomialChaos::evaluate(int n_order, double x, int phi_i) {
    if (phi_i < 0 || phi_i >= polynomialBases.size()) {
        throw std::out_of_range("Invalid dimension index phi_i.");
    }
    return polynomialBases[phi_i]->evaluatePolynomial(n_order, x);
}

// Evaluate the polynomial at kth point up to order_max
double GeneralizedPolynomialChaos::evaluate_polynomial(int order_max, int k) {
    double sum = 0.0;
    for (int i = 0; i <= order_max; ++i) {
        sum += evaluate(i, k);
    }
    return sum;
}

// Evaluate phiRan matrix
void GeneralizedPolynomialChaos::evaluatePhiRan() {
    std::cout << "Evaluating phiRan matrix..." << std::endl;
    for (int k = 0; k < totalNq; ++k) {
        for (int i = 0; i < No; ++i) {
            phiRan[k * No + i] = evaluate(i, pointsWeightsIndexList[k]);
            std::cout << phiRan[k * No + i] << " ";
            phiRan_T[i * totalNq + k] = phiRan[k * No + i];
        }
        std::cout << std::endl;
    }
}

// Helper functions
void GeneralizedPolynomialChaos::calculateMultiIndices(int d, int n, std::vector<std::vector<int>>& indices) {
    std::vector<int> index(d, 0);

    std::function<void(int, int, int)> recursiveFunction = [&](int pos, int sum, int maxOrder) {
        if (pos == d - 1) {
            index[pos] = maxOrder - sum;
            indices.push_back(index);
            return;
        }

        for (int i = 0; i <= maxOrder - sum; ++i) {
            index[pos] = i;
            recursiveFunction(pos + 1, sum + i, maxOrder);
        }
    };

    for (int order = 0; order <= n; ++order) {
        recursiveFunction(0, 0, order);
    }
}

std::vector<int> GeneralizedPolynomialChaos::findIndex(int idx, int dimension, int nq) {
    std::vector<int> index(dimension);
    for (int i = dimension - 1; i >= 0; --i) {
        index[i] = idx % nq;
        idx /= nq;
    }
    return index;
}

// Compute tensors (t2Product and t3Product)
void GeneralizedPolynomialChaos::computeTensors() {
    std::cout << "Computing tensors..." << std::endl;
    // File paths for saved matrices
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

    // Compute t2Product
    std::cout << "Computing t2Product..." << std::endl;
    if (fileExists(t2ProductFile)) {
        std::cout << "Loading t2Product from file." << std::endl;
        readVector1D(t2ProductFile, t2Product);
    } else {
        std::cout << "Calculating t2Product." << std::endl;
        // for (int i = 0; i < No; ++i) {
        //     double sum = 0.0;
        //     for (int k = 0; k < totalNq; ++k) {
        //         sum += phiRan[k * No + i] * phiRan[k * No + i] * weightsMultiplied[k];
        //     }
        //     t2Product[i] = sum;
        // }

        std::vector<std::vector<double>> tensor2d(No, std::vector<double>(No, 0.0));
        for (int i = 0; i < No; ++i) {
            for (int j = 0; j < No; ++j) {
                for (int m = 0; m < totalNq; ++m) {
                    // tensor2d[i][j] += evaluate(i, points_weights_index_list[m]) * evaluate(j, points_weights_index_list[m]) * weights_multiplied[m];
                    tensor2d[i][j] +=  phiRan[m * No + i] * phiRan[m * No + j] * weightsMultiplied[m];
                }
            }
        }
        for (int i = 0; i < No; ++i) {
            t2Product[i] = tensor2d[i][i];
        }

        saveVector1D(t2ProductFile, t2Product);
    }

    for (int i = 0; i < No; ++i) {
        t2Product_inv[i] = 1.0 / t2Product[i];
    }

    // Compute t3Product
    std::cout << "Computing t3Product..." << std::endl;
    if (fileExists(t3ProductFile)) {
        std::cout << "Loading t3Product from file." << std::endl;
        readVector1D(t3ProductFile, t3Product);
    } else {
        std::cout << "Calculating t3Product." << std::endl;
        for (int i = 0; i < No; ++i) {
            for (int j = 0; j < No; ++j) {
                for (int k = 0; k < No; ++k) {
                    double sum = 0.0;
                    for (int m = 0; m < totalNq; ++m) {
                        sum += phiRan[m * No + i] * phiRan[m * No + j] * phiRan[m * No + k] * weightsMultiplied[m];
                    }
                    t3Product[i * No * No + j * No + k] = sum;
                }
            }
        }
        saveVector1D(t3ProductFile, t3Product);
    }
}

// Transformation functions
void GeneralizedPolynomialChaos::chaosToRandom(const std::vector<double>& chaosCoefficients, std::vector<double>& randomVariables) {
    randomVariables.resize(totalNq, 0.0);

    for (int k = 0; k < totalNq; ++k) {
        auto startIt = phiRan.begin() + k * No;
        randomVariables[k] = std::inner_product(chaosCoefficients.begin(), chaosCoefficients.end(), startIt, 0.0);
    }
}

void GeneralizedPolynomialChaos::randomToChaos(const std::vector<double>& randomVariables, std::vector<double>& chaosCoefficients) {
    chaosCoefficients.resize(No, 0.0);
    std::vector<double> weightedRandomVariables(totalNq);

    // Compute weighted random variables
    for (int k = 0; k < totalNq; ++k) {
        weightedRandomVariables[k] = weightsMultiplied[k] * randomVariables[k];
    }

    // Compute chaos coefficients
    for (int i = 0; i < No; ++i) {
        auto startIt = phiRan_T.begin() + i * totalNq;
        chaosCoefficients[i] = std::inner_product(weightedRandomVariables.begin(), weightedRandomVariables.end(), startIt, 0.0);
        chaosCoefficients[i] *= t2Product_inv[i];
    }
}

// Chaos operations
void GeneralizedPolynomialChaos::chaosProduct(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& product) {
    product.resize(No, 0.0);
    std::vector<double> precomputedProductsFlat(No * No);

    for (int j = 0; j < No; ++j) {
        for (int k = 0; k < No; ++k) {
            precomputedProductsFlat[j * No + k] = chaos1[j] * chaos2[k];
        }
    }

    for (int i = 0; i < No; ++i) {
        double sum = 0.0;
        for (int j = 0; j < No; ++j) {
            for (int k = 0; k < No; ++k) {
                size_t flatIndex = i * No * No + j * No + k;
                sum += precomputedProductsFlat[j * No + k] * t3Product[flatIndex];
            }
        }
        product[i] = sum * t2Product_inv[i];
    }
}

void GeneralizedPolynomialChaos::chaosSum(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>& sum) {
    sum.resize(No);
    for (int i = 0; i < No; ++i) {
        sum[i] = chaos1[i] + chaos2[i];
    }
}

// Statistical moments
double GeneralizedPolynomialChaos::mean(const std::vector<double>& chaosCoefficients) {
    return chaosCoefficients[0];
}

double GeneralizedPolynomialChaos::std(const std::vector<double>& chaosCoefficients) {
    double variance = 0.0;
    for (int i = 1; i < No; ++i) {
        variance += t2Product[i] * chaosCoefficients[i] * chaosCoefficients[i];
    }
    return std::sqrt(variance);
}

void GeneralizedPolynomialChaos::convert2affinePCE(double para1, double para2, int polynomialType, std::vector<double>&domain)
{   
    if (polynomialType == 0) {
        double a1 = 0.5 * (para1 + para2);
        double a2 = 0.5 * (para2 - para1);
        domain[0] = a1;
        domain[1] = a2;
    }
    else if (polynomialType == 1) {
        double a1 = para1;
        double a2 = para2;
        domain[0] = a1;
        domain[1] = a2;
    }
}

// Getters
int GeneralizedPolynomialChaos::getPolynomialsOrder() const {
    return No;
}

int GeneralizedPolynomialChaos::getQuadraturePointsNumber() const {
    return totalNq;
}


void GeneralizedPolynomialChaos::getPointsAndWeights(std::vector<std::vector<double>>& points, std::vector<std::vector<double>>& weights) {
    points = this->points;
    weights = this->weights;
}

std::vector<double> GeneralizedPolynomialChaos::getWeightsMultiplied() const {
    return weightsMultiplied;
}

void GeneralizedPolynomialChaos::getTensors(std::vector<double>& t2Product, std::vector<double>& t2Product_inv, std::vector<double>& t3Product) {
    t2Product = this->t2Product;
    t2Product_inv = this->t2Product_inv;
    t3Product = this->t3Product;
}

template <typename PolynomialBasis>
std::shared_ptr<PolynomialBasis> GeneralizedPolynomialChaos::getPolynomialBasis(int dimension) const {
    if (dimension < 0 || dimension >= randomNumberDimension) {
        throw std::out_of_range("Dimension is out of bounds");
    }

    // Cast the void pointer back to the correct polynomial basis type
    return std::static_pointer_cast<PolynomialBasis>(polynomialBases[dimension]);
}

std::vector<std::vector<int>> GeneralizedPolynomialChaos::getMultiIndices() const {
    return inds;
}

void GeneralizedPolynomialChaos::getPhiRan(std::vector<double>& phiRan) {
    phiRan = this->phiRan;
}

void GeneralizedPolynomialChaos::getCoefficients(std::vector<std::vector<std::vector<double>>>& polynomialCoeffs) {
    polynomialCoeffs = this->coefficients;
}