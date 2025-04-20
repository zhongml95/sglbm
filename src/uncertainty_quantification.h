#ifndef UNCERTAINTY_QUANTIFICATION_H
#define UNCERTAINTY_QUANTIFICATION_H

#include <vector>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>

#include "distribution.h" // Include the Distribution definitions

// Include the generalized polynomial chaos header
#include "generalized_polynomial_chaos.h"

// Include the sampling methods
#include "monte_carlo.h"
#include "quasi_monte_carlo.h"
#include "latin_hypercube_sampling.h"


// namespace olb {

// namespace uq {

// Enumeration to specify the uncertainty quantification method
enum class UQMethod {
    GPC,
    MonteCarlo,
    QuasiMonteCarlo,
    LatinHypercubeSampling
};

template<typename T>
class UncertaintyQuantification {
public:
    // Constructor
    UncertaintyQuantification(UQMethod uqMethod);

    // Delete copy constructor and copy assignment operator to make the class non-copyable
    UncertaintyQuantification(const UncertaintyQuantification&) = delete;
    UncertaintyQuantification& operator=(const UncertaintyQuantification&) = delete;

    // Default move constructor and move assignment operator
    UncertaintyQuantification(UncertaintyQuantification&&) = default;
    UncertaintyQuantification& operator=(UncertaintyQuantification&&) = default;

    // Destructor
    ~UncertaintyQuantification() = default;

    // Initialization functions
    void initializeGPC(size_t order, size_t nq,
                        Distribution<T> distribution,
                        Quadrature::QuadratureMethod quadratureMethod = Quadrature::QuadratureMethod::WilkinsonShiftQR);
    void initializeGPC(size_t order, size_t nq,
                       const std::vector<Distribution<T>>& distributions,
                       Quadrature::QuadratureMethod quadratureMethod = Quadrature::QuadratureMethod::WilkinsonShiftQR);

    void initializeMonteCarlo(size_t numSamples, const std::vector<Distribution<T>>& distributions, unsigned int seed);
    void initializeMonteCarlo(size_t numSamples, Distribution<T> distribution, unsigned int seed );

    void initializeQuasiMonteCarlo(size_t numSamples,
                                    Distribution<T> distribution,
                                    const std::string& dir_file = "new-joe-kuo-6.21201",
                                    GeneratorType generator = GeneratorType::Sobol );
    void initializeQuasiMonteCarlo(size_t numSamples,
                                   const std::vector<Distribution<T>>& distributions,
                                   const std::string& dir_file = "new-joe-kuo-6.21201",
                                   GeneratorType generator = GeneratorType::Sobol );

    void initializeLatinHypercubeSampling(size_t numSamples, size_t randomNumberDimension);

    // Function to get sampling points
    std::vector<std::vector<T>> getSamplingPoints();
    size_t getSamplesNumber();

    // Statistical moments
    T mean(const std::vector<T>& samples);
    T std(const std::vector<T>& samples);

    // Other methods common to all UQ methods
    // ...

    // getters
    std::unique_ptr<GeneralizedPolynomialChaos<T>> getOps() {
        return std::move(ops);
    }

private:
    UQMethod uqMethod;
    size_t numSamples;              // Number of samples (for MC, QMC, LHS)
    size_t randomNumberDimension;   // Dimensionality of the random input
    std::vector<std::vector<T>> points;  // Sampling points
    std::vector<Distribution<T>> distributions;

    // GPC-specific members
    size_t order;                   // Order of polynomials (for GPC)
    size_t No;                      // total order of polynomials system (for GPC)
    size_t nq;                      // Number of quadrature points per dimension (for GPC)
    std::unique_ptr<GeneralizedPolynomialChaos<T>> ops;
    std::vector<std::vector<T>> weights; // Weights for quadrature (GPC)
    std::vector<T> weightsMultiplied;    // Combined weights (GPC)
    std::vector<std::vector<size_t>> multiIndices; // Multi-indices (GPC)
    std::vector<std::shared_ptr<Polynomials::PolynomialBasis<T>>> polynomialBases;
    Quadrature::QuadratureMethod quadratureMethod; // Quadrature method (GPC)

    // Sampling method instances
    std::unique_ptr<MonteCarlo<T>> monteCarlo;
    std::unique_ptr<QuasiMonteCarlo<T>> quasiMonteCarlo;
    std::unique_ptr<LatinHypercubeSampling<T>> lhs;

    // Random number generator
    std::mt19937 rng; // Random number generator

    // Sample generation functions are now in respective classes

    // Example evaluation function can be defined here or elsewhere
};

// } // namespace uq

// } // namespace olb

#endif // UNCERTAINTY_QUANTIFICATION_H
