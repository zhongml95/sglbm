#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <vector>
#include <memory>
#include <random>
#include <stdexcept>
#include <cmath>
#include <functional>
#include "distribution.h" // Include Distribution and DistributionType


// namespace olb {

// namespace uq {

template<typename T>
class MonteCarlo {
public:
    // Constructor for univariate distributions per dimension
    MonteCarlo(size_t numSamples, const std::vector<Distribution<T>>& distributions, unsigned int seed)
        : numSamples(numSamples), randomNumberDimension(distributions.size()), distributions(distributions), rng(seed) {
        if (distributions.size() != randomNumberDimension && distributions.size() != 1) {
            throw std::invalid_argument("Number of distributions must match the random number dimension or be 1 for multivariate distribution.");
        }
    }

    // Constructor for single distribution applied to all dimensions
    MonteCarlo(size_t numSamples, Distribution<T> distribution, unsigned int seed)
        : numSamples(numSamples), randomNumberDimension(1), rng(seed) {
        distributions = std::vector<Distribution<T>>(randomNumberDimension, distribution);
    }

    std::vector<std::vector<T>> generateSamples() {
        if (numSamples == 0 || randomNumberDimension == 0) {
            throw std::runtime_error("Number of samples and random number dimension must be greater than zero.");
        }

        if (distributions.empty()) {
            throw std::runtime_error("Distribution vector is empty.");
        }

        if (distributions.size() != randomNumberDimension) {
            throw std::runtime_error("Number of distributions must match the random number dimension.");
        }

        std::vector<std::vector<T>> samples(numSamples, std::vector<T>(randomNumberDimension));
        std::vector<std::function<T()>> randomGenerators(randomNumberDimension);

        for (size_t j = 0; j < randomNumberDimension; ++j) {
            const Distribution<T>& dist = distributions[j];

            switch (dist.type) {
                case DistributionType::Uniform: {
                    std::uniform_real_distribution<T> uniform_dist(dist.param1, dist.param2);
                    randomGenerators[j] = [this, uniform_dist = std::move(uniform_dist)]() mutable {
                        return uniform_dist(rng);
                    };
                    break;
                }
                case DistributionType::Normal: {
                    std::normal_distribution<T> normal_dist(dist.param1, dist.param2);
                    randomGenerators[j] = [this, normal_dist = std::move(normal_dist)]() mutable {
                        return normal_dist(rng);
                    };
                    break;
                }
                default:
                    throw std::runtime_error("Unsupported distribution type for univariate distribution.");
            }
        }

        // Generate samples
        for (size_t i = 0; i < numSamples; ++i) {
            for (size_t j = 0; j < randomNumberDimension; ++j) {
                samples[i][j] = randomGenerators[j]();
            }
        }

        return samples;
    }


    size_t getSamplesNumber() const {
        return numSamples;
    }

private:
    size_t numSamples;
    size_t randomNumberDimension;
    std::vector<Distribution<T>> distributions;
    std::mt19937 rng;  // Random number generator

};

// } // namespace uq

// } // namespace olb

#endif // MONTE_CARLO_H
