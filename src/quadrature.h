// quadrature.h
#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>
#include <stdexcept>

// Uncomment the following line if you have GSL installed and want to use it
// #define USE_GSL

#ifdef USE_GSL
#include <gsl/gsl_integration.h>
#endif

namespace Quadrature {

enum class QuadratureMethod {
    HouseholderQR,
    GSL
};

template <typename PolynomialBasis>
class Quadrature {
public:
    Quadrature(int order, QuadratureMethod method = QuadratureMethod::HouseholderQR)
        : order(order), basis(std::make_shared<PolynomialBasis>()), method(method) {
        computeQuadrature();
    }

    const std::vector<double>& getPoints() const {
        return points;
    }

    const std::vector<double>& getWeights() const {
        return weights;
    }

private:
    int order;
    std::shared_ptr<PolynomialBasis> basis;
    QuadratureMethod method;
    std::vector<double> points;
    std::vector<double> weights;

void computeQuadrature() {
#ifdef USE_GSL
    if (method == QuadratureMethod::GSL) {
        std::cout << "Using GSL quadrature method." << std::endl;
        computeQuadratureGSL();
    } else if (method == QuadratureMethod::HouseholderQR) {
        std::cout << "Using HouseholderQR quadrature method." << std::endl;
        computeQuadratureHouseholderQR();
    } else {
        std::cerr << "Warning: Unsupported quadrature method. Defaulting to HouseholderQR." << std::endl;
        computeQuadratureHouseholderQR();
    }
#else
    if (method != QuadratureMethod::HouseholderQR) {
        std::cerr << "Warning: GSL is not enabled. Falling back to HouseholderQR quadrature method." << std::endl;
    }
    std::cout << "Using HouseholderQR quadrature method." << std::endl;
    computeQuadratureHouseholderQR();
#endif
}


#ifdef USE_GSL
    void computeQuadratureGSL() {
        // Check if PolynomialBasis supports GSL
        if constexpr (std::is_same_v<PolynomialBasis, Polynomials::Legendre::LegendreBasis>) {
            computeQuadraturePointsWeightsGSL(order, points, weights);

            // Normalize weights if necessary
            double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
            for (auto& w : weights) {
                w /= sum;
            }
        } else {
            throw std::runtime_error("GSL quadrature method is only supported for Legendre polynomials.");
        }
    }

    void computeQuadraturePointsWeightsGSL(int n, std::vector<double>& points, std::vector<double>& weights) {
        if (n <= 0) {
            throw std::invalid_argument("Number of quadrature points must be positive.");
        }

        points.resize(n);
        weights.resize(n);

        gsl_integration_glfixed_table* table = gsl_integration_glfixed_table_alloc(n);

        double xi, wi;
        for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
            gsl_integration_glfixed_point(-1.0, 1.0, i, &xi, &wi, table);
            points[i] = xi;
            weights[i] = wi * 0.5; // Adjust weights if necessary
        }

        gsl_integration_glfixed_table_free(table);
    }
#endif

    void computeQuadratureHouseholderQR() {
        // Step 1: Construct the Jacobi matrix using the basis
        auto J = basis->constructJacobiMatrix(order);

        // Step 2: Compute eigenvalues (quadrature points) using Householder method
        std::vector<std::vector<double>> R;
        computeEigenvalues(J, R);

        // Extract eigenvalues (diagonal elements of R)
        points.resize(order);
        for (int i = 0; i < order; ++i) {
            points[i] = R[i][i];
        }

        // Sort the points
        std::sort(points.begin(), points.end());

        // Step 3: Compute quadrature weights
        computeWeights();
    }

    void computeWeights() {
        weights.resize(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
            double x = points[i];
            weights[i] = computeQuadratureWeight(*basis, order, x);
        }

        // Normalize weights if necessary
        double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        for (double& w : weights) {
            w /= sum;
        }
    }

    // Templated computeQuadratureWeight function
    template <typename Basis>
    double computeQuadratureWeight(const Basis& basis, int n, double x) {
        static_assert(!std::is_same_v<Basis, Basis>, "Unsupported polynomial basis for quadrature weights.");
        return 0.0;
    }

    // Specialization for LegendreBasis
    double computeQuadratureWeight(const Polynomials::Legendre::LegendreBasis& basis, int n, double x) {
        double Pn_prime = basis.derivativePolynomial(n - 1, x);
        return 2.0 / ((1.0 - x * x) * Pn_prime * Pn_prime);
    }

    // Specialization for HermiteBasis
    double computeQuadratureWeight(const Polynomials::Hermite::HermiteBasis& basis, int n, double x) {
        double Hn_minus1 = basis.evaluatePolynomial(n - 1, x);
        return std::exp(-x * x) / (Hn_minus1 * Hn_minus1);
    }

    // Householder QR decomposition and helper functions
    void computeEigenvalues(std::vector<std::vector<double>>& J, std::vector<std::vector<double>>& R) {
        int n = J.size();
        std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> J_new = J;
        int iter = 0;
        int max_iter = 1000;
        double tol = 1e-12;

        while (iter < max_iter) {
            householderQR(J_new, Q, R);
            J_new = matrixMultiplication(R, Q);

            if (iter % 10 == 0 && isConverged(J, J_new, tol)) {
                break;
            }

            J = J_new;
            iter++;
        }
    }

    void householderQR(const std::vector<std::vector<double>>& A,
                       std::vector<std::vector<double>>& Q,
                       std::vector<std::vector<double>>& R) {
        int m = A.size();
        int n = A[0].size();
        R = A;
        Q = generateIdentityMatrix(m);

        for (int k = 0; k < n; ++k) {
            std::vector<double> x(m - k);
            for (int i = k; i < m; ++i) {
                x[i - k] = R[i][k];
            }

            double norm_x = vectorNorm(x);
            double sign = (x[0] >= 0) ? 1.0 : -1.0;
            double u1 = x[0] + sign * norm_x;
            std::vector<double> v = x;
            v[0] = u1;

            double s = vectorNorm(v);
            if (s != 0.0) {
                for (auto& vi : v) {
                    vi /= s;
                }
            }

            std::vector<std::vector<double>> H = generateIdentityMatrix(m);
            for (int i = k; i < m; ++i) {
                for (int j = k; j < m; ++j) {
                    H[i][j] -= 2.0 * v[i - k] * v[j - k];
                }
            }

            R = matrixMultiplication(H, R);
            Q = matrixMultiplication(Q, H);
        }
    }

    double vectorNorm(const std::vector<double>& v) {
        double sum = 0.0;
        for (double vi : v) {
            sum += vi * vi;
        }
        return std::sqrt(sum);
    }

    std::vector<std::vector<double>> generateIdentityMatrix(int n) {
        std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            I[i][i] = 1.0;
        }
        return I;
    }

    std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>>& A,
                                                          const std::vector<std::vector<double>>& B) {
        int m = A.size();
        int n = B[0].size();
        int p = A[0].size();
        std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                double sum = 0.0;
                for (int k = 0; k < p; ++k) {
                    sum += A[i][k] * B[k][j];
                }
                C[i][j] = sum;
            }
        }
        return C;
    }

    bool isConverged(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B, double tol) {
        int n = A.size();
        double norm_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            norm_diff += std::abs(A[i][i] - B[i][i]);
        }
        return norm_diff < tol;
    }
};

} // namespace Quadrature

#endif // QUADRATURE_H
