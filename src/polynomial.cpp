// legendre_basis.cpp
#include "polynomial.h"
#include "quadrature.h"


namespace Polynomials {

// Helper function to compute Legendre coefficients recursively
static std::vector<double> computeLegendreCoefficients(int n) {
    if (n == 0) {
        return {1.0};  // P_0(x) = 1
    } else if (n == 1) {
        return {0.0, 1.0};  // P_1(x) = x
    } else {
        // Recurrence relation: P_n(x) = ((2n - 1)/n) * x * P_{n-1}(x) - ((n - 1)/n) * P_{n-2}(x)
        std::vector<double> Pn_minus1 = computeLegendreCoefficients(n - 1);
        std::vector<double> Pn_minus2 = computeLegendreCoefficients(n - 2);
        std::vector<double> Pn(n + 1, 0.0);

        double a = (2.0 * n - 1.0) / n;
        double b = (n - 1.0) / n;

        // Multiply P_{n-1} by x (shift coefficients)
        std::vector<double> xPn_minus1(n + 1, 0.0);
        for (int i = 0; i < Pn_minus1.size(); ++i) {
            xPn_minus1[i + 1] = Pn_minus1[i];
        }

        // Compute Pn = a * x * P_{n-1} - b * P_{n-2}
        for (int i = 0; i <= n; ++i) {
            double term1 = a * xPn_minus1[i];
            double term2 = (i < Pn_minus2.size()) ? b * Pn_minus2[i] : 0.0;
            Pn[i] = term1 - term2;
        }

        return Pn;
    }
}

// Compute coefficients for the Legendre polynomial of order n
std::vector<double> LegendreBasis::computeCoefficients(int n) const {
    return computeLegendreCoefficients(n);
}

// Construct the Jacobi matrix of size n for Legendre polynomials
std::vector<std::vector<double>> LegendreBasis::constructJacobiMatrix(int n) const {
    if (n <= 1) {
        throw std::invalid_argument("The order (n) must be greater than 1.");
    }

    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    // Fill the Jacobi matrix according to the recurrence relation for Legendre polynomials
    for (int i = 1; i < n; ++i) {
        double a = i / std::sqrt(4.0 * i * i - 1.0);
        J[i][i - 1] = a;
        J[i - 1][i] = a;
    }

    return J;
}

// Evaluate the Legendre polynomial of order n at point x
double LegendreBasis::evaluatePolynomial(int n, double x) const {
    std::vector<double> coeffs = computeCoefficients(n);

    // Normalize the coefficients by the largest power term
    for (int i = 0; i <= n; ++i) {
        coeffs[i] /= coeffs[n];
    }

    double result = 0.0;
    for (int i = 0; i <= n; ++i) {
        result += coeffs[i] * std::pow(x, i);
    }
    return result;
}

// Compute the derivative of the Legendre polynomial of order n at point x
double LegendreBasis::derivativePolynomial(int n, double x) const {
    std::vector<double> coeffs = computeCoefficients(n);
    double result = 0.0;
    for (int i = 1; i <= n; ++i) {
        result += i * coeffs[i] * std::pow(x, i - 1);
    }
    return result;
}

double LegendreBasis::computeQuadratureWeight(int n, double x) const {
    double Pn_prime = derivativePolynomial(n - 1, x);
    return 2.0 / ((1.0 - x * x) * Pn_prime * Pn_prime);
}

std::shared_ptr<Quadrature::QuadratureBase> LegendreBasis::getQuadrature(int order, Quadrature::QuadratureMethod method) const {
    return std::make_shared<Quadrature::Quadrature<LegendreBasis>>(order, method);
}


// Helper function to compute Hermite coefficients recursively
static std::vector<double> computeHermiteCoefficients(int n) {
    if (n == 0) {
        return {1.0};  // H_0(x) = 1
    } else if (n == 1) {
        return {0.0, 2.0};  // H_1(x) = 2x
    } else {
        // Recurrence relation: H_n(x) = 2x * H_{n-1}(x) - 2(n - 1) * H_{n-2}(x)
        std::vector<double> Hn_minus1 = computeHermiteCoefficients(n - 1);
        std::vector<double> Hn_minus2 = computeHermiteCoefficients(n - 2);
        std::vector<double> Hn(n + 1, 0.0);

        // Compute 2x * H_{n-1}
        for (int i = 0; i < Hn_minus1.size(); ++i) {
            Hn[i + 1] += 2.0 * Hn_minus1[i];
        }

        // Compute -2(n - 1) * H_{n-2}
        double c = -2.0 * (n - 1);
        for (int i = 0; i < Hn_minus2.size(); ++i) {
            Hn[i] += c * Hn_minus2[i];
        }

        // Normalize the coefficients by the largest power term (optional)
        // for (int i = 0; i <= n; ++i) {
        //     Hn[i] /= Hn[n];
        // }

        return Hn;
    }
}

// Compute coefficients for the Hermite polynomial of order n
std::vector<double> HermiteBasis::computeCoefficients(int n) const {
    return computeHermiteCoefficients(n);
}

// Construct the Jacobi matrix of size n for Hermite polynomials
std::vector<std::vector<double>> HermiteBasis::constructJacobiMatrix(int n) const {
    if (n <= 1) {
        throw std::invalid_argument("The order (n) must be greater than 1.");
    }

    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    // Fill the Jacobi matrix according to the recurrence relation for Hermite polynomials
    for (int i = 0; i < n; ++i) {
        // The recurrence relation coefficients for Hermite polynomials
        if (i > 0) {
            J[i][i - 1] = std::sqrt(i / 2.0);
            J[i - 1][i] = J[i][i - 1];
        }
    }

    return J;
}

// Evaluate the Hermite polynomial of order n at point x
double HermiteBasis::evaluatePolynomial(int n, double x) const {
    // Use recursive evaluation or a direct formula for Hermite polynomials
    if (n == 0) {
        return 1.0;
    } else if (n == 1) {
        return 2.0 * x;
    } else {
        double Hn_minus2 = 1.0;
        double Hn_minus1 = 2.0 * x;
        double Hn = 0.0;
        for (int k = 2; k <= n; ++k) {
            Hn = 2.0 * x * Hn_minus1 - 2.0 * (k - 1) * Hn_minus2;
            Hn_minus2 = Hn_minus1;
            Hn_minus1 = Hn;
        }
        return Hn;
    }
}

// Compute the derivative of the Hermite polynomial of order n at point x
double HermiteBasis::derivativePolynomial(int n, double x) const {
    // Derivative relation: H_n'(x) = 2n * H_{n-1}(x)
    if (n == 0) {
        return 0.0;
    } else {
        return 2.0 * n * evaluatePolynomial(n - 1, x);
    }
}

double HermiteBasis::computeQuadratureWeight(int n, double x) const {
    double Hn_minus1 = evaluatePolynomial(n - 1, x);
    return std::exp(-x * x) / (Hn_minus1 * Hn_minus1);
}

std::shared_ptr<Quadrature::QuadratureBase> HermiteBasis::getQuadrature(int order, Quadrature::QuadratureMethod method) const {
    return std::make_shared<Quadrature::Quadrature<HermiteBasis>>(order, method);
}


} // namespace Polynomials
