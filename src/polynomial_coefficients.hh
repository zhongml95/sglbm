// polynomial_coefficients.hh

#include <vector>
#include "polynomial_coefficients.h"
#include <gsl/gsl_integration.h>
#include "matrix.h"

    double PolynomialCoefficients::derivative(int n, double point, int order, std::vector<double> polynomialCoeffs) {
        if (n == 0) {
            return 0.0;
        } 
        else {
            double sum = 0.0;
            for (int i = 1; i < order+1; ++i) {
                sum += polynomialCoeffs[i] * i * pow(point, i-1);
            }
            return sum;
        }
    }

std::vector<double> PolynomialCoefficients::Legendre_coefficients_order_m(int m, std::vector<double> v, std::vector<double> u)
    {
        if (m == 0) 
        {               
            return {1};
        }
        if (m == 1)
        {
            return {0, 1};
        }

        // Initialize with zero, only at most (half + 1) of the terms will be changed later
        std::vector<double> coeffs(m + 1, 0.0);

        // using literals of floating point type, 'm' is promoted by the compiler
        double a = (2.0 * m - 1.0) / m;
        double b = (m - 1.0) / m;

        int first = 1;
        // If 'm' is even, so is (m - 2) and its first element is zero. It can be skipped.
        // It also avoids annoying '-0' in the output 
        if ( m % 2 == 0 )
        {
            coeffs[0] = -b * u[0];
            first = 2;
        }
        for (int i = first; i < m - 1; i += 2)
        {
            coeffs[i] = (a * v[i - 1] - b * u[i]);
        }
        coeffs[m] = a * v[m - 1];
        return coeffs;
    }

    // Function to generate Hermite polynomials up to a given order
    std::vector<double>  PolynomialCoefficients::Hermite_coefficients_order_m(int m) {

        std::vector<double> coeffs(m + 1, 0.0);
        if (m == 0) {               
            return {1};
        }
        else if (m == 1) {
            return {0, 1};
        }
        else {                
            std::vector<double> m_1 = Hermite_coefficients_order_m(m - 1);
            //std::vector<double> m_2 = Hermite_coefficients(m - 2);
            for(int i =0; i < m; ++i) {
                coeffs[i+1] = m_1[i];
            }
            for(int i = 0; i < m; ++i) {
                coeffs[i] = coeffs[i] - (i+1) * m_1[i+1];
            }
            return coeffs;
        }
    }

    void PolynomialCoefficients::Legendre_coefficients(std::vector<std::vector<double>>& polynomialCoeffs) {
        int nq = polynomialCoeffs.size() - 1;
        for (int i = 0; i < nq+1; ++i) {
            std::vector<double> coeffs_i(nq+1, 0.0);
            if (i == 0 || i == 1) {
                coeffs_i = Legendre_coefficients_order_m(i, polynomialCoeffs[i], polynomialCoeffs[i]);
            } else {
                coeffs_i = Legendre_coefficients_order_m(i, polynomialCoeffs[i-1], polynomialCoeffs[i-2]);
            }
            for (int j = 0; j < nq+1; ++j) {
                if (j < i+1) {
                    polynomialCoeffs[i][j] = coeffs_i[j];
                } else {
                    polynomialCoeffs[i][j] = 0.;
                }
            }
        }

        // Normalize the coefficients
        for (int i = 0; i < nq+1; ++i) {
            for (int j = 0; j < i+1; ++j) {
                polynomialCoeffs[i][j] = polynomialCoeffs[i][j] / polynomialCoeffs[i][i];
            }
        }
    }

    void PolynomialCoefficients::Hermite_coefficients(std::vector<std::vector<double>>& polynomialCoeffs) {
        int nq = polynomialCoeffs.size();
        for (int i = 0; i < nq; ++i) {
            std::vector<double> coeffs_i(nq, 0.0);
            coeffs_i = Hermite_coefficients_order_m(i);

            for (int j = 0; j < i+1; ++j) {
                polynomialCoeffs[i][j] = coeffs_i[j] / coeffs_i[i];
            }
        }

        for (int i = 0; i < nq; ++i) {
            for (int j = 0; j < i+1; ++j) {
                polynomialCoeffs[i][j] = polynomialCoeffs[i][j] / polynomialCoeffs[i][i];
            }
        }
    }

    // Function to return Gauss-Legendre quadrature points using GSL
    void PolynomialCoefficients::get_gauss_legendre_points_weights_GSL(int n, std::vector<double>& points, std::vector<double>& weights) {
        gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(n);

        double xi, wi;
        for (size_t i = 0; i < n; ++i) {
            gsl_integration_glfixed_point(-1.0, 1.0, i, &xi, &wi, table);
            points[i] = xi;
            weights[i] = wi * 0.5;
        }

        gsl_integration_glfixed_table_free(table);
    }

    // Function to return Gauss-Legendre quadrature points using Householder method
    void PolynomialCoefficients::get_gauss_legendre_points_weights_Householder(int nq, int order, std::vector<std::vector<double>> polynomialCoeffs, std::vector<double>& points, std::vector<double>& weights) {
        std::vector<std::vector<double>> J = constructJacobiMatrix_Legendre(nq);
        std::vector<std::vector<double>> R = J;

        Householder(J, R);

        for (int i = 0; i < nq; i++) {
            points[i] = R[i][i];
        }
        std::sort(points.begin(), points.end());

        for (int i = 0; i < nq; ++i) {
            std::cout << "derivative: " << derivative(nq, points[i], order, polynomialCoeffs[nq]) << std::endl;
            weights[i] = 1.0 / ((1.0 - points[i] * points[i]) * std::pow(derivative(nq, points[i], order, polynomialCoeffs[nq]), 2.0));
        }

        // Normalize the weights
        double sum = 0.0;
        for (int i = 0; i < nq; ++i) {
            sum += weights[i];
        }
        for (int i = 0; i < nq; ++i) {
            weights[i] = weights[i] / sum;
        }

    }

    void PolynomialCoefficients::coefficients(int polynomialType, std::vector<std::vector<double>>& polynomialCoeffs) {
        if (polynomialType == 0) {
            Legendre_coefficients(polynomialCoeffs);
        } else if (polynomialType == 1) {
            Hermite_coefficients(polynomialCoeffs);
        }
    }

    void PolynomialCoefficients::get_points_weights(int n, int order, int polynomialType, std::vector<std::vector<double>> polynomialCoeffs, std::vector<double>& points, std::vector<double>& weights) {
        if (polynomialType == 0) {
            get_gauss_legendre_points_weights_Householder(n, order, polynomialCoeffs, points, weights);
        } else if (polynomialType == 1) {
            get_gauss_legendre_points_weights_GSL(n, points, weights);
        }
    }

