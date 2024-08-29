// polynomial_coefficients.h
#ifndef POLYNOMIAL_COEFFICIENTS_H
#define POLYNOMIAL_COEFFICIENTS_H

#include <vector>

class PolynomialCoefficients {

public:
    double derivative(int n, double point, int order, std::vector<double> polynomialCoeffs);
    void coefficients(int polynomialType, std::vector<std::vector<double>>& polynomialCoeffs);
    std::vector<double> Legendre_coefficients_order_m(int m, std::vector<double> v, std::vector<double> u);    
    std::vector<double> Hermite_coefficients_order_m(int m);
    void Legendre_coefficients(std::vector<std::vector<double>>& polynomialCoeffs);
    void Hermite_coefficients(std::vector<std::vector<double>>& polynomialCoeffs);

    void get_points_weights(int n, int order, int polynomialType, std::vector<std::vector<double>> polynomialCoeffs, std::vector<double>& points, std::vector<double>& weights);
    void get_gauss_legendre_points_weights_GSL(int n, std::vector<double>& points, std::vector<double>& weights);
    void get_gauss_legendre_points_weights_Householder(int n, int order, std::vector<std::vector<double>> polynomialCoeffs, std::vector<double>& points, std::vector<double>& weights);

};

#endif // POLYNOMIAL_COEFFICIENTS_H
