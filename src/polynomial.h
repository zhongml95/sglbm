// polynomial.h
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "util.h"
#include <gsl/gsl_integration.h>
#include "polynomial_coefficients.h"


class polynomial
{    

private:
    int nq;
    int order;
    int polynomialType;
    double parameter1;
    double parameter2;

    size_t points_weights_method;// 0: from householder, 1: from scipy

    std::vector<double> points;
    std::vector<double> weights;
    //std::vector<double> alpha;
    //std::vector<double> beta;
    std::vector<double> t2Product;
    std::vector<std::vector<double>> phi;
    std::vector<std::vector<double>> points_power;
    std::vector<std::vector<double>> phiRan;
    std::vector<std::vector<double>> polynomialCoeffs;
    std::vector<std::vector<std::vector<double>>> t3Product;

    PolynomialCoefficients polyCoeff;
    
    void initializeMatrices();
    void polynomial_chaos();
    void calc_points_power();
    void evaluatePhiRan();
    void tensor();
    
public:
    polynomial(int _nq, int _n, double _parameter1, double _parameter2, int _polynomialType, size_t _points_weights_method);

    int get_polynomial_order();

    int get_quadrature_points_number();

    void getPoints(std::vector<double>& _points);

    void getWeights(std::vector<double>& _weights);   

    void getPolynomialCoeffs(std::vector<std::vector<double>>& _polynomialCoeffs);

    double evaluate(int n, int k);

    //evaluate the polynomial at kth point 
    double evaluate_polynomial(int order_max, int k);

    void chaos2ran(const std::vector<double>& uChaos, std::vector<double>&uRan);

    void ran2chaos(const std::vector<double>& ran, std::vector<double>&chaos);

    void chaos_product(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>&production);

    void convert2affinePCE(double para1, double para2, std::vector<double>&domain);

    double mean(std::vector<double> chaos);

    double std(std::vector<double> chaos);

};


#endif