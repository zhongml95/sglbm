// polynomial.hh

#include "polynomial.h"
#include "polynomial_coefficients.h"
#include <algorithm>
#include <cmath>
#include <iostream>

// Implementation
polynomial::polynomial(int _nq, int _n, double _parameter1, double _parameter2, int _polynomialType, size_t _points_weights_method)
    : nq(_nq), order(_n), parameter1(_parameter1), parameter2(_parameter2), polynomialType(_polynomialType), points_weights_method(_points_weights_method) 
    {
    
        initializeMatrices();
        
        polynomial_chaos();

        calc_points_power();

        evaluatePhiRan();

        std::cout << "calc tensor" << std::endl;
        tensor();
    };

void polynomial::initializeMatrices() {
    points.resize(nq);
    weights.resize(nq);
    phi.resize(order + 1, std::vector<double>(nq, 0.0));
    polynomialCoeffs.resize(nq + 1, std::vector<double>(nq + 1, 0.0));
    phiRan.resize(nq, std::vector<double>(order + 1, 0.0));
    t2Product.resize(order + 1);
    t3Product.resize(order + 1, std::vector<std::vector<double>>(order + 1, std::vector<double>(order + 1, 0.0)));
}

void polynomial::calc_points_power() {
    // Precompute powers of points[k]
    points_power.resize(nq, std::vector<double>(order + 1, 0.0));
    for (int i = 0; i < nq; ++i) {
        for (int j = 0; j < order + 1; ++j) {
            points_power[i][j] = std::pow(points[i], j);
        }
    }
}

void polynomial::getPoints(std::vector<double>& _points) {
    _points = points;
}

void polynomial::getWeights(std::vector<double>& _weights) {
    _weights = weights;
}

int polynomial::get_polynomial_order() {
    return order;
}

int polynomial::get_quadrature_points_number() {
    return nq;
}

void polynomial::getPolynomialCoeffs(std::vector<std::vector<double>>& _polynomialCoeffs) {
    _polynomialCoeffs = polynomialCoeffs;
}

void polynomial::evaluatePhiRan() {
    
    std::cout << "calc PhiRan" << std::endl;
    for (int k = 0; k < nq; ++k) {
        // std::vector<double> pointPowers(order + 1, 1.0);

        for (int i = 0; i < order + 1; ++i) {
            for (int j = 0; j < order + 1; ++j) {
                phiRan[k][i] += polynomialCoeffs[i][j] * points_power[k][j];
            }
        }
    }
}

double polynomial::evaluate(int n, int k)
{       
    double sum = 0.0;
    for (int j = 0; j < n+1; ++j) {
        sum += polynomialCoeffs[n][j] * points_power[k][j];
    }
    return sum;
}

    //evaluate the polynomial at kth point 
double polynomial::evaluate_polynomial(int order_max, int k)
{
    double sum = 0.0;
    for (int i = 0; i < order_max + 1; ++i){
        sum += evaluate(i, k);
    }
    return sum;
}

void polynomial::tensor(){
    std::vector<std::vector<double>> tensor2d;
    tensor2d.resize(order + 1);
    for (int i = 0; i < order + 1; ++i) {
        tensor2d[i].resize(order + 1);
    }

    for (int i = 0; i < order + 1; ++i) {
        for (int j = 0; j < order + 1; ++j) {
            tensor2d[i][j] = 0.0;
            for (int m = 0; m < order + 1; m++) {
                t3Product[i][j][m] = 0.0;
            }
        }
    }

    for (int i = 0; i < order+1; ++i) {
        for (int j = 0; j < order+1; ++j) {
            for (int k = 0; k < nq; ++k) {
                tensor2d[i][j] += evaluate(i,k) * evaluate(j,k) * weights[k];
            }
        }
    }


    for (int i = 0; i < order + 1; ++i) {
        for (int j = 0; j < order + 1; ++j) {
            for (int m = 0; m < order + 1; m++) {
                for (int k = 0; k < nq; ++k) {
                    t3Product[i][j][m] += evaluate(i,k) * evaluate(j,k) * evaluate(m,k) * weights[k];
                }
            }
        }
    }
}

void polynomial::chaos2ran(const std::vector<double>& uChaos, std::vector<double>&uRan)
{        
    std::vector<double> _ran(nq, 0.0);

    for (int k = 0; k < nq; ++k){
        _ran[k] = std::inner_product(uChaos.begin(), uChaos.end(), phiRan[k].begin(), 0.0);
    }
    uRan.swap(_ran);
}

void polynomial::ran2chaos(const std::vector<double>& ran, std::vector<double>&chaos)
{    
    std::vector<double> _chaos(order+1, 0.0);

    for (int i = 0; i < order+1; ++i){
        for (int k = 0; k < nq; ++k){
            _chaos[i] += weights[k] * ran[k] * phiRan[i][k] / t2Product[i];
        }
    }
        
    chaos.swap(_chaos);
}

void polynomial::chaos_product(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>&production)
{
    for (int i = 0; i < order + 1; ++i) {
        double sum = 0.0;

        for (int j = 0; j < order + 1; ++j) {
            for (int k = 0; k < order + 1; ++k) {
            sum += chaos1[j] * chaos2[k] * t3Product[j][k][i];
            }
        }

        production[i] = sum / t2Product[i];
    }
}

void polynomial::polynomial_chaos()
{
    polyCoeff.coefficients(polynomialType, polynomialCoeffs);
    polyCoeff.get_points_weights(nq, order, polynomialType, polynomialCoeffs, points, weights);

    double weights_sum = 0.0;
    std::cout << "points" << std::endl;
    for (int i = 0; i < nq; i++) {
        std::cout << points[i] << "\t";
        weights_sum += weights[i];
    }
    std::cout << std::endl;
    std::cout << "weights" << std::endl;
    for (int i = 0; i < nq; i++) {
        std::cout << weights[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "sum of weights: " << weights_sum << std::endl;
}

    void polynomial::convert2affinePCE(double para1, double para2, std::vector<double>&domain)
    {   
        if (polynomialType == 0) {
            double a1 = 0.5 * (para1 + para2);
            double a2 = 0.5 * (para2 - para1);
            domain[0] = a1;// + alpha[0]  * a2;
            domain[1] = a2;
        }
        else if (polynomialType == 1) {
            double a1 = para1;
            double a2 = para2;
            domain[0] = a1;
            domain[1] = a2;
        }
    }

    double polynomial::mean(std::vector<double> chaos)
    {
        return chaos[0];
    }

    double polynomial::std(std::vector<double> chaos)
    {
        double sum = 0.0;
        for (int i = 1; i < order+1; ++i){
            sum += t2Product[i] * chaos[i] * chaos[i];
        }
        return std::sqrt(sum);
    }
