// polynomials.h
#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <cmath>
#include <vector>
#include <iomanip>
#include "polynomial.h"

//#include "matrix.h"

class polynomials
{
private:
    size_t points_weights_method;
    int No;
    int nq;
    int total_nq;
    int order;
    int random_number_dimension;
    std::vector<int> polynomial_types;
    std::vector<double> parameter1;
    std::vector<double> parameter2;
    std::vector<std::vector<int>> inds;
    std::vector<std::vector<std::vector<double>>> polynomial_coefficients;
    std::vector<double> points;
    std::vector<double> weights;
    std::vector<double> weights_multiplied;
    std::vector<std::vector<int>> points_weights_index_list;
    
    std::vector<double> phiRan;
    std::vector<double> phiRan_T;
    std::vector<double> t2Product;
    std::vector<double> t2Product_inv;
    std::vector<double> t3Product;
    
    void initializeMatrices();

    void initializePolynomialCollection();

public:

    std::vector<polynomial> polynomial_collection;

    polynomials(int _nq, int _order, std::vector<double> _parameter1, std::vector<double> _parameter2, std::vector<int> _polynomial_types, size_t _points_weights_method);

    //n_order at point k
    double evaluate(int n_order, int k);
    
    //n_order at point k given polynomial basis 
    double evaluate(int n_order, int k, int phi_i);

    double evaluate(int n_order, std::vector<int> idx);

    double evaluate(int n_order, double x, int phi_i);

    //evaluate the polynomial at kth point 
    double evaluate_polynomial(int order_max, int k);

    void evaluatePhiRan();

    void tensor();

    void chaos2ran(const std::vector<double>& chaos, std::vector<double>&Ran);

    void ran2chaos(const std::vector<double>& ran, std::vector<double>&chaos);

    void chaos_product(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>&production);

    void chaos_sum(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>&sum);

    void convert2affinePCE(double para1, double para2, int polynomialType, std::vector<double>&domain);

    double mean(const std::vector<double>& chaos);

    double std(const std::vector<double>& chaos);

    int get_polynomials_order();
    int get_quadrature_points_number();
    double get_parameter1(int i);
    double get_parameter2(int i);
};


#endif // POLYNOMIALS_H