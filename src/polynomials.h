// polynomials.h
#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <cmath>
#include <vector>
#include <iomanip>
#include "polynomial.h"

//#include "matrix.h"


int numberPolynomials(int d, int n) {
    return int(std::tgamma(d+n+1) + 0.5) / int(std::tgamma(d+1) * std::tgamma(n+1)  + 0.5 );
}

std::vector<std::vector<int>> calculateMultiIndices(int d, int n) {
    int No = numberPolynomials(d, n);
    std::cout << "No = " << No << std::endl;
    std::vector<std::vector<int>> inds(No, std::vector<int>(d, 0));
    std::vector<std::vector<int>> pi(No, std::vector<int>(d, 1));
    // zeroth order polynomial
    if (n == 0) {
        return inds;
    }
    else if (n == 1) {
        // first order polynomial
        for (int j = 0; j < d; ++j) {
            inds[j+1][j] = 1;
        }

        // print index library
        // for (int i = 0; i < No; ++i) {
        //     for (int j = 0; j < d; ++j) {
        //         std::cout << inds[i][j] << "\t";
        //     }
        //     std::cout << std::endl;
        // }

        return inds;
    }
    else {
        // first order polynomial
        for (int j = 0; j < d; ++j) {
            inds[j+1][j] = 1;
        }

        for (int k = 1; k < No; ++k) {
            int g = 0;
            for (int l = 0; l < d; ++l) {
                int sum = 0;
                for (int m = 0; m < d; ++m) {
                    sum += pi[k-1][m];
                }
                pi[k][l] = sum - g;
                g = g + pi[k-1][l];
            }
        }

        int P = d+1;
        for (int k  = 1; k < n; ++k) {
            int L = P;
            for (int j = 0; j < d; ++j) {
                for (int m = L - pi[k][j]; m < L; ++m) {
                    P += 1;

                    for(int i = 0; i < d; ++i) {
                        inds[P-1][i] = inds[m][i];
                    }
                    inds[P-1][j] = inds[P-1][j] + 1;
                }
            }
        }

        // print index library
        // std::cout << "inds" << std::endl;
        // for(int i = 0; i < No; ++i) {
        //     for (int j = 0; j < d; ++j) {
        //         std::cout << inds[i][j] << "\t";
        //     }
        //     std::cout << std::endl;
        // }
        return inds;
    }
}


std::vector<int> find_index(int idx, int dimension, int nq) {
    std::vector<int> index;
    for (int i = 0; i < dimension; ++i) {
        int divisor = static_cast<int>(std::pow(nq, dimension - i - 1));
        index.push_back(idx / divisor);
        idx = idx % divisor;
    }
    return index;
}


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

    void evaluatePhiRan() {        
        for (int k = 0; k < total_nq; ++k) {
            for (int i = 0; i < No; ++i) {
                // std::vector<int> idx = find_index(k, random_number_dimension, nq);
                phiRan[k * No + i] = evaluate(i, points_weights_index_list[k]);
                phiRan_T[i * total_nq + k] = phiRan[k * No + i];
            }
        }
    }


    void tensor(){

        // File paths for saved matrices
        const std::string directoryT2Product = "../../src/t2Product/";
        if (!directoryExists(directoryT2Product)) {
            createDirectory(directoryT2Product);
        }
        const std::string directoryT3Product = "../../src/t3Product/";
        if (!directoryExists(directoryT3Product)) {
            createDirectory(directoryT3Product);
        }
        const std::string t2ProductFile = directoryT2Product + "dims_" + std::to_string(random_number_dimension) + "_order_" + std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";
        const std::string t3ProductFile = directoryT3Product + "dims_" + std::to_string(random_number_dimension) + "_order_" + std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";

        std::vector<std::vector<double>> tensor2d(No, std::vector<double>(No, 0.0));
        for (int i = 0; i < No; ++i) {
            for (int j = 0; j < No; ++j) {      
                for (int m = 0; m < No; m++) {
                    // t3Product[i][j][m] = 0.0;
                    t3Product[i*No*No + j*No + m] = 0.0;
                }
            }
        }

        // Check if t2Product results exist
        std::cout << "t2product" << std::endl;   
        if (fileExists(t2ProductFile)) {
            std::cout << "Loading t2Product from file." << std::endl;
            readVector1D(t2ProductFile, t2Product);
        } else {
            std::cout << "Calculating t2Product." << std::endl;
            for (int i = 0; i < No; ++i) {
                for (int j = 0; j < No; ++j) {
                    for (int m = 0; m < total_nq; ++m) {
                        tensor2d[i][j] += evaluate(i, points_weights_index_list[m]) * evaluate(j, points_weights_index_list[m]) * weights_multiplied[m];
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

        // print t2Product
        // for (int i = 0; i < No; ++i) {          
        //     std::cout << t2Product[i] << "\t";
        // }
        // std::cout << std::endl;


        std::cout << "t3product" << std::endl;
        // Check if t3Product results exist
        if (fileExists(t3ProductFile)) {
            std::cout << "Loading t3Product from file." << std::endl;
            readVector1D(t3ProductFile, t3Product);
        } else {
            std::cout << "Calculating t3Product." << std::endl;
            for (int i = 0; i < No; ++i) {
                for (int j = 0; j < No; ++j) {
                    for (int k = 0; k < No; k++) {
                        for (int m = 0; m < total_nq; ++m) {
                            double pi_weights = 1.0;
                            for (int n_weights = 0; n_weights < random_number_dimension; ++n_weights) {
                                pi_weights *= weights[n_weights * nq + points_weights_index_list[m][n_weights]];
                            }
                            t3Product[i*No*No + j*No + k] += evaluate(i, points_weights_index_list[m]) * evaluate(j, points_weights_index_list[m]) * evaluate(k, points_weights_index_list[m]) * pi_weights;
                        }
                    }
                }
            }
            saveVector1D(t3ProductFile, t3Product);
        }
    }

    void chaos2ran(const std::vector<double>& chaos, std::vector<double>&Ran)
    {
        std::vector<double> _ran(total_nq, 0.0);

        for (int k = 0; k < total_nq; k++) {
            auto startIt = phiRan.begin() + k * No;
            _ran[k] = std::inner_product(chaos.begin(), chaos.end(), startIt, 0.0);
        }
        Ran.swap(_ran);
    }

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