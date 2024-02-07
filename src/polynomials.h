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

        for (int i = 0; i < No; ++i) {
            for (int j = 0; j < d; ++j) {
                std::cout << inds[i][j] << "\t";
            }
            std::cout << std::endl;
        }
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

        std::cout << "inds" << std::endl;
        for(int i = 0; i < No; ++i) {
            for (int j = 0; j < d; ++j) {
                std::cout << inds[i][j] << "\t";
            }
            std::cout << std::endl;
        }
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
public:
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
    std::vector<std::vector<double>> points;
    std::vector<std::vector<double>> weights;
    std::vector<std::vector<double>> points_tensor;
    std::vector<std::vector<double>> weights_tensor;
    
    std::vector<std::vector<double>> phiRan;
    std::vector<double> t2Product;
    std::vector<std::vector<std::vector<double>>> t3Product;

    polynomials(int _nq, int _order, std::vector<double> _parameter1, std::vector<double> _parameter2, std::vector<int> _polynomial_types, size_t _points_weights_method){
        points_weights_method = _points_weights_method;
        nq = _nq;
        order = _order;
        polynomial_types = _polynomial_types;
        random_number_dimension = polynomial_types.size();
        No = numberPolynomials(random_number_dimension, order);
        inds = calculateMultiIndices(random_number_dimension, order);
        parameter1 = _parameter1;
        parameter2 = _parameter2;

        phiRan.resize(No);
        
        total_nq = std::pow(nq, random_number_dimension);
        
        for (int i = 0; i < No; ++i)
        {
            phiRan[i].resize(total_nq);
            for (int k = 0; k < total_nq; ++k){
                phiRan[i][k] = 0.0;
            }
        }
        
        t2Product.resize(No);
        t3Product.resize(No);
        for (int i = 0; i < No; ++i)
        {
            t3Product[i].resize(No);
            for (int j = 0; j < No; ++j)
            {
                t3Product[i][j].resize(No);
            }
        }

        points.resize(random_number_dimension);
        weights.resize(random_number_dimension);
        polynomial_coefficients.resize(random_number_dimension);
        for (int i = 0; i < random_number_dimension; ++i) {
            points[i].resize(nq);
            weights[i].resize(nq);
            polynomial_coefficients[i].resize(No);
            for (int j = 0; j < No; j++) {
                polynomial_coefficients[i][j].resize(No);
            }
        }

        for (int i = 0; i < random_number_dimension; ++i) {
            std::cout << "polynomial " << i+1 << std::endl;
            polynomial op(nq, order, parameter1[i], parameter2[i], polynomial_types[i], points_weights_method);
            points[i] = op.points;
            weights[i] = op.weights;
            polynomial_coefficients[i] = op.polynomialCoeffs;
        }

        points_tensor.resize(random_number_dimension);
        weights_tensor.resize(random_number_dimension);
        //std::cout << "points" << std::endl;
        for (int i = 0; i < random_number_dimension; ++i) {
            points_tensor[i].resize(total_nq);
            weights_tensor[i].resize(total_nq);
            for(int j = 0; j < total_nq; ++j) {
                std::vector<int> index = find_index(j, random_number_dimension, nq);
                points_tensor[i][j] = points[i][index[i]];
                weights_tensor[i][j] = weights[i][index[i]];
                //std::cout << points_tensor[i][j] << std::endl;
            }
        }

        evaluatePhiRan();
        tensor();

    }

    //n_order at point k
    double evaluate(int n_order, int k){
        double result = 1.0;
        for (int i = 0; i < random_number_dimension; ++i) {
            //std:: cout << "random variable " << i << " at point " << points[i][k] << " max order: " << inds[n_order][i] << " evaluate: " << evaluate(inds[n_order][i], k, i) << std::endl;
            result *= evaluate(inds[n_order][i], k, i);
        }
        return result;
    }
    
    //n_order at point k given polynomial basis 
    double evaluate(int n_order, int k, int phi_i){
        double sum = 0.0;
        for (int j = 0; j < n_order + 1; ++j) {
            sum += polynomial_coefficients[phi_i][n_order][j] * std::pow(points[phi_i][k], j);
        }
        return sum;
    }

    double evaluate(int n_order, std::vector<int> idx){
        double result = 1.0;
        for (int i = 0; i < random_number_dimension; ++i) {
            //std:: cout << "random variable " << i << " at point " << points[i][k] << " max order: " << inds[n_order][i] << " evaluate: " << evaluate(inds[n_order][i], k, i) << std::endl;
            result *= evaluate(inds[n_order][i], points[i][idx[i]], i);
        }
        return result;
    }

    double evaluate(int n_order, double x, int phi_i){
        double sum = 0.0;
        for (int j = 0; j < n_order + 1; ++j) {
            sum += polynomial_coefficients[phi_i][n_order][j] * std::pow(x, j);
        }
        return sum;
    }

    //evaluate the polynomial at kth point 
    double evaluate_polynomial(int order_max, int k)
    {
        double sum = 0.0;
        for (int i = 0; i < order_max+1; ++i){
            sum += evaluate(i, k);
        }
        return sum;
    }

    void evaluatePhiRan() {        
        //std::cout << "phiRan" << std::endl;
        for (int k = 0; k < total_nq; ++k) {
            for (int i = 0; i < No; ++i) {
                std::vector<int> idx = find_index(k, random_number_dimension, nq);
                phiRan[i][k] = evaluate(i, idx);
                //std::cout << phiRan[i][k] << "\t";
            }
            std::cout << std::endl;
        }
    }


    void tensor(){
        std::vector<std::vector<double>> tensor2d(No, std::vector<double>(No, 0.0));
        std::vector<std::vector<std::vector<double>>> tensor2d_all(random_number_dimension, std::vector<std::vector<double>>(order+1, std::vector<double>(order+1, 0.0)));

        for (int i = 0; i < No; ++i) {
            for (int j = 0; j < No; ++j) {      
                for (int m = 0; m < No; m++) {
                    t3Product[i][j][m] = 0.0;
                }
            }
        }

        std::cout << "t2product" << std::endl;
                
        for (int i = 0; i < No; ++i) {
            for (int j = 0; j < No; ++j) {
                for (int m = 0; m < total_nq; ++m) {
                    std::vector<int> index = find_index(m, random_number_dimension, nq);
                    double pi_weights = 1.0;
                    for (int n_weights = 0; n_weights < random_number_dimension; ++n_weights) {
                        pi_weights *= weights[n_weights][index[n_weights]];
                    }
                    tensor2d[i][j] += evaluate(i, index) * evaluate(j, index) * pi_weights;
                }
            }
        }
        
        for (int i = 0; i < No; ++i) {
            t2Product[i] = tensor2d[i][i];            
            //std::cout << t2Product[i] << "\t";
        }
        //std::cout << std::endl;


        std::cout << "t3product" << std::endl;
        for (int i = 0; i < No; ++i) {
            for (int j = 0; j < No; ++j) {
                for (int k = 0; k < No; k++) {
                    for (int m = 0; m < total_nq; ++m) {
                        std::vector<int> index = find_index(m, random_number_dimension, nq);
                        double pi_weights = 1.0;
                        for (int n_weights = 0; n_weights < random_number_dimension; ++n_weights) {
                            pi_weights *= weights[n_weights][index[n_weights]];
                        }
                        t3Product[i][j][k] += evaluate(i, index) * evaluate(j, index) * evaluate(k, index) * pi_weights;
                    }
                }
            }
        }

                
        /*for (int i = 0; i < No; ++i) {  
            for (int j = 0; j < No; ++j) {
                std::cout << t3Product[0][i][j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;*/
    }

    void evaluatePCE(std::vector<double> u, std::vector<double>&uRan)
    {
        
        //for (int i = 0; i < order+1; ++i){
        //    std::cout << u[i] << "\t" ;
        //}
        //std::cout << std::endl;

        std::vector<double> _ran(total_nq, 0.0);

        //for (int i = 0; i < nq; ++i){
        //    std::cout << _ran[i] << "\t" ;
        //}
        //std::cout << std::endl;

        for (int k = 0; k < total_nq; k++) {
            for (int i = 0; i < No; ++i) {
                _ran[k] += u[i] * phiRan[i][k];
            }
        }
        //for (int i = 0; i < nq; ++i){
        //    std::cout << _ran[i] << "\t" ;
        //}
        //std::cout << std::endl;
        uRan.swap(_ran);

        //for (int i = 0; i < nq; ++i){
        //    std::cout << uRan[i] << "\t" ;
        //}
        //std::cout << std::endl;
        //std::cout << std::endl;

        _ran.clear();
    }

    void ran2chaos(std::vector<double> ran, int ran_id, std::vector<double>&chaos)
    {    
        std::vector<double> _chaos(No, 0.0);

        //for (int i = 0; i < order+1; ++i){
        //    std::cout << _chaos[i] << "\t" ;
        //}
        //std::cout << std::endl;

        for (int i = 0; i < No; ++i) {
            for (int k = 0; k < total_nq; ++k) {
                std::vector<int> index = find_index(k, random_number_dimension, nq);
                double pi_weights = 1.0;
                for (int n_weights = 0; n_weights < random_number_dimension; ++n_weights) {
                    pi_weights *= weights[n_weights][index[n_weights]];
                }
                _chaos[i] += pi_weights * ran[k] * phiRan[i][k] / t2Product[i];
            }
        }
        //for (int i = 0; i < order+1; ++i){
        //    std::cout << _chaos[i] << "\t" ;
        //}
        //std::cout << std::endl;
        chaos.swap(_chaos);
        
        //for (int i = 0; i < order+1; ++i){
        //    std::cout << chaos[i] << "\t" ;
        //}
        //std::cout << std::endl;
        //std::cout << std::endl;
        _chaos.clear();
    }

    void chaos_product(std::vector<double> chaos1, std::vector<double> chaos2, std::vector<double>&production)
    {
        for (int i = 0; i < No; ++i) {
            double sum = 0.0;

            for (int j = 0; j < No; ++j) {
                for (int k = 0; k < No; ++k) {
                sum += chaos1[j] * chaos2[k] * t3Product[j][k][i];
                }
            }

            production[i] = sum / t2Product[i];
        }
    }

    void convert2affinePCE(double para1, double para2, int polynomialType, std::vector<double>&domain)
    {   
        if (polynomialType == 0) {
            double a1 = 0.5 * (para1 + para2);
            double a2 = 0.5 * (para2 - para1);
            domain[0] = a1;
            domain[1] = a2;
        }
        else if (polynomialType == 1) {
            double a1 = para1;
            double a2 = para2;
            domain[0] = a1;
            domain[1] = a2;
        }
    }

    double mean(std::vector<double> chaos)
    {
        return chaos[0];
    }

    double std(std::vector<double> chaos)
    {
        double sum = 0.0;
        for (int i = 1; i < order+1; ++i){
            sum += t2Product[i] * chaos[i] * chaos[i];
        }
        return std::sqrt(sum);
    }
};

