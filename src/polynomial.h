#include <algorithm>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "util.h"
#include <gsl/gsl_integration.h>


class polynomial
{
private:

public:
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
    std::vector<std::vector<double>> phiRan;
    std::vector<std::vector<double>> polynomialCoeffs;
    std::vector<std::vector<std::vector<double>>> t3Product;

    polynomial(int _nq, int _n, double _parameter1, double _parameter2, int _polynomialType, size_t _points_weights_method){
        nq = _nq;
        order = _n;
        parameter1 = _parameter1;
        parameter2 = _parameter2;
        polynomialType = _polynomialType;
        points_weights_method = _points_weights_method;
        
        points.resize(nq);
        weights.resize(nq);
        phi.resize(order+1);
        //alpha.resize(order+1);
        //beta.resize(order+1);


        for (int i = 0; i < order+1; ++i)
        {
            phi[i].resize(nq);
        }

        polynomialCoeffs.resize(nq+1);
        for (int i = 0; i < nq+1; ++i)
        {
            polynomialCoeffs[i].resize(nq+1);
            for (int j = 0; j < nq+1; ++j)
            {
                polynomialCoeffs[i][j] = 0.0;
            }
        }
        
        phiRan.resize(nq);
        for (int k = 0; k < nq; ++k){
        {
            phiRan[k].resize(order+1);
            for (int i = 0; i < order+1; ++i)
                phiRan[k][i] = 0.0;
            }
        }
        
        t2Product.resize(order+1);
        t3Product.resize(order+1);
        for (int i = 0; i < order+1; ++i)
        {
            t3Product[i].resize(order+1);
            for (int j = 0; j < order+1; ++j)
            {
                t3Product[i][j].resize(order+1);
            }
        }

        if(polynomialType == 0) {
            LegendrePoly();
        }
        else if (polynomialType == 1) {
            HermitePoly();
        }

        // print polynomialCoeffs
        // for (int i = 0; i < order+1; ++i) {
        //     for (int j = 0; j < order+1; ++j) {
        //         std::cout << polynomialCoeffs[i][j] << "\t";
        //     }
        //     std::cout << std::endl;
        // }
        
        // std::cout << "calc phiran" << std::endl;
        evaluatePhiRan();

        // print phiran
        // for (int i = 0; i < nq; ++i) {
        //     for (int j = 0; j < order+1; ++j) {
        //         std::cout << phiRan[i][j] << "\t";
        //     }
        //     std::cout << std::endl;
        // }

        std::cout << "calc tensor" << std::endl;
        tensor();
    };

    std::vector<double> getPoints() {
        return points;
    }

    std::vector<double> getWeights() {
        return weights;
    }

    void evaluatePhiRan() {
        for (int k = 0; k < nq; ++k) {
            double previousSumPhi = 0.0;
            for (int i = 0; i < order+1; ++i) {
                for (int j = 0; j < order+1; ++j) {
                    phiRan[k][i] += polynomialCoeffs[i][j] * std::pow(points[k], j);
                }
            }
        }
    }

    double evaluate(int n, int k)
    {       
        double sum = 0.0;
        for (int j = 0; j < n+1; ++j) {
            sum += polynomialCoeffs[n][j] * std::pow(points[k], j);
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

    void tensor(){
        std::vector<std::vector<double>> tensor2d;
        tensor2d.resize(order+1);
        for (int i = 0; i < order+1; ++i) {
            tensor2d[i].resize(order+1);
        }

        for (int i = 0; i < order+1; ++i) {
            for (int j = 0; j < order+1; ++j) {
                tensor2d[i][j] = 0.0;
                for (int m = 0; m < order+1; m++) {
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
        
        // print tensor2d
        // std::cout << "t2product" << std::endl;
        // for (int i = 0; i < order+1; ++i) {
        //     t2Product[i] = tensor2d[i][i];
            
        //     std::cout << t2Product[i] << "\t" << std::endl;
        // }
        // std::cout << std::endl;


        for (int i = 0; i < order+1; ++i) {
            for (int j = 0; j < order+1; ++j) {
                for (int m = 0; m < order+1; m++) {
                    for (int k = 0; k < nq; ++k) {
                        t3Product[i][j][m] += evaluate(i,k) * evaluate(j,k) * evaluate(m,k) * weights[k];
                    }
                }
            }
        }
    }

    void chaos2ran(const std::vector<double>& uChaos, std::vector<double>&uRan)
    {
        
        std::vector<double> _ran(nq, 0.0);

        for (int k = 0; k < nq; ++k){
            _ran[k] = std::inner_product(uChaos.begin(), uChaos.end(), phiRan[k].begin(), 0.0);
        }
        uRan.swap(_ran);
    }

    void ran2chaos(const std::vector<double>& ran, std::vector<double>&chaos)
    {    
        std::vector<double> _chaos(order+1, 0.0);

        for (int i = 0; i < order+1; ++i){
            for (int k = 0; k < nq; ++k){
                _chaos[i] += weights[k] * ran[k] * phiRan[i][k] / t2Product[i];
            }
        }
        
        chaos.swap(_chaos);
    }

    void chaos_product(const std::vector<double>& chaos1, const std::vector<double>& chaos2, std::vector<double>&production)
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

    void LegendrePoly()
    {

        //normalized the coefficients
        std::cout << "construct gpc (Legendre polynomial basis)" << std::endl;
        for (int i = 0; i < nq+1; ++i)
        {
            std::vector<double> coeffs_i(nq+1,0.0);
            if ( i == 0 || i == 1 ) {
                coeffs_i = Legendre_coefficients(i, polynomialCoeffs[i], polynomialCoeffs[i]);
            }
            else {
                  coeffs_i = Legendre_coefficients(i, polynomialCoeffs[i-1], polynomialCoeffs[i-2]);
            }
            for (int j = 0; j < nq+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                if (j < i+1) {
                    polynomialCoeffs[i][j] = coeffs_i[j];
                }
                else {
                     polynomialCoeffs[i][j] = 0.;
                }
            }
        }

        // normalize the coefficients
        for (int i = 0; i < nq+1; ++i)
        {
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = polynomialCoeffs[i][j] / polynomialCoeffs[i][i];
            }
        }

        std::cout << "quadrature points" << std::endl;
        polynomial_points_weights();
        
    }

    void HermitePoly()
    {
        std::cout << "construct gpc (Hermite polynomial basis)" << std::endl;
        for (int i = 0; i < nq; ++i)
        {
            std::vector<double> coeffs_i(nq,0.0);
            coeffs_i = Hermite_coefficients(i);
            
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = coeffs_i[j] / coeffs_i[i];
            }
        }

        // normalize the coefficients
        for (int i = 0; i < nq; ++i)
        {
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = polynomialCoeffs[i][j] / polynomialCoeffs[i][i];
            }
        }

        polynomial_points_weights();
        
    }

    // Function to calculate the derivative of Legendre polynomial of order n at x
    double legendreDerivative(int n, int x) {
        if (n == 0) {
            return 0.0;
        } 
        else {
            double sum = 0.0;
            for (int i = 1; i < order+1; ++i) {
                sum += polynomialCoeffs[n][i] * i * pow(points[x], i-1);
            }
            return sum;
        }
    }


    void polynomial_points_weights()
    {
        if (points_weights_method == 0) {
            std::cout << "householder" << std::endl;
            std::vector<std::vector<double>> J;
            if(polynomialType == 0) {
                J = constructJacobiMatrix_Legendre(nq);
            }
            else if (polynomialType == 1) {
                J = constructJacobiMatrix_Hermite(nq);
            }

            std::vector<std::vector<double>> Q = generateIdentityMatrix(nq);
            std::vector<std::vector<double>> R = J;
            std::vector<std::vector<double>> J_new = copyMatrix(J);
            int iter = 1;
            double norm = 1;

            while(1){
                householderQR(J, Q, R);
                J_new = matrixMultiplication(R,Q);
                if (iter % 100 == 0) {
                    if (isConverging(J, J_new, norm, 1E-13)){
                        std::cout << iter << " is needed." << std::endl;
                        break;
                    }
                }
                J.swap(J_new);
                iter++;
            }

            for (int i = 0; i < nq; i++) {
                points[i] = R[i][i];
            }
            std::sort(points.begin(), points.end());
            
            if(polynomialType == 0) {
                for(int i = 0; i < nq; ++i){
                    weights[i] = 1.0  / ( (1.0 - points[i]*points[i]) * std::pow(legendreDerivative(nq, i), 2.0) );
                }
            }
            else if (polynomialType == 1) {        
                for (int i = 0; i < nq; i++) {
                    weights[i] = std::pow(2, nq-1) * std::tgamma(nq+1) * std::sqrt(M_PI) / pow(nq * evaluate(nq-1, i), 2);
                }
            }
        }
        else if (points_weights_method == 1) {
            std::cout << "GSL" << std::endl;
            points = get_gauss_legendre_points(nq);
            weights = get_gauss_legendre_weights(nq);
        }

        double weights_sum = 0.0;
        std::cout << "points" << std::endl;
        for (int i = 0; i < nq; i++) {
            std::cout <<  points[i] << "\t";
            weights_sum += weights[i];
        }
        std::cout << std::endl;
        std::cout << "weights" << std::endl;
        for (int i = 0; i < nq; i++) {
            std::cout <<  weights[i] << "\t";
        }
        std::cout << std::endl;
        std::cout << "sum of weights: " << weights_sum << std::endl;        

    }

    std::vector<double> Legendre_coefficients(int m, std::vector<double> v, std::vector<double> u)
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
    std::vector<double>  Hermite_coefficients(int m) {

        std::vector<double> coeffs(m + 1, 0.0);
        if (m == 0) {               
            return {1};
        }
        else if (m == 1) {
            return {0, 1};
        }
        else {                
            std::vector<double> m_1 = Hermite_coefficients(m - 1);
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

    // Function to return Gauss-Legendre quadrature points
    std::vector<double> get_gauss_legendre_points(size_t n) {
        std::vector<double> points;
        gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(n);

        double xi, wi;
        for (size_t i = 0; i < n; ++i) {
            gsl_integration_glfixed_point(-1.0, 1.0, i, &xi, &wi, table);
            points.push_back(xi);
        }

        gsl_integration_glfixed_table_free(table);
        return points;
    }

    // Function to return Gauss-Legendre quadrature weights
    std::vector<double> get_gauss_legendre_weights(size_t n) {
        std::vector<double> weights;
        gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(n);

        double xi, wi;
        for (size_t i = 0; i < n; ++i) {
            gsl_integration_glfixed_point(-1.0, 1.0, i, &xi, &wi, table);
            weights.push_back(wi * 0.5);
        }

        gsl_integration_glfixed_table_free(table);
        return weights;
    }

    void convert2affinePCE(double para1, double para2, std::vector<double>&domain)
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
