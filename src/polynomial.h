#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iomanip>
#include "matrix.h"

class polynomial
{
public:
    int nq;
    int order;
    int polynomialType;
    double parameter1;
    double parameter2;

    std::vector<double> points;
    std::vector<double> weights;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> t2Product;
    std::vector<std::vector<double>> phi;
    std::vector<std::vector<double>> phiRan;
    std::vector<std::vector<double>> polynomialCoeffs;
    std::vector<std::vector<std::vector<double>>> t3Product;

    polynomial(int _nq, int _n, double _parameter1, double _parameter2, int _polynomialType){
        nq = _nq;
        order = _n;
        parameter1 = _parameter1;
        parameter2 = _parameter2;
        polynomialType = _polynomialType;
        
        points.resize(nq);
        weights.resize(nq);
        phi.resize(order+1);
        alpha.resize(order+1);
        beta.resize(order+1);

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
        
        phiRan.resize(order+1);
        for (int i = 0; i < order+1; ++i)
        {
            phiRan[i].resize(nq);
            for (int k = 0; k < nq; ++k){
                phiRan[i][k] = 0.0;
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
    };

    std::vector<double> getPoints() {
        return points;
    }

    std::vector<double> getWeights() {
        return weights;
    }

    void evaluatePhiRan() {        
        //std::cout << "calc phiRan" << std::endl;
        for (int k = 0; k < nq; ++k) {
            double previousSumPhi = 0.0;
            for (int i = 0; i < order+1; ++i) {
                for (int j = 0; j < order+1; ++j) {
                    phiRan[i][k] += polynomialCoeffs[i][j] * std::pow(points[k], j);
                }
                //std::cout << phiRan[i][k] << "\t";
            }
            //std::cout << std::endl;
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
        
        //std::cout << "t2product" << std::endl;
        for (int i = 0; i < order+1; ++i) {
            t2Product[i] = tensor2d[i][i];
            
            //std::cout << t2Product[i] << "\t" << std::endl;
        }
        //std::cout << std::endl;


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

    void r_jacobi(double a, double b)
    {
        double nu = (b-a) / (a+b+2);
        double mu = 0.0;
        if (a+b+2 > 128){
            mu = std::exp((a+b+1.0)*std::log(2.0) + (std::lgamma(a+1)+std::lgamma(b+1)-std::lgamma(a+b+2)));
        }
        else{
            mu = std::pow(2, (a+b+1)) * ((std::tgamma(a+1)*std::tgamma(b+1)) / std::tgamma(a+b+2));
        }

        if (order == 0){
            alpha[0] = nu;
            beta[0] = mu;
        }
        else{
            alpha[0] = nu;
            std::vector<double> nab(order+1,1.0);
            for (int i = 1; i < order+1; ++i){
                nab[i] = 2*i+a+b;
                alpha[i] = (b*b-a*a) / (nab[i]*(nab[i]+2));
            }

            beta[0] = mu;
            beta[1] = 4*(a+1)*(b+1)/((a+b+2)*(a+b+2)*(a+b+3));
            for (int i = 2; i < order+1; ++i){
                beta[i] = 4*(i+a)*(i+b)*i*(i+a+b)/(nab[i]*nab[i]*(nab[i]+1)*(nab[i]-1));
            }
        }
    }

    void computeSP2()
    {
        //Computes the nregular scalar products aka 2-norms of the orthogonal polynomials
        //Notice that only the values of β of the recurrence coefficients (α,β) are required. The computation is based on equation (1.3.7) 
        //from Gautschi, W. "Orthogonal Polynomials: Computation and Approximation". Whenever there exists an analytic expressions for β, this function should be used.
        
        r_jacobi(0, 0);
        //scale
        beta[0] = beta[0] * 0.5;

        t2Product[0] = beta[0];
        std::cout << "t2Product" << std::endl;
        for (int i = 1; i < order+1; ++i){
            t2Product[i] = t2Product[i-1]*beta[i];
            //std::cout << t2Product[i] << std::endl;
        }
        //std::cout << std::endl;
    }

    void computeSP3()
    {
        for (int i = 0; i < order+1; ++i){
            for (int j = 0; j < order+1; ++j){
                for (int m = 0; m < order+1; m++){
                    for (int k = 0; k < nq; ++k){
                        t3Product[i][j][m] += evaluate(i,k) * evaluate(j,k) * evaluate(m,k) * weights[k];
                    }
                }
            }
        }
    }

    void evaluatePCE(std::vector<double> u, std::vector<double>&uRan)
    {
        
        //for (int i = 0; i < order+1; ++i){
        //    std::cout << u[i] << "\t" ;
        //}
        //std::cout << std::endl;

        std::vector<double> _ran(nq, 0.0);

        //for (int i = 0; i < nq; ++i){
        //    std::cout << _ran[i] << "\t" ;
        //}
        //std::cout << std::endl;

        for (int k = 0; k < nq; ++k){
            for (int i = 0; i < order+1; ++i){
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

    void ran2chaos(std::vector<double> ran, std::vector<double>&chaos)
    {    
        std::vector<double> _chaos(order+1, 0.0);

        //for (int i = 0; i < order+1; ++i){
        //    std::cout << _chaos[i] << "\t" ;
        //}
        //std::cout << std::endl;

        for (int i = 0; i < order+1; ++i){
            for (int k = 0; k < nq; ++k){
                _chaos[i] += weights[k] * ran[k] * phiRan[i][k] / t2Product[i];
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
        std::cout << "construct gpc" << std::endl;
        /*for (int i = 0; i < nq+1; ++i)
        {
            std::vector<double> coeffs_i(nq+1,0.0);
            coeffs_i = Legendre_coefficients(i);
            
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = coeffs_i[j];// / coeffs_i[i];
                //std::cout << polynomialCoeffs[i][j] << "\t";
            }
            //std::cout << std::endl;
        }*/
        std::cout << "quadrature points" << std::endl;
        quadrature_rule();
        
        for (int i = 0; i < nq+1; ++i)
        {
            std::vector<double> coeffs_i(nq+1,0.0);
            coeffs_i = Legendre_coefficients(i);
            
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = coeffs_i[j] / coeffs_i[i];
                //std::cout << polynomialCoeffs[i][j] << "\t";
            }
            //std::cout << std::endl;
        }
        std::cout << "calc phiran" << std::endl;
        evaluatePhiRan();
        std::cout << "calc tensor" << std::endl;
        tensor();
        //computeSP2();
        //computeSP3();
    }

    void HermitePoly()
    {

        for (int i = 0; i < nq; ++i)
        {
            std::vector<double> coeffs_i(nq,0.0);
            coeffs_i = Hermite_coefficients(i);
            
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = coeffs_i[j];
                //std::cout << polynomialCoeffs[i][j] << "\t";
            }
            //std::cout << std::endl;
        }

        quadrature_rule();

        //normalize
        /*for (int i = 0; i < nq; ++i)
        {
            std::vector<double> coeffs_i(nq,0.0);
            coeffs_i = Hermite_coefficients(i);
            
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = coeffs_i[j] / coeffs_i[i];
                std::cout << polynomialCoeffs[i][j] << "\t";
            }
            std::cout << std::endl;
        }*/

        evaluatePhiRan();

        tensor();
        //computeSP2();
        //computeSP3();
        
    }

    // Function to calculate the derivative of Legendre polynomial of order n at x
    double legendreDerivative(int n, int x) {
        if (n == 0) {
            return 0.0;
        } 
        else {
            double sum = 0.0;
            for (int i = 1; i < nq+1; ++i) {
                sum += polynomialCoeffs[n][i] * i * pow(points[x], i-1);
            }
            return sum;
        }
    }

    void quadrature_rule()
    {
        std::vector<std::vector<double>> J;
        if(polynomialType == 0) {
            J = constructJacobiMatrix_Legendre(nq);
        }
        else if (polynomialType == 1) {
            J = constructJacobiMatrix_Hermite(nq);
        }

        std::vector<std::vector<double>> Q = generateIdentityMatrix(nq);
        std::vector<std::vector<double>> R = copyMatrix(J);
        int iter = 1;
        while(iter < 1000){
            householderQR(J, Q, R);
            J = matrixMultiplication(R,Q);
            iter++;
        }

        
        std::cout << "points" << std::endl;
        for (int i = 0; i < nq; i++) {
            points[i] = R[i][i];
        }
        
        std::sort(points.begin(), points.end());
        //for (int i = 0; i < nq; i++) {
        //    std::cout << points[i] << "\t";
        //}

        
        std::cout << std::endl;
        std::cout << "weights" << std::endl;
        if(polynomialType == 0) {
            for(int i = 0; i < nq; ++i){
                weights[i] = 1.0  / ( (1.0 - points[i]*points[i]) * std::pow(legendreDerivative(nq, i), 2.0) );
                //std::cout << std::setprecision(20);
                //std::cout << weights[i] << "\t";
            }
        }
        else if (polynomialType == 1) {        
            double weights_sum = 0.0;
            for (int i = 0; i < nq; i++) {
                weights[i] = std::pow(2, nq-1) * std::tgamma(nq+1) * std::sqrt(M_PI) / pow(nq * evaluate(nq-1, i), 2);
                weights_sum += weights[i];
            }

            //normalize weights
            for (int i = 0; i < nq; i++) {
                weights[i] /= weights_sum;
                //std::cout << weights[i] << "\t";
            }
        }
        //std::cout << std::endl;
    }

    std::vector<double> Legendre_coefficients(int m)
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

        // Consider some form of memoization instead of this recursion
        std::vector<double> v = Legendre_coefficients(m - 1);
        std::vector<double> u = Legendre_coefficients(m - 2);

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
    

    void convert2affinePCE(double para1, double para2, std::vector<double>&domain)
    {   
        if (polynomialType == 0) {
            double a1 = 0.5 * (para1 + para2);
            double a2 = 0.5 * (para2 - para1);
            domain[0] = a1 + alpha[0]  * a2;
            std::cout << "alpha0 = " << alpha[0] << std::endl;
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
