#include <iostream>
#include <cmath>
#include <vector>

class LegendrePoly
{
public:
    int nq;
    int order;
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

    LegendrePoly(int _nq, int _n, double _parameter1, double _parameter2)
    {
        nq = _nq;
        order = _n;
        parameter1 = _parameter1;
        parameter2 = _parameter2;
        points.resize(nq);
        weights.resize(nq);
        phi.resize(order+1);
        alpha.resize(order+1);
        beta.resize(order+1);
        polynomialCoeffs.resize(order+1);
        for (int i = 0; i < order+1; ++i)
        {
            phi[i].resize(nq);
            polynomialCoeffs[i].resize(order+1);
            for (int j = 0; j < order+1; ++j)
            {
                polynomialCoeffs[i][j] = 0.0;
            }
        }
        
        gauss();

        //std::cout << "coefficients" << std::endl;

        for (int i = 0; i < order+1; ++i)
        {
            std::vector<double> coeffs_i(order+1,0.0);
            coeffs_i = Legendre_coefficients(i);
            
            for (int j = 0; j < order+1; ++j)
            {
                
                polynomialCoeffs[i][j] = coeffs_i[j] / coeffs_i[i];
                //std::cout << polynomialCoeffs[i][j] << "\t";
            }
            //std::cout << std::endl;
        }

        phiRan.resize(order+1);
        for (int i = 0; i < order+1; ++i)
        {
            phiRan[i].resize(nq);
            for (int k = 0; k < nq; ++k){
                phiRan[i][k] = 0.0;
            }
        }

        evaluatePhiRan();
        

        r_jacobi(0, 0);
        //scale
        beta[0] = beta[0] * 0.5;

        /*std::cout <<"alpha\tbeta" << std::endl;
        for (int i = 0; i < order+1; ++i)
        {
            std::cout << alpha[i] << "\t" << beta[i] << std::endl;
        }*/

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

        computeSP2();
        computeSP3();

        /*std::cout << "t2Product" << std::endl;
        for (int i = 0; i < order+1; ++i)
        {
            std::cout << t2Product[i] << "\t";
        }
        std::cout<<std::endl;

        std::cout << "t3Product, m = 1" << std::endl;
        for (int i = 0; i < order+1; ++i){
            for (int j = 0; j < order+1; ++j){
                std::cout << t3Product[i][j][1] << "\t";
            }
            std::cout<<std::endl;
        }*/
    }

    double eval(int n, double x)
    {   
        if (n == 0)
            return 1;
        else if (n == 1)
            return x;
        else
            return ((2.0 * n - 1.0) * x * eval(n - 1, x) - (n - 1) * eval(n - 2, x)) / n;
            //return ((2.0 * n + 1.0) * x * eval(n, x) - n * eval(n -1, x)) / (n + 1);
    }

    // Compute maximum pointwise difference
    double dist_max(const std::vector<double>& a, const std::vector<double>& b)
    {
        int n = a.size();
        double dist = 0.0;
        double max_dist = 0.0;

        for(int i=0; i<n; ++i) {
            dist = std::abs(a[i]-b[i]);
            if(dist>max_dist) {
                max_dist = dist; 
            }  
        }
        return max_dist;
    }

    // Compute Chebyshev Gauss Nodes
    void chebnodes() 
    {
        for(int i = 0; i < nq; ++i) 
        {
            points[i] = cos(M_PI * double(2 * i + 1) / double(2 * nq));    
        } 
    } 

    void initguess() 
    {   
        for(int i = 0; i < nq; ++i) 
        {   
            points[i] = -cos(M_PI * double(i + .75)/double(nq + 0.5));    
        }     
    }

    void rec_legendre(std::vector<double>& a, std::vector<double>& b)
    {
        int n = a.size();
        for(int i=0;i<n;++i) 
        {
            a[i] = (2.0*i+1.0)/(i+1.0);
            b[i] = double(i)/(i+1.0);
        }
    }

    // Evaluate the nth order Legendre polynomial and its derivative
    void legendre(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& L, std::vector<double>& Lp) 
    {
        std::vector<double> L0(nq,1.0);
        std::vector<double> L1(nq,0.0);
        
        // Iterate over grid points
        for(int j = 0; j < nq; ++j) 
        {
            L1[j] = points[j];
            // Iterate over polynomials  
            for(int k = 1; k < nq; ++k) 
            {
                L[j] = a[k] * points[j] * L1[j] - b[k] * L0[j];
                L0[j] = L1[j];
                L1[j] = L[j];
            }
            Lp[j] = nq * (L0[j] - points[j] * L[j]) / (1.0 - points[j] * points[j]); 
        } 
    } 

    // Update grid points as well as the nth Legendre polynomial and its derivative  
    void newton_step(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& L, std::vector<double>& Lp) 
    {
        legendre(a,b,L,Lp);
        for(int i = 0; i < nq; ++i) 
        {
            points[i] -= L[i]/Lp[i];
        } 
    }

    int gauss() 
    {
        double dist = 1;
        double tol = 1e-15;
        int iter = 0;

        // Use Chebyshev-Gauss nodes and initial guess
        
        //initguess();
        chebnodes();

        std::vector<double> x0(points);
        std::vector<double> L(nq,0.0);
        std::vector<double> Lp(nq,0.0);
        std::vector<double> a(nq,0.0);
        std::vector<double> b(nq,0.0);


        rec_legendre(a,b);

        // Iteratively correct nodes with Newton's method
        while(dist>tol) {     
            newton_step(a,b,L,Lp);
            dist = dist_max(points, x0); 
            ++iter;
            x0 = points;
        } 

        // Compute Weights
        for(int i = 0; i < nq; ++i){
            weights[i] = 1.0/((1-points[i]*points[i])*(Lp[i]*Lp[i]));
            //std::cout << points[i] << std::endl;
            //points[i] = points[i];
        }

        return iter;
    }

    std::vector<double> getPoints()
    {
        return points;
    }

    std::vector<double> getWeights()
    {
        return weights;
    }

    void evaluatePhiRan()
    {        
        //std::cout << "phiRan" << std::endl;
        for (int k = 0; k < nq; ++k)
        {
            double previousSumPhi = 0.0;
            for (int i = 0; i < order+1; ++i)
            {
                for (int j = 0; j < order+1; ++j)
                {
                    phiRan[i][k] += polynomialCoeffs[i][j] * std::pow(points[k], j);
                    
                    //std::cout << points[k] << "\t";
                }
                //std::cout << std::endl;
                //std::cout << phiRan[i][k] << "\t";
            }
            //std::cout << std::endl;
        }
    }

    double evaluate(int n, int k)
    {       
        double sum = 0.0;
        for (int j = 0; j < n+1; ++j){
            sum += polynomialCoeffs[n][j] * std::pow(points[k], j);
        }
        return sum;
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
        std::vector<double> coeffs(m + 1);

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
        t2Product[0] = beta[0];
        for (int i = 1; i < order+1; ++i){
            t2Product[i] = t2Product[i-1]*beta[i];
        }
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

    void tensor()
    {
        std::vector<std::vector<double>> tensor2d;
        std::vector<std::vector<std::vector<double>>> tensor3d;
        tensor2d.resize(order+1);
        tensor3d.resize(order+1);
        for (int i = 0; i < order+1; ++i){
            tensor2d[i].resize(order+1);
            tensor3d[i].resize(order+1);
            for (int j = 0; j < order+1; ++j){
                tensor3d[i][j].resize(order+1);
            }
        }


        for (int i = 0; i < order+1; ++i){
            for (int j = 0; j < order+1; ++j){
                tensor2d[i][j] = 0.0;
                for (int m = 0; m < order+1; m++){
                    tensor3d[i][j][m] = 0.0;
                }
            }
        }

        for (int i = 0; i < order+1; ++i){
            for (int j = 0; j < order+1; ++j){
                for (int k = 0; k < nq; ++k){
                    tensor2d[i][j] += evaluate(i,k) * evaluate(j,k) * weights[k];
                }
            }
        }

        for (int i = 0; i < order+1; ++i){
            for (int j = 0; j < order+1; ++j){
                for (int m = 0; m < order+1; m++){
                    for (int k = 0; k < nq; ++k){
                        tensor3d[i][j][m] += evaluate(i,k) * evaluate(j,k) * evaluate(m,k) * weights[k];
                    }
                }
            }
        }

        /*std::cout << "t2Product2" << std::endl;
        for (int i = 0; i < order+1; ++i){
            for (int j = 0; j < order+1; ++j){
                std::cout << tensor2d[i][j] << "\t";
            }
            std::cout<<std::endl;
        }

        std::cout << "t3Product, m = 2" << std::endl;
        for (int i = 0; i < order+1; ++i){
            for (int j = 0; j < order+1; ++j){
                std::cout << tensor3d[i][j][2] << "\t";
            }
            std::cout<<std::endl;
        }*/
    }

    void convert2affinePCE(double para1, double para2, std::vector<double>&domain)
    {
        double a1 = 0.5 * (para1 + para2);
        double a2 = 0.5 * (para2 - para1);
        domain[0] = a1 + alpha[0]  * a2;
        domain[1] = a2;
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