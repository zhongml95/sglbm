//#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "util.h"


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
        
        std::cout << "calc phiran" << std::endl;
        evaluatePhiRan();
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

    void evaluatePCE(std::vector<double> uChaos, std::vector<double>&uRan)
    {
        
        std::vector<double> _ran(nq, 0.0);

        for (int k = 0; k < nq; ++k){
            for (int i = 0; i < order+1; ++i){
                _ran[k] += uChaos[i] * phiRan[i][k];
            }
        }
        uRan.swap(_ran);
    }

    void ran2chaos(std::vector<double> ran, std::vector<double>&chaos)
    {    
        std::vector<double> _chaos(order+1, 0.0);

        for (int i = 0; i < order+1; ++i){
            for (int k = 0; k < nq; ++k){
                _chaos[i] += weights[k] * ran[k] * phiRan[i][k] / t2Product[i];
            }
        }
        
        chaos.swap(_chaos);
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
<<<<<<< HEAD
        for (int i = 0; i < nq+1; ++i)
        {
            std::vector<double> coeffs_i(nq+1,0.0);
            if ( i == 0 || i == 1 ) {
                coeffs_i = Legendre_coefficients(i, polynomialCoeffs[i], polynomialCoeffs[i]);
            }
            else {
                  coeffs_i = Legendre_coefficients(i, polynomialCoeffs[i-1], polynomialCoeffs[i-2]);
            }
            //std::cout << i << std::endl;
            //for (int j = 0; j < i+1; ++j)
            for (int j = 0; j < nq+1; ++j)
=======
        for (int i = 0; i < order+1; ++i)
        {
            std::vector<double> coeffs_i(order+1,0.0);
            coeffs_i = Legendre_coefficients(i);
            
            for (int j = 0; j < i+1; ++j)
            {
                //normalize the coefficient with the zeroth term
                polynomialCoeffs[i][j] = coeffs_i[j];// / coeffs_i[i];
                //std::cout << polynomialCoeffs[i][j] << "\t";
            }
            //std::cout << std::endl;
        }
        std::cout << "quadrature points" << std::endl;
        quadrature_rule();
        
        for (int i = 0; i < order+1; ++i)
        {
            std::vector<double> coeffs_i(order+1,0.0);
            coeffs_i = Legendre_coefficients(i);
            
            for (int j = 0; j < i+1; ++j)
>>>>>>> bf2a71be27f7927d2a4b59af8409955525ea971c
            {
                //normalize the coefficient with the zeroth term
                if (j < i+1) {
                    polynomialCoeffs[i][j] = coeffs_i[j];// / coeffs_i[i];
                }
                else {
                     polynomialCoeffs[i][j] = 0.;
                }
                //std::cout << polynomialCoeffs[i][j] << "\t";
            }
            //std::cout << std::endl;
        }
        std::cout << "quadrature points" << std::endl;
        polynomial_points_weights();
        
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
                //std::cout << polynomialCoeffs[n][i] << std::endl;
            }
            return sum;
        }
    }

<<<<<<< HEAD

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
            std::cout << "scipy" << std::endl;
            get_polynomial_quadrature_points_and_weights_from_scipy();
        }

        double weights_sum = 0.0;
=======
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
            points[i] = cos(M_PI * double(2.0 * i + 1.0) / double(2.0 * nq));    
        } 
    } 

    void initguess() 
    {   
        for(int i = 0; i < nq; ++i) 
        {   
            points[i] = -cos(M_PI * double(i + 0.75)/double(nq + 0.5));    
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

    void quadrature_rule()
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

        /*std::vector<std::vector<double>> J;
        if(polynomialType == 0) {
            J = constructJacobiMatrix_Legendre(nq);
        }
        else if (polynomialType == 1) {
            J = constructJacobiMatrix_Hermite(nq);
        }

        std::vector<std::vector<double>> Q = generateIdentityMatrix(nq);
        std::vector<std::vector<double>> R = copyMatrix(J);
        int iter = 1;
        while(iter < 10000){
            householderQR(J, Q, R);
            J = matrixMultiplication(R,Q);
            iter++;
            //std::cout << iter << std::endl;
        }

>>>>>>> bf2a71be27f7927d2a4b59af8409955525ea971c
        std::cout << "points" << std::endl;
        for (int i = 0; i < nq; i++) {
            std::cout <<  points[i] << "\t";
            weights_sum += weights[i];
        }
<<<<<<< HEAD
        std::cout << "sum of weights: " << weights_sum << std::endl;        
=======
        
        std::sort(points.begin(), points.end());
        for (int i = 0; i < nq; i++) {
            std::cout << points[i] << "\t";
        }*/

        
        std::cout << std::endl;
        std::cout << "weights" << std::endl;
        if(polynomialType == 0) {
            for(int i = 0; i < nq; ++i){
                weights[i] = 1.0/((1-points[i]*points[i])*(Lp[i]*Lp[i]));
                //std::cout << std::setprecision(20);
                //std::cout << weights[i] << "\t";
            }
            /*for(int i = 0; i < nq; ++i){
                weights[i] = 1.0  / ( (1.0 - points[i]*points[i]) * std::pow(legendreDerivative(nq, i), 2.0) );
                //std::cout << std::setprecision(20);
                std::cout << weights[i] << "\t";
            }
            std::cout << std::endl;*/
        }
        else if (polynomialType == 1) {        
            double weights_sum = 0.0;
            for (int i = 0; i < nq; i++) {
                weights[i] = std::pow(2, nq-1) * std::tgamma(nq+1) * std::sqrt(M_PI) / pow(nq * evaluate(nq-1, i), 2);
                weights_sum += weights[i];
            }
>>>>>>> bf2a71be27f7927d2a4b59af8409955525ea971c

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

    void get_polynomial_quadrature_points_from_scipy()
    {
        // Get the directory of the executable and append "/src/"
        std::string srcPath = findRelativePathToSrc();

        // Add the src directory to sys.path
        PyObject* sysPath = PySys_GetObject("path");
        PyList_Append(sysPath, PyUnicode_FromString(srcPath.c_str()));

        // Import the util module
        PyObject* moduleName = PyUnicode_FromString("util");
        PyObject* utilModule = PyImport_Import(moduleName);
        
        if (!utilModule) {
            PyErr_Print();
            std::cerr << "Failed to import util module" << std::endl;
            return;
        }

        PyObject* func;
        // Get the function from the module
        if (polynomialType == 0) {
            func = PyObject_GetAttrString(utilModule, "get_legendre_quadrature_points");
        }
        else if (polynomialType == 1) {
            func = PyObject_GetAttrString(utilModule, "get_hermite_quadrature_points");
        }
        
        if (!func || !PyCallable_Check(func)) {
            PyErr_Print();
            std::cerr << "Failed to get Python function" << std::endl;
            Py_DECREF(utilModule);
            return;
        }

    // Prepare arguments and call the function
        PyObject* args = PyTuple_Pack(1, PyLong_FromLong(nq));
        PyObject* resultList = PyObject_CallObject(func, args);
        Py_DECREF(args);
        Py_DECREF(func);
        Py_DECREF(utilModule);

        if (!resultList || !PyList_Check(resultList)) {
            PyErr_Print();
            std::cerr << "Function call failed or did not return a list" << std::endl;
            Py_XDECREF(resultList);
            Py_Finalize();
            return;
        }

        // Extract points from the Python list and store them in outPoints
        Py_ssize_t size = PyList_Size(resultList);
        for (Py_ssize_t i = 0; i < size; ++i) {
            PyObject* item = PyList_GetItem(resultList, i); // Borrowed reference, no need to DECREF
            if (PyFloat_Check(item)) {
                double value = PyFloat_AsDouble(item);
                points[i] = value;
            }
        }

        Py_DECREF(resultList);
    }

    void get_polynomial_quadrature_weights_from_scipy() {
        // Get the directory of the executable and append "/src/"
        std::string srcPath = findRelativePathToSrc();

        std::cout << "srcPath: " << srcPath << std::endl;

        // Add the src directory to sys.path
        PyObject* sysPath = PySys_GetObject("path");
        PyList_Append(sysPath, PyUnicode_FromString(srcPath.c_str()));

        // Import the util module
        PyObject* moduleName = PyUnicode_FromString("util");
        PyObject* utilModule = PyImport_Import(moduleName);
        if (!utilModule) {
            PyErr_Print();
            std::cerr << "Failed to import util module" << std::endl;
            return;
        }

        PyObject* func;
        // Get the function from the module
        if (polynomialType == 0) {
            func = PyObject_GetAttrString(utilModule, "get_legendre_quadrature_weights");
        }
        else if (polynomialType == 1) {
            func = PyObject_GetAttrString(utilModule, "get_Hermite_quadrature_weights");
        }
        
        if (!func || !PyCallable_Check(func)) {
            PyErr_Print();
            std::cerr << "Failed to get Python function" << std::endl;
            Py_DECREF(utilModule);
            return;
        }

// Prepare arguments and call the function
    PyObject* args = PyTuple_Pack(1, PyLong_FromLong(nq));
    PyObject* resultList = PyObject_CallObject(func, args);
    Py_DECREF(args);
    Py_DECREF(func);
    Py_DECREF(utilModule);

    if (!resultList || !PyList_Check(resultList)) {
        PyErr_Print();
        std::cerr << "Function call failed or did not return a list" << std::endl;
        Py_XDECREF(resultList);
        Py_Finalize();
        return;
    }

    // Extract points from the Python list and store them in outPoints
    Py_ssize_t size = PyList_Size(resultList);
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyList_GetItem(resultList, i); // Borrowed reference, no need to DECREF
        if (PyFloat_Check(item)) {
            double value = PyFloat_AsDouble(item);
            weights[i] = value;
        }
    }

    Py_DECREF(resultList);
}
    
    void get_polynomial_quadrature_points_and_weights_from_scipy()
    {        
        Py_Initialize();
        get_polynomial_quadrature_points_from_scipy();
        get_polynomial_quadrature_weights_from_scipy();
        Py_Finalize();
    }

    void convert2affinePCE(double para1, double para2, std::vector<double>&domain)
    {   
        if (polynomialType == 0) {
            double a1 = 0.5 * (para1 + para2);
            double a2 = 0.5 * (para2 - para1);
            domain[0] = a1;// + alpha[0]  * a2;
            //std::cout << "alpha0 = " << alpha[0] << std::endl;
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
