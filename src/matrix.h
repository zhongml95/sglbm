#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>




double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

// Function to subtract one vector from another
std::vector<double> vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size(), 0.0);
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Function to multiply a vector by a scalar
std::vector<double> vectorScalarProduct(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size(), 0.0);
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] * scalar;
    }
    return result;
}

// Function to normalize a vector
void normalize(std::vector<double>& v) {
    double norm = std::sqrt(dotProduct(v, v));
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] /= norm;
    }
}

// Function to perform Gram-Schmidt orthogonalization
void gramSchmidt(std::vector<std::vector<double>>& vectors) {
    size_t n = vectors.size();
    
    // Initialize the orthogonalized vectors
    std::vector<std::vector<double>> orthogonalized(n, std::vector<double>(vectors[0].size(), 0.0));
    
    for (size_t i = 0; i < n; ++i) {
        orthogonalized[i] = vectors[i];
        for (size_t j = 0; j < i; ++j) {
            double projection = dotProduct(vectors[i], orthogonalized[j]) / dotProduct(orthogonalized[j], orthogonalized[j]);
            orthogonalized[i] = vectorSubtraction(orthogonalized[i], vectorScalarProduct(orthogonalized[j], projection));
        }
        normalize(orthogonalized[i]);
    }
    
    // Print the orthogonalized vectors
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < vectors[i].size(); ++j) {
            vectors[i][j] = orthogonalized[i][j];
        }
    }
}



std::vector<std::vector<double>> generateIdentityMatrix(int n) {
    std::vector<std::vector<double>> identityMatrix(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        identityMatrix[i][i] = 1.0;
    }

    return identityMatrix;
}

std::vector<std::vector<double>> copyMatrix(const std::vector<std::vector<double>>& A) {
    int numRows = A.size();
    int numCols = A[0].size();
    
    // Initialize a new matrix R with the same dimensions as A
    std::vector<std::vector<double>> R(numRows, std::vector<double>(numCols, 0.0));
    
    // Copy values from A to R
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            R[i][j] = A[i][j];
        }
    }
    
    return R;
}


// Function to add one vector to another
std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size(), 0.0);
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}


// Function to multiply two matrices
std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    int m = A.size();
    int n = A[0].size();
    int p = B[0].size();

    std::vector<std::vector<double>> C(m, std::vector<double>(p, 0.0));

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}


// Function to transpose a matrix
std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& A) {
    int m = A.size();
    int n = A[0].size();

    std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            result[j][i] = A[i][j];
        }
    }

    return result;
}





// Function to perform Householder QR decomposition
void householderQR(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Q, std::vector<std::vector<double>>& R) {
    int m = A.size(); // Number of rows
    int n = (A.empty() ? 0 : A[0].size()); // Number of columns (assuming all rows have the same number of columns)
    R = copyMatrix(A);
    Q = generateIdentityMatrix(m);
    for (int k = 0; k < n; ++k) {
        std::vector<double> x(m - k, 0.0);
        for (int i = k; i < m; ++i) {
            x[i - k] = R[i][k];
        }

        double sign = 1.0;
        if(x[0] < 0) {
            sign = -1.0;
        }            
        double norm = std::sqrt(dotProduct(x, x));

        std::vector<double> v = vectorScalarProduct(x, 1.0 / (x[0] + sign * norm));
        v[0] = 1;        
        double tau = 2.0 / dotProduct(v, v);
        
        std::vector<std::vector<double>> H = generateIdentityMatrix(m);

        for (int i = k; i < m; ++i) {
            for (int j = k; j < m; ++j) {
                H[i][j] -= tau * v[i - k] * v[j - k];
                //std::cout << v[i - k] << std::endl;
            }
        }

        // Update R with Householder transformation
        R = matrixMultiplication(H, R);

        // Update Q with Householder transformation
        Q = matrixMultiplication(Q, transposeMatrix(H));
    }

    //std::vector<std::vector<double>> H_1 = generateIdentityMatrix(m);
    //H_1[m-1][m-1] = -1;
    
    // Update R with Householder transformation
    //R = matrixMultiplication(H_1, R);

        // Update Q with Householder transformation
    //Q = matrixMultiplication(Q, transposeMatrix(H_1));
}

// Function to construct the Jacobi matrix for Gaussian-Hermite quadrature
std::vector<std::vector<double>> constructJacobiMatrix_Hermite(int n) {
if (n <= 1) {
        throw std::invalid_argument("The order (n) must be greater than 1.");
    }

    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    // Set the values for the first two diagonals
    //J[0][1] = std::sqrt(0.5);
    //J[1][0] = std::sqrt(0.5);

    // Calculate the rest of the matrix
    for (int i = 1; i < n; ++i) {
        //if (i == 1) {
        //    J[i][i - 1] = std::sqrt(std::sqrt(M_PI));
        //}
        //else {
            J[i][i - 1] = std::sqrt(i);
        //}
        J[i - 1][i] = J[i][i - 1];
    }

    /*double mu = 0.0;

    for (int i = 1; i < n; ++i) {
        if (i == 1) {
            J[i][i - 1] = std::tgamma(mu+0.5);
            std::cout << std::tgamma(mu+0.5) << std::endl;
        }
        else if ((i-1) % 2 == 0) {
            J[i][i - 1] = 0.5 * (i-1);
        }
        else if ((i-1) % 2 != 0) {
            J[i][i - 1] = 0.5 * (i-1) + mu;
        }

        J[i - 1][i] = J[i][i - 1];
    }*/

    /*for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << J[i][j] << "\t";
        }
        std::cout << std::endl;
    }*/

    return J;
}

// Function to construct the Jacobi matrix for Gaussian-Legendre quadrature
std::vector<std::vector<double>> constructJacobiMatrix_Legendre(int n) {
    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));
    
    for (int i = 1; i < n; ++i) {
        J[i][i - 1] = i / std::sqrt(4*i*i-1);
        J[i - 1][i] = J[i][i - 1];
    }


    return J;
}