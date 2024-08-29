// matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <numeric> // for std::inner_product



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

// Function to get the colomn of a matrix
std::vector<double> getColomn(std::vector<std::vector<double>> matrix, int j) {
    std::vector<double> colomn(matrix[0].size(), 0.0);
    for (int i = 0; i < matrix.size(); ++i) {
        colomn[i] = matrix[i][j];
    }
    return colomn;
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

// Function to compute the Frobenius norm of a matrix
double frobeniusNorm(const std::vector<std::vector<double>>& matrix) {
    double sum = 0.0;
    for (const auto& row : matrix) {
        for (double val : row) {
            sum += val * val;
        }
    }
    return std::sqrt(sum);
}

// Function to check if matrix A is converging to matrix B
bool isConverging(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B, double& former_norm, double tolerance) {
    // Ensure A and B are the same size
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrices must be the same size");
    }

    std::vector<std::vector<double>> diff(A.size(), std::vector<double>(A[0].size()));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            diff[i][j] = A[i][j] - B[i][j];
        }
    }

    double norm = frobeniusNorm(diff);
    bool converge = (std::fabs(norm-former_norm)/former_norm) <= tolerance;
    //std::cout << norm << std::endl;
    former_norm = norm;
    return converge;
}


// Function to perform Householder QR decomposition
void householderQR(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Q, std::vector<std::vector<double>>& R) {
    int m = A.size(); // Number of rows
    int n = (A.empty() ? 0 : A[0].size()); // Number of columns (assuming all rows have the same number of columns)
    //R = copyMatrix(A);
    R = A;
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
}

void Householder(std::vector<std::vector<double>>& J, std::vector<std::vector<double>>& R)
{
    int nq = J.size();
    std::vector<std::vector<double>> Q = generateIdentityMatrix(nq);
    std::vector<std::vector<double>> J_new = copyMatrix(J);
    int iter = 1;
    double norm = 1;

    while (true) {
        householderQR(J, Q, R);
        J_new = matrixMultiplication(R, Q);
        if (iter % 100 == 0) {
            if (isConverging(J, J_new, norm, 1E-13)) {
                std::cout << iter << " iterations needed." << std::endl;
                break;
            }
        }
        J.swap(J_new);
        iter++;
    }
}


// Function to construct the Jacobi matrix for Gaussian-Hermite quadrature
std::vector<std::vector<double>> constructJacobiMatrix_Hermite(int n) {
    if (n <= 1) {
        throw std::invalid_argument("The order (n) must be greater than 1.");
    }

    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    // Calculate the rest of the matrix
    for (int i = 1; i < n; ++i) {
        J[i][i - 1] = std::sqrt(i);
        J[i - 1][i] = J[i][i - 1];
    }

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

void constructJacobiMatrix(int nq, int polynomialType, std::vector<std::vector<double>>& J) {
    if (polynomialType == 0) {
        J = constructJacobiMatrix_Legendre(nq);
    } else if (polynomialType == 1) {
        J = constructJacobiMatrix_Hermite(nq);
    }
}



#endif // MATRIX_H