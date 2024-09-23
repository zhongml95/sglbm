// matrix.cpp
#include "matrix.h"

// Function implementations

double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be the same size for dot product.");
    }
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

// Function to subtract one vector from another
std::vector<double> vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be the same size for subtraction.");
    }
    std::vector<double> result(a.size());
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
    if (norm == 0.0) {
        throw std::runtime_error("Cannot normalize a zero vector.");
    }
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] /= norm;
    }
}

// Function to get the column of a matrix
std::vector<double> getColumn(const std::vector<std::vector<double>>& matrix, int j) {
    if (matrix.empty() || j < 0 || j >= matrix[0].size()) {
        throw std::out_of_range("Invalid column index.");
    }
    std::vector<double> column(matrix.size());
    for (size_t i = 0; i < matrix.size(); ++i) {
        column[i] = matrix[i][j];
    }
    return column;
}

std::vector<std::vector<double>> generateIdentityMatrix(int n) {
    if (n <= 0) {
        throw std::invalid_argument("Size of identity matrix must be positive.");
    }
    std::vector<std::vector<double>> identityMatrix(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        identityMatrix[i][i] = 1.0;
    }
    return identityMatrix;
}

std::vector<std::vector<double>> copyMatrix(const std::vector<std::vector<double>>& A) {
    if (A.empty()) {
        return {};
    }
    int numRows = A.size();
    int numCols = A[0].size();
    std::vector<std::vector<double>> R(numRows, std::vector<double>(numCols));
    for (int i = 0; i < numRows; ++i) {
        if (A[i].size() != numCols) {
            throw std::invalid_argument("All rows must have the same number of columns.");
        }
        for (int j = 0; j < numCols; ++j) {
            R[i][j] = A[i][j];
        }
    }
    return R;
}

// Function to add one vector to another
std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be the same size for addition.");
    }
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Function to multiply two matrices
std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>>& A,
                                                      const std::vector<std::vector<double>>& B) {
    if (A.empty() || B.empty()) {
        throw std::invalid_argument("Input matrices cannot be empty.");
    }
    int m = A.size();
    int n = A[0].size();
    int p = B[0].size();

    // Check dimensions
    if (n != B.size()) {
        throw std::invalid_argument("Number of columns in A must equal number of rows in B.");
    }

    // Initialize result matrix
    std::vector<std::vector<double>> C(m, std::vector<double>(p, 0.0));

    // Perform multiplication
    for (int i = 0; i < m; ++i) {
        if (A[i].size() != n) {
            throw std::invalid_argument("All rows in A must have the same number of columns.");
        }
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
    if (A.empty()) {
        return {};
    }
    int m = A.size();
    int n = A[0].size();
    std::vector<std::vector<double>> result(n, std::vector<double>(m));
    for (int i = 0; i < m; ++i) {
        if (A[i].size() != n) {
            throw std::invalid_argument("All rows in A must have the same number of columns.");
        }
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
bool isConverging(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B,
                  double& former_norm, double tolerance) {
    // Ensure A and B are the same size
    if (A.size() != B.size() || A.empty() || B.empty() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrices must be the same size and non-empty.");
    }
    std::vector<std::vector<double>> diff(A.size(), std::vector<double>(A[0].size()));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            diff[i][j] = A[i][j] - B[i][j];
        }
    }
    double norm = frobeniusNorm(diff);

    // Handle the first iteration
    if (former_norm == 0.0) {
        former_norm = norm;
        return false;
    }

    double relative_change = std::fabs(norm - former_norm) / former_norm;

    bool converged = relative_change <= tolerance;
    former_norm = norm;
    return converged;
}

// Function to perform Householder QR decomposition
void householderQR(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Q,
                   std::vector<std::vector<double>>& R) {
    int m = A.size(); // Number of rows
    int n = (A.empty() ? 0 : A[0].size()); // Number of columns
    R = A;
    Q = generateIdentityMatrix(m);
    for (int k = 0; k < n; ++k) {
        std::vector<double> x(m - k);
        for (int i = k; i < m; ++i) {
            x[i - k] = R[i][k];
        }

        double norm_x = std::sqrt(dotProduct(x, x));
        if (norm_x == 0.0) {
            continue; // Skip the transformation
        }

        double sign = (x[0] < 0) ? -1.0 : 1.0;
        x[0] += sign * norm_x;
        normalize(x);

        // Form the Householder matrix H = I - 2*v*v^T
        std::vector<std::vector<double>> H = generateIdentityMatrix(m - k);
        for (size_t i = 0; i < x.size(); ++i) {
            for (size_t j = 0; j < x.size(); ++j) {
                H[i][j] -= 2 * x[i] * x[j];
            }
        }

        // Extend H to full size
        std::vector<std::vector<double>> H_full = generateIdentityMatrix(m);
        for (size_t i = k; i < m; ++i) {
            for (size_t j = k; j < m; ++j) {
                H_full[i][j] = H[i - k][j - k];
            }
        }

        // Update R and Q
        R = matrixMultiplication(H_full, R);
        Q = matrixMultiplication(Q, transposeMatrix(H_full));
    }
}

void Householder(std::vector<std::vector<double>>& J, std::vector<std::vector<double>>& R) {
    int nq = J.size();
    std::vector<std::vector<double>> Q;
    std::vector<std::vector<double>> J_new = copyMatrix(J);
    int iter = 1;
    double norm = 0.0;

    const int max_iter = 10000; // Set a reasonable maximum
    const double tolerance = 1E-13;

    while (iter <= max_iter) {
        householderQR(J_new, Q, R);
        J_new = matrixMultiplication(R, Q);
        if (iter % 100 == 0) {
            if (isConverging(J, J_new, norm, tolerance)) {
                std::cout << iter << " iterations needed." << std::endl;
                break;
            }
            J.swap(J_new);
        }
        iter++;
    }
    if (iter > max_iter) {
        std::cerr << "Householder method did not converge within the maximum number of iterations." << std::endl;
    }
}

// Function to construct the Jacobi matrix for Gaussian-Hermite quadrature
std::vector<std::vector<double>> constructJacobiMatrix_Hermite(int n) {
    if (n <= 1) {
        throw std::invalid_argument("The order (n) must be greater than 1.");
    }

    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    for (int i = 1; i < n; ++i) {
        J[i][i - 1] = std::sqrt(i / 2.0);
        J[i - 1][i] = J[i][i - 1];
    }

    return J;
}

// Function to construct the Jacobi matrix for Gaussian-Legendre quadrature
std::vector<std::vector<double>> constructJacobiMatrix_Legendre(int n) {
    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));

    for (int i = 1; i < n; ++i) {
        double a = i / std::sqrt(4 * i * i - 1);
        J[i][i - 1] = a;
        J[i - 1][i] = a;
    }

    return J;
}

void constructJacobiMatrix(int nq, int polynomialType, std::vector<std::vector<double>>& J) {
    if (polynomialType == 0) {
        J = constructJacobiMatrix_Legendre(nq);
    } else if (polynomialType == 1) {
        J = constructJacobiMatrix_Hermite(nq);
    } else {
        throw std::invalid_argument("Unsupported polynomial type.");
    }
}
