// matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>

// Function declarations

double dotProduct(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> vectorScalarProduct(const std::vector<double>& v, double scalar);

void normalize(std::vector<double>& v);

std::vector<double> getColumn(const std::vector<std::vector<double>>& matrix, int j);

std::vector<std::vector<double>> generateIdentityMatrix(int n);

std::vector<std::vector<double>> copyMatrix(const std::vector<std::vector<double>>& A);

std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b);

std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>>& A,
                                                      const std::vector<std::vector<double>>& B);

std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& A);

double frobeniusNorm(const std::vector<std::vector<double>>& matrix);

bool isConverging(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B,
                  double& former_norm, double tolerance);

void householderQR(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& Q,
                   std::vector<std::vector<double>>& R);

void Householder(std::vector<std::vector<double>>& J, std::vector<std::vector<double>>& R);

std::vector<std::vector<double>> constructJacobiMatrix_Hermite(int n);

std::vector<std::vector<double>> constructJacobiMatrix_Legendre(int n);

void constructJacobiMatrix(int nq, int polynomialType, std::vector<std::vector<double>>& J);

#endif // MATRIX_H
