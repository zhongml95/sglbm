#pragma once

#include <unistd.h>
#include <limits.h>
#include <libgen.h>    // For dirname
#include <sys/stat.h>  // For stat, mkdir
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <iostream>

// Utility function to check if a directory exists
bool directoryExists(const std::string& path);

// Utility function to create a directory
bool createDirectory(const std::string& path);

// Utility function to delete a directory
bool deleteDirectory(const std::string& path);

// Utility function to find the relative path to the "src" directory
std::string findRelativePathToSrc();

// Structure to hold various parameters
struct Parameters {
    size_t points_weights_method;
    unsigned int order;
    unsigned int nq;
    std::vector<int> polynomialType;
    int resolution;
    std::vector<double> parameter1;
    std::vector<double> parameter2;
    double L;
    double lx;
    double ly;
    double Re;
    double physVelocity;
    double physViscosity;
    double tau;
    double Ma;
};

// Utility function to read parameters from a file
bool readParameters(const std::string& filePath, Parameters& params);

// Utility function to check if a file exists
bool fileExists(const std::string& name);

// Utility function to save a 1D vector to a binary file
void saveVector1D(const std::string& filePath, const std::vector<double>& vec);

// Utility function to read a 1D vector from a binary file
void readVector1D(const std::string& filePath, std::vector<double>& vec);

// Utility function to save a 3D vector to a binary file
void saveVector3D(const std::string& filePath, const std::vector<std::vector<std::vector<double>>>& vec);

// Utility function to read a 3D vector from a binary file
void readVector3D(const std::string& filePath, std::vector<std::vector<std::vector<double>>>& vec);
