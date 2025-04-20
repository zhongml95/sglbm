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
#include <cmath>
#include <stdexcept>
#include <vector>
#include <memory>
#include <algorithm>
#include <numeric> // for std::inner_product

// #include "matrix.h"

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

bool readParameters(const std::string& filePath, Parameters& params) {
    std::ifstream file(filePath);
    if (!file) {
        std::cerr << "Unable to open file: " << filePath << std::endl;
        return false;
    }

    std::unordered_map<std::string, std::string> paramMap;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            paramMap[key] = value;
        }
    }

    try {
        params.points_weights_method = std::stoul(paramMap["points_weights_method"]);
        params.order = std::stoul(paramMap["order"]);
        params.nq = std::stoul(paramMap["nq"]);
        params.resolution = std::stoi(paramMap["resolution"]);
        params.L = std::stod(paramMap["L"]);
        params.lx = std::stod(paramMap["lx"]);
        params.ly = std::stod(paramMap["ly"]);
        params.Re = std::stod(paramMap["Re"]);
        params.physVelocity = std::stod(paramMap["physVelocity"]);
        params.physViscosity = std::stod(paramMap["physViscosity"]);
        params.tau = std::stod(paramMap["tau"]);
        params.Ma = std::stod(paramMap["Ma"]);

        auto parseVectorDouble = [](const std::string& s) -> std::vector<double> {
            std::vector<double> result;
            std::istringstream iss(s);
            std::string item;
            while (std::getline(iss, item, ' ')) {
                result.push_back(std::stod(item));
            }
            return result;
        };

        auto parseVectorInt = [](const std::string& s) -> std::vector<int> {
            std::vector<int> result;
            std::istringstream iss(s);
            std::string item;
            while (std::getline(iss, item, ' ')) {
                if (!item.empty()) {
                    result.push_back(std::stoi(item));
                }
            }
            return result;
        };

        params.polynomialType = parseVectorInt(paramMap["polynomialType"]);
        params.parameter1 = parseVectorDouble(paramMap["parameter1"]);
        params.parameter2 = parseVectorDouble(paramMap["parameter2"]);

    } catch (const std::exception& e) {
        std::cerr << "Error parsing parameters: " << e.what() << std::endl;
        return false;
    }

    return true;
}


// Utility function to check if a directory exists
bool directoryExists(const std::string& path);

// Utility function to create a directory
bool createDirectory(const std::string& path);

// Utility function to delete a directory
bool deleteDirectory(const std::string& path);

// Utility function to find the relative path to the "src" directory
std::string findRelativePathToSrc();

// Utility function to check if a file exists
bool fileExists(const std::string& name);

// Utility function to save a 1D vector to a binary file
template<typename T>
void saveVector1D(const std::string& filePath, const std::vector<T>& vec);

// Utility function to read a 1D vector from a binary file
template<typename T>
void readVector1D(const std::string& filePath, std::vector<T>& vec);

// Utility function to save a 3D vector to a binary file
template<typename T>
void saveVector3D(const std::string& filePath, const std::vector<std::vector<std::vector<T>>>& vec);

// Utility function to read a 3D vector from a binary file
template<typename T>
void readVector3D(const std::string& filePath, std::vector<std::vector<std::vector<T>>>& vec);
