#pragma once
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <stdexcept>

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

std::string trim(const std::string& str) {
    const auto begin = str.find_first_not_of(" \t\n\r");
    if (begin == std::string::npos) return "";
    const auto end = str.find_last_not_of(" \t\n\r");
    return str.substr(begin, end - begin + 1);
}

std::unordered_map<std::string, std::string> parseSimpleXML(const std::string& filePath) {
    std::unordered_map<std::string, std::string> paramMap;
    std::ifstream file(filePath);
    std::string line;

    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] != '<' || line[1] == '/') continue;

        size_t tagStart = line.find('<') + 1;
        size_t tagEnd = line.find('>', tagStart);
        std::string key = line.substr(tagStart, tagEnd - tagStart);

        size_t valueStart = tagEnd + 1;
        size_t valueEnd = line.find('<', valueStart);
        std::string value = line.substr(valueStart, valueEnd - valueStart);

        paramMap[key] = trim(value);
    }

    return paramMap;
}

bool readParametersFromXML(const std::string& filePath, Parameters& params) {
    try {
        auto map = parseSimpleXML(filePath);

        // auto toInt = [](const std::string& s) { return std::stoi(s); };
        // auto toUInt = [](const std::string& s) { return static_cast<unsigned int>(std::stoul(s)); };
        // auto toDouble = [](const std::string& s) { return std::stod(s); };

        auto parseVectorDouble = [](const std::string& s) {
            std::vector<double> vec;
            std::istringstream iss(s);
            std::string item;
            while (iss >> item) vec.push_back(std::stod(item));
            return vec;
        };

        auto parseVectorInt = [](const std::string& s) {
            std::vector<int> vec;
            std::istringstream iss(s);
            std::string item;
            while (iss >> item) vec.push_back(std::stoi(item));
            return vec;
        };

        params.points_weights_method = std::stoul(map.at("points_weights_method"));
        params.order                 = std::stoul(map.at("order"));
        params.nq                    = std::stoul(map.at("nq"));
        params.resolution            = std::stoi(map.at("resolution"));
        params.L                     = std::stod(map.at("L"));
        params.lx                    = std::stod(map.at("lx"));
        params.ly                    = std::stod(map.at("ly"));
        params.Re                    = std::stod(map.at("Re"));
        params.physVelocity          = std::stod(map.at("physVelocity"));
        params.physViscosity         = std::stod(map.at("physViscosity"));
        params.tau                   = std::stod(map.at("tau"));
        params.Ma                    = std::stod(map.at("Ma"));
        params.polynomialType        = parseVectorInt(map.at("polynomialType"));
        params.parameter1            = parseVectorDouble(map.at("parameter1"));
        params.parameter2            = parseVectorDouble(map.at("parameter2"));

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing XML parameters: " << e.what() << std::endl;
        return false;
    }
}

#endif // PARAMETERS_H