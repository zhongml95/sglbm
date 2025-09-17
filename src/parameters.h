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
    std::string quadratureRule;
    std::vector<std::string> distributionType;
    std::vector<double> parameter1;
    std::vector<double> parameter2;
    unsigned int order;
    unsigned int nq;
    bool sparse = false;
    int resolution;
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

// ---------- helpers (header-scope) ----------
inline std::string trim_copy(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\n\r");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\n\r");
    return s.substr(b, e - b + 1);
}

// Accepts whitespace/comma/semicolon; ignores bracket chars [](){}; supports quotes for strings
inline std::vector<std::string> parseVectorString(const std::string& s) {
    std::vector<std::string> out;
    std::string token;
    bool inQuote = false;
    char quoteChar = '\0';

    auto flush = [&]() {
        std::string t = trim_copy(token);
        if (!t.empty()) out.push_back(t);
        token.clear();
    };

    for (char c : s) {
        if (inQuote) {
            if (c == quoteChar) { inQuote = false; }
            else { token.push_back(c); }
            continue;
        }
        if (c == '"' || c == '\'') { inQuote = true; quoteChar = c; continue; }

        // separators
        if (std::isspace(static_cast<unsigned char>(c)) || c == ',' || c == ';') { flush(); continue; }

        // ignore bracket-like grouping characters
        if (c == '[' || c == ']' || c == '(' || c == ')' || c == '{' || c == '}') { flush(); continue; }

        token.push_back(c);
    }
    flush();
    return out;
}

inline std::vector<double> parseVectorDouble(const std::string& s) {
    std::vector<double> vec;
    std::string buf; buf.reserve(s.size());
    // Normalize delimiters, strip brackets
    for (char c : s) {
        if (c == ',' || c == ';') buf.push_back(' ');
        else if (c == '[' || c == ']' || c == '(' || c == ')' || c == '{' || c == '}') continue;
        else buf.push_back(c);
    }
    std::istringstream iss(buf);
    std::string item;
    while (iss >> item) vec.push_back(std::stod(item));
    return vec;
}

inline std::vector<int> parseVectorInt(const std::string& s) {
    std::vector<int> vec;
    std::string buf; buf.reserve(s.size());
    for (char c : s) {
        if (c == ',' || c == ';') buf.push_back(' ');
        else if (c == '[' || c == ']' || c == '(' || c == ')' || c == '{' || c == '}') continue;
        else buf.push_back(c);
    }
    std::istringstream iss(buf);
    std::string item;
    while (iss >> item) vec.push_back(std::stoi(item));
    return vec;
}

// ---------- rewritten reader ----------
inline bool readParametersFromXML(const std::string& filePath, Parameters& params) {
    try {
        auto map = parseSimpleXML(filePath);

        params.quadratureRule   = map.at("quadratureRule");
        params.distributionType = parseVectorString(map.at("distributionType"));
        params.parameter1       = parseVectorDouble(map.at("parameter1"));
        params.parameter2       = parseVectorDouble(map.at("parameter2"));
        params.order            = static_cast<unsigned int>(std::stoul(map.at("order")));
        params.nq               = static_cast<unsigned int>(std::stoul(map.at("nq")));
        params.sparse           = (map.find("sparse") != map.end()) ? (map.at("sparse") == "true" || map.at("sparse") == "1") : false;
        params.resolution       = std::stoi(map.at("resolution"));
        params.L                = std::stod(map.at("L"));
        params.lx               = std::stod(map.at("lx"));
        params.ly               = std::stod(map.at("ly"));
        params.Re               = std::stod(map.at("Re"));
        params.physVelocity     = std::stod(map.at("physVelocity"));
        params.physViscosity    = std::stod(map.at("physViscosity"));
        params.tau              = std::stod(map.at("tau"));
        params.Ma               = std::stod(map.at("Ma"));

        return true;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: missing XML tag: " << e.what() << std::endl;
        return false;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing XML parameters: " << e.what() << std::endl;
        return false;
    }
};


#endif // PARAMETERS_H