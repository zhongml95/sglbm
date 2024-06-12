#include <unistd.h>
#include <limits.h>
#include <libgen.h> // For dirname
#include <sys/stat.h> // For stat
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>


// Function to calculate factorial of a number using double
double factorial(int n) {
    if (n == 0) {
        return 1.0;
    }
    double result = 1.0;
    for (int i = 1; i <= n; ++i) {
        result *= i;
    }
    return result;
}

bool directoryExists(const std::string& path) {
    struct stat statbuf;
    if (stat(path.c_str(), &statbuf) != 0) {
        return false;
    }
    return S_ISDIR(statbuf.st_mode);
}

/*bool directoryExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
        return false;
    else if (info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}*/

bool createDirectory(const std::string& path) {
    int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status == 0)
        return true;
    else
        return false;
}

bool deleteDirectory(const std::string& path) {
    std::string command = "rm -rf " + path;
    int status = std::system(command.c_str());
    if (status == 0)
        return true;
    else
        return false;
}

std::string findRelativePathToSrc() {
    char exePath[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", exePath, PATH_MAX);
    if (count == -1) {
        std::cerr << "Failed to determine the path of the executable" << std::endl;
        return "";
    }
    exePath[count] = '\0'; // Null-terminate the read path

    std::string currentPath = dirname(exePath); // Starting directory
    std::string relativePath = ""; // Initialize empty relative path

    // Try moving up the directory tree
    for (int i = 0; i < 100; ++i) { // Limit the depth to prevent potential infinite loop
        std::string testPath = currentPath + "/src"; // Construct test path to src
        if (directoryExists(testPath)) {
            // If the src directory is found, return the relative path
            return relativePath + "src";
        }

        // Update the paths for the next iteration
        relativePath += "../";
        currentPath = dirname((char*)currentPath.c_str()); // Move up one directory
    }

    std::cerr << "Failed to locate the 'src' directory" << std::endl;
    return "";
}

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

    // Parse scalar values directly
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

    // Parse vector values
    auto parseVectorDouble = [](const std::string& s) -> std::vector<double> {
        std::vector<double> result;
        std::istringstream iss(s);
        std::string item;
        while (std::getline(iss, item, ' ')) { // Space-separated
            result.push_back(std::stod(item));
        }
        return result;
    };

    auto parseVectorInt = [](const std::string& s) -> std::vector<int> {
        std::vector<int> result;
        std::istringstream iss(s);
        std::string item;
        while (std::getline(iss, item, ' ')) { // Space-separated
            if (!item.empty()) { // Check if the item is not just whitespace
                result.push_back(std::stoi(item));
            }
        }
        return result;
    };

    // Adjusted for parsing space-separated values correctly
    params.polynomialType = parseVectorInt(paramMap["polynomialType"]);
    params.parameter1 = parseVectorDouble(paramMap["parameter1"]);
    params.parameter2 = parseVectorDouble(paramMap["parameter2"]);

    return true;
}

// Utility function to check if a file exists.
bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

void saveVector1D(const std::string& filePath, const std::vector<double>& vec) {
    std::ofstream out(filePath, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("Cannot open file for writing: " + filePath);
    
    size_t size = vec.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(double));
}

void readVector1D(const std::string& filePath, std::vector<double>& vec) {
    std::ifstream in(filePath, std::ios::binary);
    if (!in.is_open()) throw std::runtime_error("Cannot open file for reading: " + filePath);
    
    size_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);
    in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(double));
}

void saveVector3D(const std::string& filePath, const std::vector<std::vector<std::vector<double>>>& vec) {
    std::ofstream out(filePath, std::ios::binary);
    if (!out.is_open()) throw std::runtime_error("Cannot open file for writing: " + filePath);

    size_t outerSize = vec.size();
    out.write(reinterpret_cast<const char*>(&outerSize), sizeof(outerSize));
    
    for (const auto& midVec : vec) {
        size_t midSize = midVec.size();
        out.write(reinterpret_cast<const char*>(&midSize), sizeof(midSize));
        
        for (const auto& innerVec : midVec) {
            size_t innerSize = innerVec.size();
            out.write(reinterpret_cast<const char*>(&innerSize), sizeof(innerSize));
            out.write(reinterpret_cast<const char*>(innerVec.data()), innerSize * sizeof(double));
        }
    }
}

void readVector3D(const std::string& filePath, std::vector<std::vector<std::vector<double>>>& vec) {
    std::ifstream in(filePath, std::ios::binary);
    if (!in.is_open()) throw std::runtime_error("Cannot open file for reading: " + filePath);
    
    size_t outerSize;
    in.read(reinterpret_cast<char*>(&outerSize), sizeof(outerSize));
    vec.resize(outerSize);
    
    for (auto& midVec : vec) {
        size_t midSize;
        in.read(reinterpret_cast<char*>(&midSize), sizeof(midSize));
        midVec.resize(midSize);
        
        for (auto& innerVec : midVec) {
            size_t innerSize;
            in.read(reinterpret_cast<char*>(&innerSize), sizeof(innerSize));
            innerVec.resize(innerSize);
            in.read(reinterpret_cast<char*>(innerVec.data()), innerSize * sizeof(double));
        }
    }
}

