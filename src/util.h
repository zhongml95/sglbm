#include <unistd.h>
#include <limits.h>
#include <libgen.h> // For dirname
#include <sys/stat.h> // For stat
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>


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
    int polynomialType;
    int resolution;
    double parameter1;
    double parameter2;
    double L;
    double lx;
    double ly;
    double Re;
    double physVelocity;
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

    // Assuming all keys exist and are correctly formatted in the file
    params.points_weights_method = std::stoul(paramMap["points_weights_method"]);
    params.order = std::stoul(paramMap["order"]);
    params.nq = std::stoul(paramMap["nq"]);
    params.polynomialType = std::stoi(paramMap["polynomialType"]);
    params.resolution = std::stoi(paramMap["resolution"]);
    params.parameter1 = std::stod(paramMap["parameter1"]);
    params.parameter2 = std::stod(paramMap["parameter2"]);
    params.L = std::stod(paramMap["L"]);
    params.lx = std::stod(paramMap["lx"]);
    params.ly = std::stod(paramMap["ly"]);
    params.Re = std::stod(paramMap["Re"]);
    params.physVelocity = std::stod(paramMap["physVelocity"]);

    return true;
}


