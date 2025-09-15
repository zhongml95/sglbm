#ifndef UTILS_H
#define UTILS_H

#include "utils.h"
#include <cstdlib>  // For std::system



std::string findRelativePathToSrc() {
    char exePath[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", exePath, PATH_MAX);
    if (count == -1) {
        std::cerr << "Failed to determine the path of the executable" << std::endl;
        return "";
    }
    exePath[count] = '\0';

    std::string currentPath = dirname(exePath);
    std::string relativePath = "";

    for (int i = 0; i < 100; ++i) {
        std::string testPath = currentPath + "/src";
        if (directoryExists(testPath)) {
            return relativePath + "src";
        }

        relativePath += "../";
        currentPath = dirname((char*)currentPath.c_str());
    }

    std::cerr << "Failed to locate the 'src' directory" << std::endl;
    return "";
}



#endif  // utils_hh