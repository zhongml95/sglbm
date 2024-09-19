# Compiler settings
CXX = g++
CXXFLAGS = -O3 -fopenmp -DUSE_GSL  # Add -DUSE_GSL here for global definition
INCLUDES = -I/home/zhongml95/gsl/include
LDFLAGS = -L/home/zhongml95/gsl/lib
LIBS = -lgsl -lgslcblas  # Link against GSL libraries

# Target executable
TARGET = tgv2d

# Source directory
SRC_DIR = ../../src

# Build directory for object files
BUILD_DIR = build

# Source files
SOURCES = tgv2d.cpp $(wildcard $(SRC_DIR)/*.cpp)

# Object files will be placed in the build directory
OBJECTS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(notdir $(SOURCES)))

# Include directories
VPATH = $(SRC_DIR)

# Default target
all: $(BUILD_DIR) $(TARGET)

# Create the build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Linking the target executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LIBS)

# Compiling source files, placing the output in the build directory
$(BUILD_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Cleaning up
clean:
	rm -rf $(TARGET) $(BUILD_DIR)

# Phony targets
.PHONY: all clean
