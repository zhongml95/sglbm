# Compiler settings
CXX = g++
CXXFLAGS = -O3 -fopenmp
INCLUDES = -I/home/zhongml95/gsl/include
LIBS = -lgsl -lgslcblas

# Target executable
TARGET = test

# Source directory
SRC_DIR = ../../src

# Find all .cpp files in the source directory
SOURCES = test.cpp $(wildcard $(SRC_DIR)/*.cpp)

# Create a list of object files
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Linking the target executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Compiling source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Cleaning up
clean:
	rm -f $(TARGET) $(OBJECTS)

# Phony targets
.PHONY: all clean
