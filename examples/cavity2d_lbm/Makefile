# Makefile for compiling cavity2d.cpp with OpenMP and GSL

# Compiler settings
CXX = g++
CXXFLAGS = -O3 -fopenmp
INCLUDES = -I/home/zhongml95/gsl/include
LIBS = -lgsl -lgslcblas # Link against GSL libraries

# Target executable
TARGET = cavity2d

# Source files
SOURCES = cavity2d.cpp
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
