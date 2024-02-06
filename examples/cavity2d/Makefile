# Makefile for compiling cavity.cpp with OpenMP and Python3.8 support

# Compiler settings
CXX = g++
CXXFLAGS = -O3 -fopenmp
INCLUDES = -I/usr/include/python3.8/
LIBS = -lpython3.8

# Target executable
TARGET = cavity

# Source files
SOURCES = cavity.cpp
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
