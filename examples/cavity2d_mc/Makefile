# Makefile for compiling tgv2d_lbm.cpp with OpenMP support

# Compiler settings
CXX = g++
CXXFLAGS = -O3 -fopenmp

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
