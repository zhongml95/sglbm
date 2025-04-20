# ────────────────────────────────────────────────────────────────────────────
#  Build settings
# ────────────────────────────────────────────────────────────────────────────
CXX      ?= g++
STD        = -std=c++20

# Speed‑oriented flags (tune as needed)
OPTFLAGS  = -O3 -march=native -ffast-math -flto
OMPFLAG   = -fopenmp                           # remove if you don't use OpenMP
WARNFLAGS = -Wall -Wextra -Wpedantic

CXXFLAGS += $(STD) $(OPTFLAGS) $(OMPFLAG) $(WARNFLAGS) -MMD -MP
LDFLAGS  += $(OMPFLAG)    # link with OpenMP runtime

# Headers that live in ../../src
INCLUDES  = -I../../src
# Extra libs (e.g. -lm, -lfftw3) can go here
LIBS      =

# ────────────────────────────────────────────────────────────────────────────
#  File lists
# ────────────────────────────────────────────────────────────────────────────
SRC_DIR   := ../../src
BUILD_DIR := build
TARGET    := tgv2d

# Driver plus every *.cpp inside SRC_DIR
SOURCES   := tgv2d.cpp $(wildcard $(SRC_DIR)/*.cpp)

# Map each source file to an object file under BUILD_DIR mirroring the path
OBJECTS   := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPS      := $(OBJECTS:.o=.d)

# ────────────────────────────────────────────────────────────────────────────
#  Default target
# ────────────────────────────────────────────────────────────────────────────
all: $(TARGET)

# ────────────────────────────────────────────────────────────────────────────
#  Link the final executable
# ────────────────────────────────────────────────────────────────────────────
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

# ────────────────────────────────────────────────────────────────────────────
#  Compile every .cpp –> build/…/.o (auto‑create sub‑dirs)
# ────────────────────────────────────────────────────────────────────────────
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# ────────────────────────────────────────────────────────────────────────────
#  House‑keeping
# ────────────────────────────────────────────────────────────────────────────
clean:
	$(RM) -r $(BUILD_DIR) $(TARGET)

# If any dependency file is missing just ignore the error
-include $(DEPS)

.PHONY: all clean
