CC = g++
CFLAGS = -std=c++11 -Wall -Wextra
SRCDIR = ./src
INCLUDES = -I$(SRCDIR)

SOURCES = main.cpp
OBJECTS = $(SOURCES:%.cpp=%.o)
EXECUTABLE = main

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@

%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
