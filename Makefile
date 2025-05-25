# Makefile for ROOT to CSV converter

# ROOT configuration
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS := $(shell root-config --libs)

# Compiler settings
CXX = g++
CXXFLAGS = -Wall -O2 $(ROOTCFLAGS)
LDFLAGS = $(ROOTLDFLAGS)
LIBS = $(ROOTLIBS)

# Target
TARGET = rootToCSV
SOURCE = rootToCSV.C

# Build rule
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS) $(LIBS)

# Clean rule
clean:
	rm -f $(TARGET) jet_data.csv

# Run rule
run: $(TARGET)
	./$(TARGET)

.PHONY: clean run