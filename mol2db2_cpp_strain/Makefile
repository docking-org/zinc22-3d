# Makefile for mol2db2 C++ conversion

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -I./src/core -I./src/algorithms -I./src/molecular -I./src/hierarchy
#INCLUDES = -I/usr/include/eigen3
INCLUDES = -I./eigen-3.4.0
LIBS = -lz

# Output directory
BUILD_DIR = build
BIN_DIR = bin

# Source directories
CORE_DIR = src/core
ALG_DIR = src/algorithms
MOL_DIR = src/molecular
HIER_DIR = src/hierarchy
MAIN_DIR = src/main
TEST_DIR = tests

# Object files
CORE_OBJS = $(BUILD_DIR)/geometry.o $(BUILD_DIR)/unionfind.o $(BUILD_DIR)/floydwarshall.o
ALG_OBJS = $(BUILD_DIR)/shortestpaths.o $(BUILD_DIR)/buckets.o $(BUILD_DIR)/buckets2.o $(BUILD_DIR)/pca.o $(BUILD_DIR)/divisive_clustering.o
MOL_OBJS = $(BUILD_DIR)/clash.o $(BUILD_DIR)/mol2.o $(BUILD_DIR)/hydrogens.o
HIER_OBJS = $(BUILD_DIR)/hierarchy.o

ALL_OBJS = $(CORE_OBJS) $(ALG_OBJS) $(MOL_OBJS) $(HIER_OBJS)

# Install prefix — override with: make install PREFIX=/your/path
PREFIX = $(HOME)/.local

# Targets
.PHONY: all clean test install dirs

all: dirs $(BIN_DIR)/mol2db2

dirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

# Main executable
$(BIN_DIR)/mol2db2: $(ALL_OBJS) $(BUILD_DIR)/mol2db2.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	@echo "Build complete: $(BIN_DIR)/mol2db2"

# Core objects
$(BUILD_DIR)/geometry.o: $(CORE_DIR)/geometry.cpp $(CORE_DIR)/geometry.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/unionfind.o: $(CORE_DIR)/unionfind.cpp $(CORE_DIR)/unionfind.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

#$(BUILD_DIR)/floydwarshall.o: $(CORE_DIR)/floydwarshall.cpp $(CORE_DIR)/floydwarshall.h
$(BUILD_DIR)/floydwarshall.o: $(ALG_DIR)/floydwarshall.cpp $(ALG_DIR)/floydwarshall.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Algorithm objects
$(BUILD_DIR)/shortestpaths.o: $(ALG_DIR)/shortestpaths.cpp $(ALG_DIR)/shortestpaths.h $(CORE_DIR)/priodict.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/buckets.o: $(ALG_DIR)/buckets.cpp $(ALG_DIR)/buckets.h $(CORE_DIR)/geometry.h $(CORE_DIR)/unionfind.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/buckets2.o: $(ALG_DIR)/buckets2.cpp $(ALG_DIR)/buckets2.h $(CORE_DIR)/geometry.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/pca.o: $(ALG_DIR)/pca.cpp $(ALG_DIR)/pca.h $(CORE_DIR)/geometry.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/divisive_clustering.o: $(ALG_DIR)/divisive_clustering.cpp $(ALG_DIR)/divisive_clustering.h $(ALG_DIR)/pca.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Molecular objects
$(BUILD_DIR)/clash.o: $(MOL_DIR)/clash.cpp $(MOL_DIR)/clash.h $(MOL_DIR)/mol2.h $(CORE_DIR)/geometry.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/mol2.o: $(MOL_DIR)/mol2.cpp $(MOL_DIR)/mol2.h $(CORE_DIR)/geometry.h $(ALG_DIR)/floydwarshall.h $(ALG_DIR)/shortestpaths.h $(ALG_DIR)/divisive_clustering.h $(CORE_DIR)/unionfind.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/hydrogens.o: $(MOL_DIR)/hydrogens.cpp $(MOL_DIR)/hydrogens.h $(MOL_DIR)/mol2.h $(CORE_DIR)/geometry.h $(CORE_DIR)/combinatorics.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Hierarchy objects
$(BUILD_DIR)/hierarchy.o: $(HIER_DIR)/hierarchy.cpp $(HIER_DIR)/hierarchy.h $(MOL_DIR)/mol2.h $(MOL_DIR)/clash.h $(CORE_DIR)/geometry.h $(ALG_DIR)/buckets.h $(ALG_DIR)/buckets2.h $(CORE_DIR)/unionfind.h $(ALG_DIR)/divisive_clustering.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Main program
$(BUILD_DIR)/mol2db2.o: $(MAIN_DIR)/mol2db2.cpp $(MOL_DIR)/mol2.h $(MOL_DIR)/clash.h $(MOL_DIR)/hydrogens.h $(HIER_DIR)/hierarchy.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Tests
test: dirs $(BIN_DIR)/test_geometry $(BIN_DIR)/test_unionfind $(BIN_DIR)/test_mol2
	@echo "Running tests..."
	@$(BIN_DIR)/test_geometry
	@$(BIN_DIR)/test_unionfind
	@$(BIN_DIR)/test_mol2

$(BIN_DIR)/test_geometry: $(BUILD_DIR)/test_geometry.o $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/test_geometry.o: $(TEST_DIR)/test_geometry.cpp $(CORE_DIR)/geometry.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BIN_DIR)/test_unionfind: $(BUILD_DIR)/test_unionfind.o $(BUILD_DIR)/unionfind.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/test_unionfind.o: $(TEST_DIR)/test_unionfind.cpp $(CORE_DIR)/unionfind.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BIN_DIR)/test_mol2: $(BUILD_DIR)/test_mol2.o $(ALL_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/test_mol2.o: $(TEST_DIR)/test_mol2.cpp $(MOL_DIR)/mol2.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	@echo "Clean complete"

# Install — copies binary to PREFIX/bin (default: ~/.local/bin, no sudo needed)
install: all
	@mkdir -p $(PREFIX)/bin
	cp $(BIN_DIR)/mol2db2 $(PREFIX)/bin/mol2db2
	@echo "Installed to $(PREFIX)/bin/mol2db2"
	@echo "Ensure $(PREFIX)/bin is in your PATH."

# Help
help:
	@echo "Available targets:"
	@echo "  all      - Build mol2db2 executable (default)"
	@echo "  test     - Build and run unit tests"
	@echo "  clean    - Remove all build files"
	@echo "  install  - Copy binary to PREFIX/bin (default: ~/.local/bin)"
	@echo "             Override: make install PREFIX=/your/path"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Dependencies:"
	@echo "  - Eigen3 library (place eigen-3.4.0/ inside mol2db2_cpp_strain/)"
	@echo "  - zlib library (install: sudo apt-get install zlib1g-dev)"
