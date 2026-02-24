#!/bin/bash
# mol2db2 C++ Quick Setup Script

set -e  # Exit on error

echo "=========================================="
echo "mol2db2 C++ Setup Script"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if we're in the right directory
if [ ! -f "Makefile" ]; then
    echo -e "${RED}Error: Makefile not found!${NC}"
    echo "Please run this script from the mol2db2-cpp directory"
    exit 1
fi

echo "Step 1: Checking dependencies..."
echo "-----------------------------------"

# Check for C++ compiler
if ! command -v g++ &> /dev/null; then
    echo -e "${RED}g++ not found!${NC}"
    echo "Install with: sudo apt-get install build-essential"
    exit 1
else
    echo -e "${GREEN}✓ g++ found:${NC} $(g++ --version | head -1)"
fi

# Check for Eigen3 — first check local copy (from get_eigen_library.sh),
# then fall back to system locations
if [ -d "./eigen-3.4.0/Eigen" ]; then
    echo -e "${GREEN}✓ Eigen3 found:${NC} ./eigen-3.4.0 (local)"
elif [ -d "/usr/include/eigen3" ]; then
    echo -e "${GREEN}✓ Eigen3 found:${NC} /usr/include/eigen3"
elif [ -d "/usr/local/include/eigen3" ]; then
    echo -e "${GREEN}✓ Eigen3 found:${NC} /usr/local/include/eigen3"
elif [ -d "/opt/homebrew/include/eigen3" ]; then
    echo -e "${GREEN}✓ Eigen3 found:${NC} /opt/homebrew/include/eigen3"
else
    echo -e "${RED}Eigen3 not found!${NC}"
    echo "Run ./get_eigen_library.sh to download it locally, or:"
    echo "  sudo apt-get install libeigen3-dev"
    exit 1
fi

# Check for zlib
if ! ldconfig -p | grep -q libz.so; then
    echo -e "${RED}zlib not found!${NC}"
    echo "Install with: sudo apt-get install zlib1g-dev"
    exit 1
else
    echo -e "${GREEN}✓ zlib found${NC}"
fi

echo ""
echo "Step 2: Creating directory structure..."
echo "-----------------------------------"

mkdir -p build bin test_data
echo -e "${GREEN}✓ Directories created${NC}"

echo ""
echo "Step 3: Compiling mol2db2..."
echo "-----------------------------------"

make clean 2>/dev/null || true
if make; then
    echo -e "${GREEN}✓ Build successful!${NC}"
else
    echo -e "${RED}✗ Build failed!${NC}"
    exit 1
fi

echo ""
echo "Step 4: Running unit tests..."
echo "-----------------------------------"

if make test; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
else
    echo -e "${YELLOW}⚠ Some tests failed${NC}"
fi

echo ""
echo "Step 5: Creating test data..."
echo "-----------------------------------"

# Create minimal test MOL2 file
cat > test_data/test.mol2 << 'EOF'
@<TRIPOS>MOLECULE
BENZENE
12 12 0 0 0
SMALL
USER_CHARGES


mmff94s_NoEstat =   0.00
@<TRIPOS>ATOM
      1 C1      0.0000    1.4000    0.0000 C.ar    1  <0>      -0.1500
      2 C2      1.2124    0.7000    0.0000 C.ar    1  <0>      -0.1500
      3 C3      1.2124   -0.7000    0.0000 C.ar    1  <0>      -0.1500
      4 C4      0.0000   -1.4000    0.0000 C.ar    1  <0>      -0.1500
      5 C5     -1.2124   -0.7000    0.0000 C.ar    1  <0>      -0.1500
      6 C6     -1.2124    0.7000    0.0000 C.ar    1  <0>      -0.1500
      7 H1      0.0000    2.5000    0.0000 H       1  <0>       0.1500
      8 H2      2.1651    1.2500    0.0000 H       1  <0>       0.1500
      9 H3      2.1651   -1.2500    0.0000 H       1  <0>       0.1500
     10 H4      0.0000   -2.5000    0.0000 H       1  <0>       0.1500
     11 H5     -2.1651   -1.2500    0.0000 H       1  <0>       0.1500
     12 H6     -2.1651    1.2500    0.0000 H       1  <0>       0.1500
@<TRIPOS>BOND
     1     1     2   ar
     2     2     3   ar
     3     3     4   ar
     4     4     5   ar
     5     5     6   ar
     6     6     1   ar
     7     1     7    1
     8     2     8    1
     9     3     9    1
    10     4    10    1
    11     5    11    1
    12     6    12    1
EOF

# Create test SOLV file
cat > test_data/test.solv << 'EOF'
BENZENE 12 0.0000 -5.234 285.432 -0.234 -5.468
-0.1500 -0.234 28.543 -0.023 -0.257
-0.1500 -0.234 28.543 -0.023 -0.257
-0.1500 -0.234 28.543 -0.023 -0.257
-0.1500 -0.234 28.543 -0.023 -0.257
-0.1500 -0.234 28.543 -0.023 -0.257
-0.1500 -0.234 28.543 -0.023 -0.257
0.1500 -0.123 18.234 -0.012 -0.135
0.1500 -0.123 18.234 -0.012 -0.135
0.1500 -0.123 18.234 -0.012 -0.135
0.1500 -0.123 18.234 -0.012 -0.135
0.1500 -0.123 18.234 -0.012 -0.135
0.1500 -0.123 18.234 -0.012 -0.135
EOF

# Create test name file
cat > test_data/name.txt << 'EOF'
name.txt dummy BENZENE c1ccccc1 dummy Benzene_Test_Molecule
EOF

echo -e "${GREEN}✓ Test data created${NC}"

echo ""
echo "Step 6: Running mol2db2 on test data..."
echo "-----------------------------------"

if ./bin/mol2db2 \
    -m test_data/test.mol2 \
    -s test_data/test.solv \
    -n test_data/name.txt \
    -o test_data/output.db2.gz \
    -v; then
    echo -e "${GREEN}✓ mol2db2 ran successfully!${NC}"
else
    echo -e "${RED}✗ mol2db2 failed!${NC}"
    exit 1
fi

echo ""
echo "Step 7: Verifying output..."
echo "-----------------------------------"

if [ -f "test_data/output.db2.gz" ]; then
    SIZE=$(du -h test_data/output.db2.gz | cut -f1)
    echo -e "${GREEN}✓ Output file created:${NC} test_data/output.db2.gz ($SIZE)"
    
    # Check if it's actually gzipped
    if file test_data/output.db2.gz | grep -q "gzip"; then
        echo -e "${GREEN}✓ Output is properly gzipped${NC}"
    else
        echo -e "${YELLOW}⚠ Output may not be gzipped correctly${NC}"
    fi
    
    # Show first few lines
    echo ""
    echo "First 10 lines of output:"
    zcat test_data/output.db2.gz | head -10
else
    echo -e "${RED}✗ Output file not created!${NC}"
    exit 1
fi

echo ""
echo "=========================================="
echo -e "${GREEN}✓✓✓ Setup Complete! ✓✓✓${NC}"
echo "=========================================="
echo ""
echo "mol2db2 is installed and working!"
echo ""
echo "Next steps:"
echo "  1. Process your own data:"
echo "     ./bin/mol2db2 -m your_file.mol2.gz -s your_file.solv -o output.db2.gz -v"
echo ""
echo "  2. See all options:"
echo "     ./bin/mol2db2 --help"
echo ""
echo "  3. Install to ~/.local/bin (optional, no sudo needed):"
echo "     make install"
echo ""
