#!/bin/bash
# Download Eigen 3.4.0 into the current directory.
# After running this, mol2db2 can be built with: make

set -e

EIGEN_VERSION="3.4.0"
EIGEN_DIR="eigen-${EIGEN_VERSION}"
EIGEN_ZIP="${EIGEN_DIR}.zip"
EIGEN_URL="https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/${EIGEN_ZIP}"

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

if [ -d "$EIGEN_DIR" ]; then
    echo -e "${GREEN}✓ ${EIGEN_DIR}/ already exists — nothing to do.${NC}"
    exit 0
fi

echo "Downloading Eigen ${EIGEN_VERSION}..."

if command -v wget &> /dev/null; then
    wget -q --show-progress "$EIGEN_URL" -O "$EIGEN_ZIP"
elif command -v curl &> /dev/null; then
    curl -L --progress-bar "$EIGEN_URL" -o "$EIGEN_ZIP"
else
    echo -e "${RED}Error: neither wget nor curl is available.${NC}"
    exit 1
fi

echo "Extracting..."
unzip -q "$EIGEN_ZIP"
rm "$EIGEN_ZIP"

echo -e "${GREEN}✓ Eigen ${EIGEN_VERSION} ready at ./${EIGEN_DIR}/${NC}"
echo "You can now run: make"
