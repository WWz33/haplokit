#!/bin/bash
# Build script for haplokit network algorithms

set -e

echo "=== Building haplokit with network algorithms ==="

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    echo "Error: CMakeLists.txt not found. Run from project root."
    exit 1
fi

# Create build directory
BUILD_DIR="build-wsl"
mkdir -p "$BUILD_DIR"

# Configure CMake
echo "Configuring CMake..."
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17

# Build
echo "Building..."
cmake --build "$BUILD_DIR" --parallel $(nproc)

# Test
echo "Running tests..."
cd "$BUILD_DIR"
ctest --output-on-failure

echo ""
echo "=== Build complete ==="
echo "Executables:"
echo "  - haplokit_cpp: Main VCF viewer"
echo "  - haplokit_network_backend: Network computation backend"
echo ""
echo "To use network algorithms:"
echo "  export HAPLOKIT_NETWORK_BIN=$PWD/haplokit_network_backend"
echo "  python -m haplokit.network"
