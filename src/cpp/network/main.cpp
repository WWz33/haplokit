#include "network_backend.h"

// Main entry point for network backend executable
int main(int argc, char* argv[]) {
    return haplokit::network::NetworkBackend::run(argc, argv);
}
