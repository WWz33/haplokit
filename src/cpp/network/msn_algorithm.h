#ifndef HAPLOKIT_NETWORK_MSN_H
#define HAPLOKIT_NETWORK_MSN_H

#include "graph.h"
#include "haplotype.h"
#include "hamming_simd.h"
#include <vector>
#include <string>

namespace haplokit {
namespace network {

// Minimum Spanning Network (MSN) algorithm
// Constructs MST using Kruskal's algorithm with union-find
class MSNAlgorithm {
public:
    // Compute MSN from haplotypes
    static Graph compute(const std::vector<Haplotype>& haplotypes, int epsilon = 0);

private:
    // Union-Find data structure for Kruskal's algorithm
    class UnionFind {
    public:
        explicit UnionFind(int n);
        int find(int x);
        bool unite(int x, int y);
        int component_count() const { return components_; }

    private:
        std::vector<int> parent_;
        std::vector<int> rank_;
        int components_;
    };

    struct EdgeCandidate {
        int u, v, weight;
        bool operator<(const EdgeCandidate& other) const {
            return weight < other.weight;
        }
    };
};

} // namespace network
} // namespace haplokit

#endif // HAPLOKIT_NETWORK_MSN_H
