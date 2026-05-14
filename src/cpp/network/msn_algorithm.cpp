#include "msn_algorithm.h"
#include <algorithm>
#include <queue>

namespace haplokit {
namespace network {

// UnionFind implementation
MSNAlgorithm::UnionFind::UnionFind(int n) : parent_(n), rank_(n, 0), components_(n) {
    for (int i = 0; i < n; ++i) {
        parent_[i] = i;
    }
}

int MSNAlgorithm::UnionFind::find(int x) {
    if (parent_[x] != x) {
        parent_[x] = find(parent_[x]); // path compression
    }
    return parent_[x];
}

bool MSNAlgorithm::UnionFind::unite(int x, int y) {
    int root_x = find(x);
    int root_y = find(y);

    if (root_x == root_y) return false;

    // Union by rank
    if (rank_[root_x] < rank_[root_y]) {
        parent_[root_x] = root_y;
    } else if (rank_[root_x] > rank_[root_y]) {
        parent_[root_y] = root_x;
    } else {
        parent_[root_y] = root_x;
        rank_[root_x]++;
    }

    components_--;
    return true;
}

// MSN computation using Kruskal's algorithm
Graph MSNAlgorithm::compute(const std::vector<Haplotype>& haplotypes, int epsilon) {
    Graph graph;
    int n = haplotypes.size();

    // Add vertices
    for (const auto& hap : haplotypes) {
        graph.add_vertex(hap.sequence, hap.samples);
    }

    // Compute all pairwise distances and create edge candidates
    std::vector<EdgeCandidate> candidates;
    candidates.reserve(n * (n - 1) / 2);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int dist = HammingSIMD::hamming_distance(
                haplotypes[i].sequence,
                haplotypes[j].sequence
            );
            candidates.push_back({i, j, dist});
        }
    }

    // Sort edges by weight
    std::sort(candidates.begin(), candidates.end());

    // Kruskal's algorithm with epsilon relaxation
    UnionFind uf(n);
    int max_threshold = std::numeric_limits<int>::max();

    for (const auto& edge : candidates) {
        // Stop if threshold exceeded
        if (edge.weight > max_threshold) break;

        // Add edge if it connects different components
        if (uf.unite(edge.u, edge.v)) {
            graph.add_edge(edge.u, edge.v, edge.weight);

            // Once connected, allow epsilon relaxation
            if (uf.component_count() == 1 && max_threshold == std::numeric_limits<int>::max()) {
                max_threshold = edge.weight + epsilon;
            }
        }
    }

    return graph;
}

} // namespace network
} // namespace haplokit
