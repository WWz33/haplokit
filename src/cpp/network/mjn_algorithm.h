#ifndef HAPLOKIT_NETWORK_MJN_H
#define HAPLOKIT_NETWORK_MJN_H

#include "graph.h"
#include "haplotype.h"
#include "msn_algorithm.h"
#include "hamming_simd.h"
#include <vector>
#include <string>
#include <set>

namespace haplokit {
namespace network {

// Median-Joining Network (MJN) algorithm
// Iteratively adds median vectors to reduce network length
class MJNAlgorithm {
public:
    // Compute MJN from haplotypes
    static Graph compute(const std::vector<Haplotype>& haplotypes, int epsilon = 0);

private:
    // Compute quasi-median sequences from three sequences
    static std::set<std::string> compute_quasi_medians(
        const std::string& seq1,
        const std::string& seq2,
        const std::string& seq3
    );

    // Compute cost of adding a median vertex
    static int compute_median_cost(
        const std::string& seq1,
        const std::string& seq2,
        const std::string& seq3,
        const std::string& median
    );

    // Remove obsolete vertices (degree < 3 median vertices)
    static bool remove_obsolete_vertices(Graph& graph);

    // Check if adding median reduces network length
    static bool is_feasible_median(
        const Graph& graph,
        const std::string& median,
        int u, int v, int w
    );
};

} // namespace network
} // namespace haplokit

#endif // HAPLOKIT_NETWORK_MJN_H
