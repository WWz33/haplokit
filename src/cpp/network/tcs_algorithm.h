#ifndef HAPLOKIT_NETWORK_TCS_H
#define HAPLOKIT_NETWORK_TCS_H

#include "graph.h"
#include "haplotype.h"
#include "hamming_simd.h"
#include <vector>
#include <string>

namespace haplokit {
namespace network {

// TCS (Templeton-Crandall-Sing) Statistical Parsimony Network
// Connects haplotypes at increasing distance thresholds
// Inserts intermediate vertices when distance > 1
class TCSAlgorithm {
public:
    // Compute TCS network from haplotypes
    static Graph compute(const std::vector<Haplotype>& haplotypes);

private:
    // Find or create intermediate vertices between u and v
    static void insert_intermediates(Graph& graph, int u, int v, int distance);

    // Generate intermediate sequence at step i of distance steps
    static std::string generate_intermediate(const std::string& seq1,
                                             const std::string& seq2,
                                             int step, int total_steps);
};

} // namespace network
} // namespace haplokit

#endif // HAPLOKIT_NETWORK_TCS_H
