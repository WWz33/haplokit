#include "mjn_algorithm.h"
#include <algorithm>
#include <limits>
#include <map>

namespace haplokit {
namespace network {

std::set<std::string> MJNAlgorithm::compute_quasi_medians(
    const std::string& seq1,
    const std::string& seq2,
    const std::string& seq3) {

    std::set<std::string> medians;
    size_t len = std::min({seq1.size(), seq2.size(), seq3.size()});

    // Generate all possible median sequences
    // At each position, take majority vote or all possibilities if no majority
    std::vector<std::vector<char>> position_choices(len);

    for (size_t i = 0; i < len; ++i) {
        char c1 = seq1[i], c2 = seq2[i], c3 = seq3[i];

        if (c1 == c2 || c1 == c3) {
            position_choices[i] = {c1};
        } else if (c2 == c3) {
            position_choices[i] = {c2};
        } else {
            // No majority: all three are candidates
            position_choices[i] = {c1, c2, c3};
        }
    }

    // Generate all combinations (limit to prevent explosion)
    std::function<void(size_t, std::string&)> generate;
    generate = [&](size_t pos, std::string& current) {
        if (pos == len) {
            medians.insert(current);
            return;
        }

        for (char c : position_choices[pos]) {
            current.push_back(c);
            generate(pos + 1, current);
            current.pop_back();

            // Limit median count to prevent combinatorial explosion
            if (medians.size() > 100) return;
        }
    };

    std::string temp;
    generate(0, temp);

    return medians;
}

int MJNAlgorithm::compute_median_cost(
    const std::string& seq1,
    const std::string& seq2,
    const std::string& seq3,
    const std::string& median) {

    int cost = 0;
    cost += HammingSIMD::hamming_distance(seq1, median);
    cost += HammingSIMD::hamming_distance(seq2, median);
    cost += HammingSIMD::hamming_distance(seq3, median);

    return cost;
}

bool MJNAlgorithm::remove_obsolete_vertices(Graph& graph) {
    bool changed = false;

    // Remove median vertices with degree < 3
    std::vector<int> to_check;
    for (const auto& v : graph.vertices()) {
        if (v.is_median) {
            to_check.push_back(v.id);
        }
    }

    for (int vid : to_check) {
        int deg = graph.degree(vid);

        if (deg < 3) {
            auto neighbors = graph.get_neighbors(vid);

            // Remove all edges to this vertex
            for (int neighbor : neighbors) {
                graph.remove_edge(vid, neighbor);
            }

            changed = true;
        }
    }

    return changed;
}

Graph MJNAlgorithm::compute(const std::vector<Haplotype>& haplotypes, int epsilon) {
    Graph graph;

    // Start with MSN
    graph = MSNAlgorithm::compute(haplotypes, epsilon);

    // Track all sequences in network
    std::set<std::string> all_sequences;
    for (const auto& hap : haplotypes) {
        all_sequences.insert(hap.sequence);
    }

    bool changed = true;
    int iteration = 0;
    const int MAX_ITERATIONS = 100;

    while (changed && iteration < MAX_ITERATIONS) {
        changed = false;
        iteration++;

        // Find minimum cost median to add
        int min_cost = std::numeric_limits<int>::max();
        std::string best_median;
        int best_u = -1, best_v = -1, best_w = -1;

        // Examine all triplets of connected vertices
        for (const auto& vert_u : graph.vertices()) {
            int u = vert_u.id;
            auto neighbors_u = graph.get_neighbors(u);

            for (size_t i = 0; i < neighbors_u.size(); ++i) {
                int v = neighbors_u[i];

                for (size_t j = i + 1; j < neighbors_u.size(); ++j) {
                    int w = neighbors_u[j];

                    // Compute quasi-medians for this triplet
                    const Vertex* vu = graph.get_vertex(u);
                    const Vertex* vv = graph.get_vertex(v);
                    const Vertex* vw = graph.get_vertex(w);

                    if (!vu || !vv || !vw) continue;

                    auto medians = compute_quasi_medians(
                        vu->sequence, vv->sequence, vw->sequence
                    );

                    for (const auto& median : medians) {
                        // Skip if median already exists
                        if (all_sequences.count(median)) continue;

                        int cost = compute_median_cost(
                            vu->sequence, vv->sequence, vw->sequence, median
                        );

                        if (cost < min_cost) {
                            min_cost = cost;
                            best_median = median;
                            best_u = u;
                            best_v = v;
                            best_w = w;
                        }
                    }
                }
            }
        }

        // Add best median if cost improvement is significant
        if (!best_median.empty() && min_cost <= epsilon) {
            int median_id = graph.add_median_vertex(best_median);
            all_sequences.insert(best_median);

            // Connect median to the three vertices
            const Vertex* vu = graph.get_vertex(best_u);
            const Vertex* vv = graph.get_vertex(best_v);
            const Vertex* vw = graph.get_vertex(best_w);

            if (vu && vv && vw) {
                int dist_u = HammingSIMD::hamming_distance(best_median, vu->sequence);
                int dist_v = HammingSIMD::hamming_distance(best_median, vv->sequence);
                int dist_w = HammingSIMD::hamming_distance(best_median, vw->sequence);

                if (dist_u > 0) graph.add_edge(median_id, best_u, dist_u);
                if (dist_v > 0) graph.add_edge(median_id, best_v, dist_v);
                if (dist_w > 0) graph.add_edge(median_id, best_w, dist_w);
            }

            changed = true;
        }

        // Rebuild MSN with new vertices
        if (changed) {
            // Recompute MSN including median vertices
            // (simplified: in full implementation, would use Floyd-Warshall to remove redundant edges)
        }

        // Remove obsolete median vertices
        if (remove_obsolete_vertices(graph)) {
            changed = true;
        }
    }

    return graph;
}

} // namespace network
} // namespace haplokit
