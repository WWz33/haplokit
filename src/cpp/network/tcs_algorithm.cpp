#include "tcs_algorithm.h"
#include <algorithm>
#include <map>
#include <set>

namespace haplokit {
namespace network {

std::string TCSAlgorithm::generate_intermediate(const std::string& seq1,
                                                 const std::string& seq2,
                                                 int step, int total_steps) {
    // Haplotype state strings carry one '|'-separated token per variant site,
    // so swap whole tokens: a character-wise walk would emit impossible
    // genotypes such as "1/0" midway between "0/0" and "1/1".
    if (seq1.find('|') != std::string::npos || seq2.find('|') != std::string::npos) {
        auto tokens1 = split_state_tokens(seq1);
        const auto tokens2 = split_state_tokens(seq2);
        std::vector<std::size_t> diff_positions;
        for (std::size_t i = 0; i < tokens1.size() && i < tokens2.size(); ++i) {
            if (tokens1[i] != tokens2[i]) {
                diff_positions.push_back(i);
            }
        }
        const std::size_t changes_to_make =
            (diff_positions.size() * static_cast<std::size_t>(step)) / static_cast<std::size_t>(total_steps);
        for (std::size_t i = 0; i < changes_to_make && i < diff_positions.size(); ++i) {
            tokens1[diff_positions[i]] = tokens2[diff_positions[i]];
        }
        std::string joined;
        for (std::size_t i = 0; i < tokens1.size(); ++i) {
            if (i > 0) {
                joined += '|';
            }
            joined += tokens1[i];
        }
        return joined;
    }

    std::string result = seq1;
    std::vector<int> diff_positions;

    // Find all differing positions
    for (size_t i = 0; i < std::min(seq1.size(), seq2.size()); ++i) {
        if (seq1[i] != seq2[i]) {
            diff_positions.push_back(i);
        }
    }

    // Gradually change bases from seq1 to seq2
    int changes_to_make = (diff_positions.size() * step) / total_steps;
    for (int i = 0; i < changes_to_make && i < static_cast<int>(diff_positions.size()); ++i) {
        result[diff_positions[i]] = seq2[diff_positions[i]];
    }

    return result;
}

void TCSAlgorithm::insert_intermediates(Graph& graph, int u, int v, int distance) {
    const Vertex* vert_u = graph.get_vertex(u);
    const Vertex* vert_v = graph.get_vertex(v);

    if (!vert_u || !vert_v) return;

    // add_median_vertex below can reallocate Graph::vertices_, which would
    // invalidate these pointers, so copy the sequences out first.
    const std::string seq_u = vert_u->sequence;
    const std::string seq_v = vert_v->sequence;

    int prev_id = u;

    // Create intermediate vertices
    for (int step = 1; step < distance; ++step) {
        std::string intermediate_seq = generate_intermediate(seq_u, seq_v, step, distance);

        int intermediate_id = graph.add_median_vertex(intermediate_seq);

        // Connect previous vertex to this intermediate
        graph.add_edge(prev_id, intermediate_id, 1);
        prev_id = intermediate_id;
    }

    // Connect last intermediate to target
    graph.add_edge(prev_id, v, 1);
}

Graph TCSAlgorithm::compute(const std::vector<Haplotype>& haplotypes) {
    Graph graph;
    int n = haplotypes.size();

    // Add vertices
    std::vector<int> vertex_ids;
    for (const auto& hap : haplotypes) {
        int id = graph.add_vertex(hap.sequence, hap.samples);
        vertex_ids.push_back(id);
    }

    // Track connected components
    std::vector<int> component_id(n);
    for (int i = 0; i < n; ++i) {
        component_id[i] = i;
    }

    // Group vertex pairs by distance
    std::map<int, std::vector<std::pair<int, int>>> pairs_by_distance;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int dist = HammingSIMD::hamming_distance(
                haplotypes[i].sequence,
                haplotypes[j].sequence
            );
            pairs_by_distance[dist].push_back({i, j});
        }
    }

    // Process pairs in order of increasing distance
    for (auto& [distance, pairs] : pairs_by_distance) {
        int comp_a = -1, comp_b = -1;

        for (auto [u, v] : pairs) {
            int comp_u = component_id[u];
            int comp_v = component_id[v];

            // Skip if already in same component
            if (comp_u == comp_v) continue;

            // Ensure comp_u < comp_v
            if (comp_u > comp_v) {
                std::swap(comp_u, comp_v);
                std::swap(u, v);
            }

            // First pair sets the components to merge
            if (comp_a < 0) {
                comp_a = comp_u;
                comp_b = comp_v;
            }

            // Only connect pairs between comp_a and comp_b
            if (comp_u == comp_a && comp_v == comp_b) {
                if (distance == 1) {
                    // Direct connection
                    graph.add_edge(u, v, 1);
                } else {
                    // Insert intermediate vertices
                    insert_intermediates(graph, u, v, distance);
                }
            }
        }

        // Merge components
        if (comp_a >= 0) {
            for (int i = 0; i < n; ++i) {
                if (component_id[i] == comp_b) {
                    component_id[i] = comp_a;
                } else if (component_id[i] > comp_b) {
                    component_id[i]--;
                }
            }
        }
    }

    // Remove degree-2 intermediate vertices (collapse linear paths). Collapsing
    // one vertex can drop a neighbour to degree 2, so repeat until stable.
    bool collapsed = true;
    while (collapsed) {
        collapsed = false;
        std::vector<int> to_remove;
        for (const auto& vertex : graph.vertices()) {
            if (vertex.is_median && graph.degree(vertex.id) == 2) {
                to_remove.push_back(vertex.id);
            }
        }
        for (int vid : to_remove) {
            // An earlier collapse in this batch may have taken this vertex out.
            if (graph.get_vertex(vid) == nullptr) {
                continue;
            }
            auto neighbors = graph.get_neighbors(vid);
            if (neighbors.size() != 2) {
                continue;
            }
            int n1 = neighbors[0];
            int n2 = neighbors[1];
            int w1 = graph.get_edge_weight(vid, n1);
            int w2 = graph.get_edge_weight(vid, n2);

            graph.remove_vertex(vid);

            // Add direct edge
            graph.add_edge(n1, n2, w1 + w2);
            collapsed = true;
        }
    }

    return graph;
}

} // namespace network
} // namespace haplokit
