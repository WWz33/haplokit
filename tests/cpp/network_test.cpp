#include "graph.h"
#include "hamming_simd.h"
#include "haplotype.h"
#include "tcs_algorithm.h"
#include "msn_algorithm.h"
#include "mjn_algorithm.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <string>

using namespace haplokit::network;

void test_hamming_distance() {
    std::cout << "Testing Hamming distance...\n";

    std::string seq1 = "ATCGATCG";
    std::string seq2 = "ATCGATCG";
    assert(HammingSIMD::hamming_distance(seq1, seq2) == 0);

    seq2 = "ATCGTTCG";
    assert(HammingSIMD::hamming_distance(seq1, seq2) == 1);

    seq2 = "TTTTTTTT";
    assert(HammingSIMD::hamming_distance(seq1, seq2) == 6);

    std::cout << "  ✓ Hamming distance tests passed\n";
}

void test_graph_operations() {
    std::cout << "Testing graph operations...\n";

    Graph graph;

    // Add vertices
    int v1 = graph.add_vertex("ATCG", {"S1", "S2"});
    int v2 = graph.add_vertex("ATCG", {"S3"});
    int v3 = graph.add_median_vertex("ATCC");

    assert(graph.vertex_count() == 3);
    assert(graph.get_vertex(v1)->samples.size() == 2);
    assert(graph.get_vertex(v3)->is_median == true);

    // Add edges
    graph.add_edge(v1, v2, 1);
    graph.add_edge(v2, v3, 2);

    assert(graph.has_edge(v1, v2));
    assert(graph.has_edge(v2, v3));
    assert(!graph.has_edge(v1, v3));
    assert(graph.get_edge_weight(v1, v2) == 1);

    // Test neighbors
    auto neighbors = graph.get_neighbors(v2);
    assert(neighbors.size() == 2);

    // Remove edge
    graph.remove_edge(v1, v2);
    assert(!graph.has_edge(v1, v2));

    std::cout << "  ✓ Graph operations tests passed\n";
}

void test_msn_algorithm() {
    std::cout << "Testing MSN algorithm...\n";

    std::vector<Haplotype> haplotypes = {
        {"ATCGATCG", {"S1", "S2"}},
        {"ATCGTTCG", {"S3"}},
        {"TTCGATCG", {"S4"}},
        {"ATCGATCC", {"S5"}}
    };

    Graph graph = MSNAlgorithm::compute(haplotypes, 0);

    // MSN should connect all vertices
    assert(graph.vertex_count() == 4);

    // MSN should have n-1 edges for n vertices (minimum spanning tree)
    assert(graph.edges().size() >= 3);

    std::cout << "  ✓ MSN algorithm tests passed\n";
}

void test_tcs_algorithm() {
    std::cout << "Testing TCS algorithm...\n";

    std::vector<Haplotype> haplotypes = {
        {"ATCG", {"S1"}},
        {"ATCC", {"S2"}},
        {"TTCG", {"S3"}}
    };

    Graph graph = TCSAlgorithm::compute(haplotypes);

    // TCS should create network with original vertices
    assert(graph.vertex_count() >= 3);

    std::cout << "  ✓ TCS algorithm tests passed\n";
}

void test_mjn_algorithm() {
    std::cout << "Testing MJN algorithm...\n";

    std::vector<Haplotype> haplotypes = {
        {"ATCGATCG", {"S1"}},
        {"ATCGTTCG", {"S2"}},
        {"TTCGTTCG", {"S3"}}
    };

    Graph graph = MJNAlgorithm::compute(haplotypes, 0);

    // MJN should create network with original vertices
    // May add median vertices
    assert(graph.vertex_count() >= 3);

    std::cout << "  ✓ MJN algorithm tests passed\n";
}

void test_edge_storage() {
    std::cout << "Testing EdgeStorage (hashmap + linked list)...\n";

    EdgeStorage storage;

    // Add edges
    storage.add_edge(0, 1, 5);
    storage.add_edge(1, 2, 3);
    storage.add_edge(0, 2, 8);

    assert(storage.size() == 3);
    assert(storage.has_edge(0, 1));
    assert(storage.has_edge(1, 0)); // symmetric
    assert(storage.get_weight(0, 1) == 5);

    // Remove edge - O(1) operation
    storage.remove_edge(0, 1);
    assert(!storage.has_edge(0, 1));
    assert(storage.size() == 2);

    // Clear
    storage.clear();
    assert(storage.size() == 0);

    std::cout << "  ✓ EdgeStorage tests passed\n";
}

void test_simd_availability() {
    std::cout << "Testing SIMD availability...\n";

    bool has_avx2 = HammingSIMD::has_avx2_support();
    std::cout << "  AVX2 support: " << (has_avx2 ? "YES" : "NO") << "\n";

    // Test packed sequence encoding
    std::string seq = "ATCGATCGATCGATCG";
    auto packed = HammingSIMD::pack_sequence(seq);
    assert(packed.size() == 1); // 16 bases fit in 1 uint32_t

    std::cout << "  ✓ SIMD tests passed\n";
}

int main() {
    std::cout << "=== Haplokit Network Algorithm Tests ===\n\n";

    try {
        test_hamming_distance();
        test_graph_operations();
        test_edge_storage();
        test_simd_availability();
        test_msn_algorithm();
        test_tcs_algorithm();
        test_mjn_algorithm();

        std::cout << "\n✓ All tests passed!\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n✗ Test failed: " << e.what() << "\n";
        return 1;
    }
}
