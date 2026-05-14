#ifndef HAPLOKIT_NETWORK_GRAPH_H
#define HAPLOKIT_NETWORK_GRAPH_H

#include <string>
#include <vector>
#include <unordered_map>
#include <list>
#include <memory>
#include <limits>

namespace haplokit {
namespace network {

// Forward declarations
class Edge;
class Vertex;

// Optimized edge storage using hashmap + linked list
// O(1) insertion, deletion, and lookup vs PopART's vector O(n) deletion
class EdgeStorage {
public:
    void add_edge(int u, int v, int weight);
    void remove_edge(int u, int v);
    bool has_edge(int u, int v) const;
    int get_weight(int u, int v) const;

    const std::list<Edge>& edges() const { return edges_; }
    size_t size() const { return edge_map_.size(); }
    void clear();

private:
    struct EdgeKey {
        int u, v;
        bool operator==(const EdgeKey& other) const {
            return (u == other.u && v == other.v) || (u == other.v && v == other.u);
        }
    };

    struct EdgeKeyHash {
        size_t operator()(const EdgeKey& k) const {
            // Symmetric hash for undirected edges
            int min_v = std::min(k.u, k.v);
            int max_v = std::max(k.u, k.v);
            return std::hash<int>()(min_v) ^ (std::hash<int>()(max_v) << 1);
        }
    };

    std::unordered_map<EdgeKey, std::list<Edge>::iterator, EdgeKeyHash> edge_map_;
    std::list<Edge> edges_;
};

// Edge representation
struct Edge {
    int source;
    int target;
    int weight;

    Edge(int s, int t, int w) : source(s), target(t), weight(w) {}
};

// Vertex representation
struct Vertex {
    int id;
    std::string sequence;
    std::vector<std::string> samples; // samples with this haplotype
    bool is_median;

    Vertex(int id_, const std::string& seq, const std::vector<std::string>& samp = {})
        : id(id_), sequence(seq), samples(samp), is_median(false) {}
};

// Graph class for haplotype networks
class Graph {
public:
    Graph() = default;

    // Vertex operations
    int add_vertex(const std::string& sequence, const std::vector<std::string>& samples = {});
    int add_median_vertex(const std::string& sequence);
    const Vertex* get_vertex(int id) const;
    size_t vertex_count() const { return vertices_.size(); }
    const std::vector<Vertex>& vertices() const { return vertices_; }

    // Edge operations
    void add_edge(int u, int v, int weight);
    void remove_edge(int u, int v);
    bool has_edge(int u, int v) const;
    int get_edge_weight(int u, int v) const;
    const std::list<Edge>& edges() const { return edge_storage_.edges(); }

    // Graph algorithms
    std::vector<std::vector<int>> compute_distance_matrix() const;
    void floyd_warshall_subset(const std::vector<int>& subgraph_nodes,
                               std::vector<std::vector<int>>& dist) const;
    int shortest_path_length(int u, int v) const;

    // Adjacency
    std::vector<int> get_neighbors(int v) const;
    int degree(int v) const;

    void clear();

private:
    std::vector<Vertex> vertices_;
    EdgeStorage edge_storage_;
    int next_vertex_id_ = 0;
};

} // namespace network
} // namespace haplokit

#endif // HAPLOKIT_NETWORK_GRAPH_H
