#include "graph.h"
#include <algorithm>
#include <queue>
#include <stdexcept>

namespace haplokit {
namespace network {

// EdgeStorage implementation
void EdgeStorage::add_edge(int u, int v, int weight) {
    EdgeKey key{u, v};
    if (edge_map_.find(key) != edge_map_.end()) {
        return; // edge already exists
    }

    edges_.emplace_back(u, v, weight);
    edge_map_[key] = std::prev(edges_.end());
}

void EdgeStorage::remove_edge(int u, int v) {
    EdgeKey key{u, v};
    auto it = edge_map_.find(key);
    if (it != edge_map_.end()) {
        edges_.erase(it->second);
        edge_map_.erase(it);
    }
}

bool EdgeStorage::has_edge(int u, int v) const {
    EdgeKey key{u, v};
    return edge_map_.find(key) != edge_map_.end();
}

int EdgeStorage::get_weight(int u, int v) const {
    EdgeKey key{u, v};
    auto it = edge_map_.find(key);
    if (it != edge_map_.end()) {
        return it->second->weight;
    }
    return -1;
}

void EdgeStorage::clear() {
    edges_.clear();
    edge_map_.clear();
}

// Graph implementation
int Graph::add_vertex(const std::string& sequence, const std::vector<std::string>& samples) {
    int id = next_vertex_id_++;
    vertices_.emplace_back(id, sequence, samples);
    return id;
}

int Graph::add_median_vertex(const std::string& sequence) {
    int id = next_vertex_id_++;
    vertices_.emplace_back(id, sequence);
    vertices_.back().is_median = true;
    return id;
}

const Vertex* Graph::get_vertex(int id) const {
    for (const auto& v : vertices_) {
        if (v.id == id) return &v;
    }
    return nullptr;
}

void Graph::add_edge(int u, int v, int weight) {
    edge_storage_.add_edge(u, v, weight);
}

void Graph::remove_edge(int u, int v) {
    edge_storage_.remove_edge(u, v);
}

bool Graph::has_edge(int u, int v) const {
    return edge_storage_.has_edge(u, v);
}

int Graph::get_edge_weight(int u, int v) const {
    return edge_storage_.get_weight(u, v);
}

std::vector<int> Graph::get_neighbors(int v) const {
    std::vector<int> neighbors;
    for (const auto& edge : edge_storage_.edges()) {
        if (edge.source == v) {
            neighbors.push_back(edge.target);
        } else if (edge.target == v) {
            neighbors.push_back(edge.source);
        }
    }
    return neighbors;
}

int Graph::degree(int v) const {
    return get_neighbors(v).size();
}

std::vector<std::vector<int>> Graph::compute_distance_matrix() const {
    int n = vertices_.size();
    const int INF = std::numeric_limits<int>::max() / 2;
    std::vector<std::vector<int>> dist(n, std::vector<int>(n, INF));

    // Initialize diagonal
    for (int i = 0; i < n; ++i) {
        dist[i][i] = 0;
    }

    // Initialize edges
    for (const auto& edge : edge_storage_.edges()) {
        dist[edge.source][edge.target] = edge.weight;
        dist[edge.target][edge.source] = edge.weight;
    }

    // Floyd-Warshall
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }

    return dist;
}

void Graph::floyd_warshall_subset(const std::vector<int>& subgraph_nodes,
                                   std::vector<std::vector<int>>& dist) const {
    // Subset Floyd-Warshall: O(k³) where k = subgraph size << n
    for (int k : subgraph_nodes) {
        for (int i : subgraph_nodes) {
            for (int j : subgraph_nodes) {
                if (dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
}

int Graph::shortest_path_length(int u, int v) const {
    // BFS for unweighted shortest path
    if (u == v) return 0;

    std::queue<std::pair<int, int>> q;
    std::unordered_map<int, bool> visited;

    q.push({u, 0});
    visited[u] = true;

    while (!q.empty()) {
        auto [curr, dist] = q.front();
        q.pop();

        for (int neighbor : get_neighbors(curr)) {
            if (neighbor == v) return dist + 1;
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.push({neighbor, dist + 1});
            }
        }
    }

    return -1; // not connected
}

void Graph::clear() {
    vertices_.clear();
    edge_storage_.clear();
    next_vertex_id_ = 0;
}

} // namespace network
} // namespace haplokit
