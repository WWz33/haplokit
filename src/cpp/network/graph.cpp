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

void Graph::remove_vertex(int id) {
    const auto neighbors = get_neighbors(id);
    for (int neighbor : neighbors) {
        edge_storage_.remove_edge(id, neighbor);
    }
    vertices_.erase(
        std::remove_if(vertices_.begin(), vertices_.end(),
                       [id](const Vertex& v) { return v.id == id; }),
        vertices_.end());
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

} // namespace network
} // namespace haplokit
