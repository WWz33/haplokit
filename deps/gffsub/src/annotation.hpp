#ifndef ANNOTATION_HPP
#define ANNOTATION_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace gffsub {

struct GffRecord {
    std::string seqid;
    std::string source;
    std::string type;
    int64_t start;
    int64_t end;
    std::optional<double> score;
    char strand;
    char phase;
    std::string attr_raw;
    std::optional<std::string> id;
    std::optional<std::string> parent_id;
    std::optional<std::string> gene_id;
    std::optional<std::string> transcript_id;
    int line_idx;
    bool kept;
};

struct GeneModel {
    GffRecord gene;
    std::vector<GffRecord> records;
};

class GffData {
public:
    std::vector<GffRecord> records;

    void append(const GffRecord& rec) { records.push_back(rec); }
    auto size() const { return records.size(); }
    auto begin() { return records.begin(); }
    auto end() { return records.end(); }
    auto begin() const { return records.begin(); }
    auto end() const { return records.end(); }
    void clear() { records.clear(); }
    void reserve(size_t n) { records.reserve(n); }
};

class IdIndex {
public:
    std::unordered_map<std::string, std::vector<int>> index;

    void add(const std::string& id, int idx) {
        index[id].push_back(idx);
    }

    std::optional<int> lookup(const std::string& id) const {
        auto it = index.find(id);
        if (it != index.end() && !it->second.empty()) {
            return it->second.front();
        }
        return std::nullopt;
    }

    void clear() { index.clear(); }
};

class AnnotationIndex {
public:
    static AnnotationIndex from_gff3(const std::string& path);

    std::optional<GffRecord> find_by_id(std::string_view id) const;
    std::optional<GffRecord> find_gene(std::string_view id) const;
    std::vector<GffRecord> parents_of(std::string_view id) const;
    std::vector<GffRecord> children_of(std::string_view parent_id) const;
    std::vector<GffRecord> descendants_of(std::string_view parent_id) const;
    std::vector<GffRecord> overlap(std::string_view seqid, int64_t start, int64_t end) const;
    std::optional<GffRecord> nearest_gene(std::string_view seqid, int64_t start, int64_t end) const;
    std::vector<GffRecord> with_attribute(std::string_view key, std::string_view value) const;
    std::optional<GeneModel> gene_model(std::string_view id) const;

private:
    GffData data_;
    std::unordered_map<std::string, int> id_to_record_;
    std::unordered_map<std::string, std::vector<int>> gene_lookup_;
    std::unordered_map<std::string, std::vector<int>> parents_by_child_id_;
    std::unordered_map<std::string, std::vector<int>> children_by_parent_id_;

    explicit AnnotationIndex(GffData data);
};

}  // namespace gffsub

#endif  // ANNOTATION_HPP
