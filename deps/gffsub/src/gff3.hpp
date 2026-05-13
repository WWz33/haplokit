#ifndef GFF3_HPP
#define GFF3_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <optional>
#include <unordered_map>
#include <string_view>
#include <thread>
#include <functional>

namespace gffsub {

enum class OutputFormat { GFF3, GTF2, GTF3, BED };
enum class InputFormat { GFF3, GTF, BED };

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

int parse_file(const std::string& filename, GffData& data, IdIndex& idx, InputFormat format);

struct Region {
    std::string seqid;
    int64_t start;
    int64_t end;
};

std::optional<Region> parse_region(std::string_view region_str);

void filter_by_region(GffData& data, const Region& region);
void filter_by_regions_from_file(GffData& data, const std::string& bed_file);
void filter_by_feature(GffData& data, std::string_view feature_type);
void filter_longest(GffData& data, IdIndex& idx, std::string_view feature_type, size_t num_threads = 1);

void print_gff3(std::ostream& out, const GffData& data);
void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt);
void print_bed(std::ostream& out, const GffData& data);

}  // namespace gffsub

#endif  // GFF3_HPP
