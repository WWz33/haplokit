#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace haptools {

struct GeneAnnotation {
    std::string mode;  // "overlap", "nearest", "none"
    std::string id;
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
    char strand = '.';
    std::optional<int64_t> distance;
};

class GffAnnotator {
public:
    bool load(const std::string& gff_path);
    GeneAnnotation annotate(const std::string& chrom, int64_t start, int64_t end) const;

private:
    struct GeneRecord {
        std::string id;
        int64_t start;
        int64_t end;
        char strand;
    };
    std::unordered_map<std::string, std::vector<GeneRecord>> genes_;
};

}  // namespace haptools
