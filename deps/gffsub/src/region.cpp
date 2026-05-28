#include "gff3.hpp"

namespace gffsub {

std::optional<Region> parse_region(std::string_view region_str) {
    size_t colon = region_str.find(':');
    if (colon == std::string_view::npos) return std::nullopt;

    std::string seqid(region_str.substr(0, colon));
    auto range_part = region_str.substr(colon + 1);

    size_t dash = range_part.find('-');
    if (dash == std::string_view::npos) return std::nullopt;

    int64_t start = 0;
    int64_t end = 0;
    try {
        start = std::stoll(std::string(range_part.substr(0, dash)));
        end = std::stoll(std::string(range_part.substr(dash + 1)));
    } catch (const std::exception&) {
        return std::nullopt;
    }

    return Region{seqid, start, end};
}

BedRegion to_bed_region(const GffRecord& rec) {
    return BedRegion{rec.seqid, rec.start - 1, rec.end};
}

Region from_bed_region(const BedRegion& region) {
    return Region{region.seqid, region.start + 1, region.end};
}

Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware) {
    int64_t left_extension = upstream;
    int64_t right_extension = downstream;
    if (strand_aware && rec.strand == '-') {
        left_extension = downstream;
        right_extension = upstream;
    }

    int64_t start = rec.start - left_extension;
    if (start < 1) {
        start = 1;
    }
    return Region{rec.seqid, start, rec.end + right_extension};
}

}  // namespace gffsub
