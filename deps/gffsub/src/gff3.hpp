#ifndef GFF3_HPP
#define GFF3_HPP

#include "annotation.hpp"

#include <ostream>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace gffsub {

enum class OutputFormat { GFF3, GTF2, GTF3, BED };
enum class InputFormat { GFF3, GTF, BED };

int parse_file(const std::string& filename, GffData& data, IdIndex& idx, InputFormat format);
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);

struct Region {
    std::string seqid;
    int64_t start;
    int64_t end;
};

struct BedRegion {
    std::string seqid;
    int64_t start;
    int64_t end;
};

std::optional<Region> parse_region(std::string_view region_str);
BedRegion to_bed_region(const GffRecord& rec);
Region from_bed_region(const BedRegion& region);
Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware);

void filter_by_region(GffData& data, const Region& region);
void filter_by_regions_from_file(GffData& data, const std::string& bed_file);
void filter_by_feature(GffData& data, std::string_view feature_type);
void filter_by_seqid(GffData& data, std::string_view seqid);
void filter_by_source(GffData& data, std::string_view source);
void filter_by_score(GffData& data, std::optional<double> score);
void filter_by_strand(GffData& data, char strand);
void filter_by_phase(GffData& data, char phase);
void filter_longest_isoform(GffData& data, IdIndex& idx, std::string_view feature_type, size_t num_threads = 1);
void filter_longest(GffData& data, IdIndex& idx, std::string_view feature_type, size_t num_threads = 1);

void print_gff3(std::ostream& out, const GffData& data);
void print_gtf3(std::ostream& out, const GffData& data);
void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt);
void print_bed(std::ostream& out, const GffData& data);

}  // namespace gffsub

#endif  // GFF3_HPP
