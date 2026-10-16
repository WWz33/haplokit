#pragma once

#include <optional>
#include <map>
#include <string>
#include <vector>

#include "gff_annotator.h"
#include "vcf_reader.h"

namespace haplokit {

enum class GroupBy {
    Region,
    Site,
};

enum class OutputMode {
    Summary,
    Detail,
    Both,
};

struct ViewOptions {
    GroupBy by = GroupBy::Region;
    bool impute = false;
    std::optional<double> max_diff;
    OutputMode output_mode = OutputMode::Summary;
    std::string hap_prefix = "Hap";
    int hap_pad = 2;
    std::map<std::string, std::string> sample_populations;
};

struct PopulationBreakdownRow {
    std::string population;
    std::vector<std::string> samples;
    int count = 0;
    int total = 0;
    double frequency = 0.0;
    std::string frequency_label;
};

struct HaplotypeSummaryRow {
    std::string id;
    // Grouping key: the exact allele pattern in strict mode, the cluster label
    // ("A001") in --max-diff mode.
    std::string hap;
    // Pipe-delimited allele pattern this row represents. Equals hap in strict
    // mode; in --max-diff mode it is the cluster representative's pattern, so
    // per-site states and network distances stay meaningful.
    std::string pattern;
    std::vector<std::string> states;
    // Composite CDS consequence for this haplotype, filled only when --cds is
    // active. Format: transcript:effect:protein:dna:nmd, joined with ';' across
    // transcripts.
    std::string cds;
    std::vector<std::string> samples;
    std::vector<PopulationBreakdownRow> populations;
    int count = 0;
    int total = 0;
    double frequency = 0.0;
    std::string frequency_label;
};

struct HaplotypeDetailRow {
    std::string sample;
    std::string hap;
    std::string population;
};

struct SiteRow {
    std::string chrom;
    int pos = 0;
    std::string allele;
    // Per-site CDS annotation produced when --cds is active; empty otherwise.
    std::string cds;
};

struct ViewResult {
    std::string grouping_mode;
    std::string grouping_method;
    std::string output_mode;
    bool imputed_ref = false;
    std::optional<double> max_diff;
    std::string hap_prefix = "Hap";
    int hap_pad = 2;
    int variant_count = 0;
    // Samples resolved from the VCF (after --samples-file).
    int sample_count = 0;
    // Samples that produced a haplotype, i.e. the denominator behind the
    // per-haplotype counts. Lower than sample_count when strict grouping drops
    // heterozygous or missing calls.
    int analyzed_sample_count = 0;
    int haplotype_count = 0;
    std::vector<SiteRow> sites;
    std::vector<HaplotypeSummaryRow> haplotypes;
    std::vector<HaplotypeDetailRow> accessions;
    std::optional<GeneAnnotation> annotation;
};

ViewResult build_view_result(const RegionData& data, const ViewOptions& options);
std::string serialize_view_result_json(const ViewResult& result);
std::vector<std::string> load_sample_list(const std::string& path);
std::map<std::string, std::string> load_population_groups(const std::string& path);

}  // namespace haplokit
