#include <iostream>
#include <optional>
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "cds_annotator.h"
#include "gff_annotator.h"
#include "selector.h"
#include "vcf_reader.h"
#include "view_backend.h"

namespace {

struct ParsedViewJsonCommand {
    std::string input_vcf;
    haplokit::Region region;
    std::vector<std::string> samples;
    haplokit::ViewOptions options;
    std::string gff3_path;
    std::string cds_fasta_path;
};

struct ParsedViewTargetsJsonCommand {
    std::string input_vcf;
    std::vector<haplokit::Region> targets;
    std::vector<std::string> samples;
    haplokit::ViewOptions options;
    std::string gff3_path;
    std::string cds_fasta_path;
};

double parse_max_diff(const std::string& value) {
    std::size_t consumed = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("max diff must be a finite number in [0,1]");
    }
    if (consumed != value.size() || !std::isfinite(parsed) || parsed < 0.0 || parsed > 1.0) {
        throw std::runtime_error("max diff must be in [0,1]");
    }
    return parsed;
}

int parse_hap_pad(const std::string& value) {
    std::size_t consumed = 0;
    int parsed = 0;
    try {
        parsed = std::stoi(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("hap pad must be a positive integer");
    }
    if (consumed != value.size() || parsed < 1) {
        throw std::runtime_error("hap pad must be a positive integer");
    }
    return parsed;
}

int64_t parse_nonnegative_window(const std::string& value, const std::string& name) {
    std::size_t consumed = 0;
    int64_t parsed = 0;
    try {
        parsed = std::stoll(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(name + " must be a non-negative integer");
    }
    if (consumed != value.size() || parsed < 0) {
        throw std::runtime_error(name + " must be a non-negative integer");
    }
    return parsed;
}

std::string trim_copy(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<haplokit::Region> parse_target_regions(const std::string& value) {
    std::vector<haplokit::Region> targets;
    std::string chrom;
    std::size_t start = 0;
    while (start <= value.size()) {
        const auto end = value.find(',', start);
        const auto token = trim_copy(end == std::string::npos
            ? value.substr(start)
            : value.substr(start, end - start));
        if (token.empty()) {
            throw std::runtime_error("target list contains an empty entry");
        }

        auto region = haplokit::parse_region(token);
        if (chrom.empty()) {
            chrom = region.chrom;
        } else if (region.chrom != chrom) {
            throw std::runtime_error("targets must be on the same chromosome");
        }
        targets.push_back(std::move(region));

        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    if (targets.empty()) {
        throw std::runtime_error("target list must not be empty");
    }
    return targets;
}

std::vector<haplokit::Region> load_target_regions_file(const std::string& path) {
    if (path == "-") {
        throw std::runtime_error("targets file does not read from stdin");
    }

    std::ifstream handle(path);
    if (!handle) {
        throw std::runtime_error("failed to open targets file: " + path);
    }

    std::vector<haplokit::Region> targets;
    std::string line;
    while (std::getline(handle, line)) {
        const auto stripped = trim_copy(line);
        if (stripped.empty() || stripped.front() == '#') {
            continue;
        }
        auto parsed = parse_target_regions(stripped);
        targets.insert(targets.end(), parsed.begin(), parsed.end());
    }
    if (targets.empty()) {
        throw std::runtime_error("targets file did not contain any target regions");
    }
    return targets;
}

haplokit::Region covering_region(const std::vector<haplokit::Region>& targets) {
    if (targets.empty()) {
        throw std::runtime_error("target list must not be empty");
    }
    haplokit::Region span = targets.front();
    for (const auto& target : targets) {
        if (target.chrom != span.chrom) {
            throw std::runtime_error("targets must be on the same chromosome");
        }
        if (target.start < span.start) {
            span.start = target.start;
        }
        if (target.end > span.end) {
            span.end = target.end;
        }
    }
    return span;
}

ParsedViewJsonCommand parse_view_json_command(int argc, char** argv) {
    if (argc < 4) {
        throw std::runtime_error(
            "usage: haplokit_cpp view-json <vcf> <region> [--by region|site] [--samples-file path] "
            "[--impute] [--output summary|detail] [--max-diff x] [--population-file path]");
    }

    ParsedViewJsonCommand parsed;
    parsed.input_vcf = argv[2];
    parsed.region = haplokit::parse_region(argv[3]);

    for (int idx = 4; idx < argc; ++idx) {
        const std::string arg = argv[idx];
        if (arg == "--by") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--by requires a value");
            }
            const std::string value = argv[++idx];
            if (value == "region") {
                parsed.options.by = haplokit::GroupBy::Region;
            } else if (value == "site") {
                parsed.options.by = haplokit::GroupBy::Site;
            } else {
                throw std::runtime_error("unsupported --by value: " + value);
            }
            continue;
        }
        if (arg == "--samples-file") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--samples-file requires a path");
            }
            parsed.samples = haplokit::load_sample_list(argv[++idx]);
            continue;
        }
        if (arg == "--population-file") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--population-file requires a path");
            }
            parsed.options.sample_populations = haplokit::load_population_groups(argv[++idx]);
            continue;
        }
        if (arg == "--impute") {
            parsed.options.impute = true;
            continue;
        }
        if (arg == "--output") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--output requires a value");
            }
            const std::string value = argv[++idx];
            if (value == "summary") {
                parsed.options.output_mode = haplokit::OutputMode::Summary;
            } else if (value == "detail") {
                parsed.options.output_mode = haplokit::OutputMode::Detail;
            } else if (value == "both") {
                parsed.options.output_mode = haplokit::OutputMode::Both;
            } else {
                throw std::runtime_error("unsupported --output value: " + value);
            }
            continue;
        }
        if (arg == "--max-diff") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--max-diff requires a value");
            }
            parsed.options.max_diff = parse_max_diff(argv[++idx]);
            continue;
        }
        if (arg == "--gff3") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--gff3 requires a path");
            }
            parsed.gff3_path = argv[++idx];
            continue;
        }
        if (arg == "--cds") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--cds requires a path");
            }
            parsed.cds_fasta_path = argv[++idx];
            continue;
        }
        if (arg == "--hap-prefix") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-prefix requires a value");
            }
            parsed.options.hap_prefix = argv[++idx];
            continue;
        }
        if (arg == "--hap-pad") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-pad requires a value");
            }
            parsed.options.hap_pad = parse_hap_pad(argv[++idx]);
            continue;
        }
        throw std::runtime_error("unsupported argument: " + arg);
    }

    if (parsed.options.by == haplokit::GroupBy::Site && !haplokit::is_site_region(parsed.region)) {
        throw std::runtime_error("--by site requires a single-position region");
    }

    return parsed;
}

haplokit::ViewOptions parse_common_view_options(
    int argc,
    char** argv,
    int start_index,
    std::vector<std::string>* samples,
    std::string* gff3_path,
    std::string* cds_fasta_path = nullptr) {
    haplokit::ViewOptions options;
    if (samples == nullptr) {
        throw std::runtime_error("samples output is required");
    }
    for (int idx = start_index; idx < argc; ++idx) {
        const std::string arg = argv[idx];
        if (arg == "--samples-file") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--samples-file requires a path");
            }
            *samples = haplokit::load_sample_list(argv[++idx]);
            continue;
        }
        if (arg == "--population-file") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--population-file requires a path");
            }
            options.sample_populations = haplokit::load_population_groups(argv[++idx]);
            continue;
        }
        if (arg == "--impute") {
            options.impute = true;
            continue;
        }
        if (arg == "--output") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--output requires a value");
            }
            const std::string value = argv[++idx];
            if (value == "summary") {
                options.output_mode = haplokit::OutputMode::Summary;
            } else if (value == "detail") {
                options.output_mode = haplokit::OutputMode::Detail;
            } else if (value == "both") {
                options.output_mode = haplokit::OutputMode::Both;
            } else {
                throw std::runtime_error("unsupported --output value: " + value);
            }
            continue;
        }
        if (arg == "--max-diff") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--max-diff requires a value");
            }
            options.max_diff = parse_max_diff(argv[++idx]);
            continue;
        }
        if (arg == "--gff3") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--gff3 requires a path");
            }
            if (gff3_path) {
                *gff3_path = argv[++idx];
            }
            continue;
        }
        if (arg == "--cds") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--cds requires a path");
            }
            if (cds_fasta_path) {
                *cds_fasta_path = argv[++idx];
            }
            continue;
        }
        if (arg == "--hap-prefix") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-prefix requires a value");
            }
            options.hap_prefix = argv[++idx];
            continue;
        }
        if (arg == "--hap-pad") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-pad requires a value");
            }
            options.hap_pad = parse_hap_pad(argv[++idx]);
            continue;
        }
        throw std::runtime_error("unsupported argument: " + arg);
    }
    return options;
}

ParsedViewTargetsJsonCommand parse_view_targets_json_command(int argc, char** argv, bool targets_file = false) {
    if (argc < 4) {
        throw std::runtime_error(
            "usage: haplokit_cpp view-targets-json <vcf> <targets> [--samples-file path] "
            "[--population-file path] [--impute] [--output summary|detail|both] [--max-diff x] [--gff3 path]");
    }

    ParsedViewTargetsJsonCommand parsed;
    parsed.input_vcf = argv[2];
    parsed.targets = targets_file ? load_target_regions_file(argv[3]) : parse_target_regions(argv[3]);
    parsed.options = parse_common_view_options(argc, argv, 4, &parsed.samples, &parsed.gff3_path, &parsed.cds_fasta_path);
    parsed.options.by = haplokit::GroupBy::Region;
    return parsed;
}

std::vector<haplokit::Region> load_bed_regions(const std::string& path) {
    std::ifstream handle(path);
    if (!handle) {
        throw std::runtime_error("failed to open regions file: " + path);
    }

    std::vector<haplokit::Region> regions;
    std::string line;
    while (std::getline(handle, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            const auto end = line.find('\t', start);
            if (end == std::string::npos) {
                fields.push_back(line.substr(start));
                break;
            }
            fields.push_back(line.substr(start, end - start));
            start = end + 1;
        }
        if (fields.size() < 3) {
            throw std::runtime_error("BED rows must have at least 3 columns");
        }

        haplokit::Region region;
        region.chrom = fields[0];
        // BED is 0-based half-open; Region and the VCF query are 1-based
        // inclusive, so the left edge shifts by one. The right edge already
        // matches (a BED end of N is exclusive, meaning 1-based base N).
        region.start = std::stoi(fields[1]) + 1;
        region.end = std::stoi(fields[2]);
        if (region.start > region.end) {
            throw std::runtime_error("BED region start must be <= end: " + line);
        }
        regions.push_back(region);
    }
    return regions;
}

// Annotate result.sites with CDS consequences when --cds is active. The
// allele field holds "REF/ALT1,ALT2,..."; each ALT is annotated against every
// overlapping CDS transcript and the results are joined with ';'.
void apply_cds_annotation(haplokit::ViewResult& result, const haplokit::RegionData& data,
                          const haplokit::CdsAnnotator& annotator) {
    if (data.variants.size() != result.sites.size()) {
        return;
    }
    for (std::size_t idx = 0; idx < result.sites.size(); ++idx) {
        const auto& variant = data.variants[idx];
        auto& site = result.sites[idx];
        std::string joined;
        for (const auto& alt : variant.alt) {
            for (const auto& ann : annotator.annotate_snp(variant.chrom, variant.pos, variant.ref, alt)) {
                if (!joined.empty()) {
                    joined += ";";
                }
                joined += ann.transcript_id + ":c." + std::to_string(ann.cds_pos)
                    + ann.ref_base + ">" + ann.alt_base
                    + ":" + ann.protein_change + "(" + ann.effect + ")";
            }
        }
        site.cds = joined;
    }
}

// Split "REF/ALT1,ALT2" into its parts.
void split_site_allele(const std::string& allele, std::string* ref, std::vector<std::string>* alts) {
    const auto slash = allele.find('/');
    if (slash == std::string::npos) {
        *ref = allele;
        return;
    }
    *ref = allele.substr(0, slash);
    std::istringstream alt_input(allele.substr(slash + 1));
    std::string alt;
    while (std::getline(alt_input, alt, ',')) {
        alts->push_back(alt);
    }
}

// Rebuild the variants carried by one haplotype from its encoded states
// ("0/0|1/1|2/2" for "ref|alt1|alt2") and the site allele strings.
std::vector<haplokit::HaplotypeVariant> haplotype_variants(const std::string& hap,
                                                          const std::vector<haplokit::SiteRow>& sites) {
    std::vector<haplokit::HaplotypeVariant> variants;
    std::size_t start = 0;
    std::size_t site_idx = 0;
    while (start <= hap.size() && site_idx < sites.size()) {
        const auto end = hap.find('|', start);
        const std::string state = hap.substr(start, end == std::string::npos ? std::string::npos : end - start);
        const auto slash = state.find('/');
        const int allele_index = slash == std::string::npos ? 0 : std::stoi(state.substr(0, slash));
        if (allele_index > 0) {
            std::string ref;
            std::vector<std::string> alts;
            split_site_allele(sites[site_idx].allele, &ref, &alts);
            const std::size_t alt_index = static_cast<std::size_t>(allele_index - 1);
            if (alt_index < alts.size()) {
                variants.push_back(haplokit::HaplotypeVariant{
                    sites[site_idx].pos, ref, alts[alt_index]});
            }
        }
        ++site_idx;
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return variants;
}

// Annotate every haplotype of the result with its composite CDS consequence.
void apply_haplotype_cds_annotation(haplokit::ViewResult& result,
                                    const haplokit::CdsAnnotator& annotator,
                                    const std::string& chrom) {
    for (auto& hap : result.haplotypes) {
        const auto variants = haplotype_variants(hap.hap, result.sites);
        if (variants.empty()) continue;
        std::string joined;
        for (const auto& ann : annotator.annotate_haplotype(chrom, variants)) {
            if (!joined.empty()) joined += ";";
            joined += ann.transcript_id + ":" + ann.effect + ":" + ann.protein_change
                + ":" + ann.dna_change + ":" + ann.nmd;
            for (const auto& skipped : ann.skipped) {
                joined += ":skip(" + skipped + ")";
            }
        }
        hap.cds = joined;
    }
}

// Load the CDS annotator when --cds is given; requires --gff3 for the CDS
// interval structure.
std::optional<haplokit::CdsAnnotator> load_cds_annotator(
    const std::string& gff3_path, const std::string& cds_fasta_path) {
    if (cds_fasta_path.empty()) {
        return std::nullopt;
    }
    if (gff3_path.empty()) {
        throw std::runtime_error("--cds requires --gff3 for the CDS interval structure");
    }
    haplokit::CdsAnnotator annotator;
    if (!annotator.load(gff3_path, cds_fasta_path)) {
        throw std::runtime_error("failed to load CDS annotation: " + annotator.last_error());
    }
    return annotator;
}

int run_view_json(int argc, char** argv) {
    const auto parsed = parse_view_json_command(argc, argv);

    std::optional<haplokit::GffAnnotator> annotator;
    if (!parsed.gff3_path.empty()) {
        haplokit::GffAnnotator ann;
        if (!ann.load(parsed.gff3_path)) {
            throw std::runtime_error("failed to parse GFF3: " + parsed.gff3_path);
        }
        annotator = std::move(ann);
    }

    auto cds_annotator = load_cds_annotator(parsed.gff3_path, parsed.cds_fasta_path);

    haplokit::VcfReader reader(parsed.input_vcf);
    const auto data = reader.fetch(parsed.region, parsed.samples);
    auto result = haplokit::build_view_result(data, parsed.options);

    if (annotator.has_value()) {
        result.annotation = annotator->annotate(parsed.region.chrom, parsed.region.start, parsed.region.end);
    }
    if (cds_annotator.has_value()) {
        apply_cds_annotation(result, data, *cds_annotator);
        apply_haplotype_cds_annotation(result, *cds_annotator, data.variants.empty() ? std::string() : data.variants.front().chrom);
    }

    std::cout << haplokit::serialize_view_result_json(result) << "\n";
    return 0;
}

int run_view_targets_json(int argc, char** argv, bool targets_file = false) {
    const auto parsed = parse_view_targets_json_command(argc, argv, targets_file);

    std::optional<haplokit::GffAnnotator> annotator;
    if (!parsed.gff3_path.empty()) {
        haplokit::GffAnnotator ann;
        if (!ann.load(parsed.gff3_path)) {
            throw std::runtime_error("failed to parse GFF3: " + parsed.gff3_path);
        }
        annotator = std::move(ann);
    }

    auto cds_annotator = load_cds_annotator(parsed.gff3_path, parsed.cds_fasta_path);

    const auto span = covering_region(parsed.targets);
    haplokit::VcfReader reader(parsed.input_vcf);
    const auto data = reader.fetch_targets(parsed.targets, parsed.samples);
    auto result = haplokit::build_view_result(data, parsed.options);

    if (annotator.has_value()) {
        result.annotation = annotator->annotate(span.chrom, span.start, span.end);
    }
    if (cds_annotator.has_value()) {
        apply_cds_annotation(result, data, *cds_annotator);
        apply_haplotype_cds_annotation(result, *cds_annotator, data.variants.empty() ? std::string() : data.variants.front().chrom);
    }

    std::cout << haplokit::serialize_view_result_json(result) << "\n";
    return 0;
}

int run_view_bed_jsonl(int argc, char** argv) {
    if (argc < 4) {
        throw std::runtime_error(
            "usage: haplokit_cpp view-bed-jsonl <vcf> <regions.bed> [--samples-file path] "
            "[--population-file path] [--impute] [--output summary|detail] [--max-diff x] [--gff3 path]");
    }

    const std::string input_vcf = argv[2];
    const std::string regions_file = argv[3];
    std::vector<std::string> samples;
    std::string gff3_path;
    std::string cds_fasta_path;
    auto options = parse_common_view_options(argc, argv, 4, &samples, &gff3_path, &cds_fasta_path);
    options.by = haplokit::GroupBy::Region;

    std::optional<haplokit::GffAnnotator> annotator;
    if (!gff3_path.empty()) {
        haplokit::GffAnnotator ann;
        if (!ann.load(gff3_path)) {
            throw std::runtime_error("failed to parse GFF3: " + gff3_path);
        }
        annotator = std::move(ann);
    }

    auto cds_annotator = load_cds_annotator(gff3_path, cds_fasta_path);

    haplokit::VcfReader reader(input_vcf);
    const auto regions = load_bed_regions(regions_file);
    for (std::size_t idx = 0; idx < regions.size(); ++idx) {
        const auto data = reader.fetch(regions[idx], samples);
        auto result = haplokit::build_view_result(data, options);

        if (cds_annotator.has_value()) {
            apply_cds_annotation(result, data, *cds_annotator);
        apply_haplotype_cds_annotation(result, *cds_annotator, data.variants.empty() ? std::string() : data.variants.front().chrom);
        }
        if (annotator.has_value()) {
            // Region coordinates are already 1-based inclusive (load_bed_regions
            // converts the BED left edge), which is what GFF annotation expects.
            result.annotation = annotator->annotate(regions[idx].chrom, regions[idx].start, regions[idx].end);
        }

        std::cout << haplokit::serialize_view_result_json(result) << "\n";
    }
    return 0;
}

int run_debug_fetch(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: haplokit_cpp <vcf> <region>\n";
        return 1;
    }

    const haplokit::Region region = haplokit::parse_region(argv[2]);
    haplokit::VcfReader reader(argv[1]);
    const auto data = reader.fetch(region);
    std::cout << "opened " << argv[1] << " for " << region.chrom << ":" << region.start << "-" << region.end << "\n";
    std::cout << "samples " << data.samples.size() << "\n";
    std::cout << "variants " << data.variants.size() << "\n";
    return 0;
}

int run_resolve_gene(int argc, char** argv) {
    if (argc < 4) {
        throw std::runtime_error(
            "usage: haplokit_cpp resolve-gene <gff3> <gene-id> "
            "[--upstream n] [--downstream n] [--strand-aware]");
    }

    haplokit::GffAnnotator annotator;
    const std::string gff3_path = argv[2];
    if (!annotator.load(gff3_path)) {
        throw std::runtime_error("failed to parse GFF3: " + gff3_path);
    }

    const std::string gene_id = argv[3];
    int64_t upstream = 0;
    int64_t downstream = 0;
    bool strand_aware = false;
    for (int idx = 4; idx < argc; ++idx) {
        const std::string arg = argv[idx];
        if (arg == "--upstream") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--upstream requires a value");
            }
            upstream = parse_nonnegative_window(argv[++idx], "upstream");
            continue;
        }
        if (arg == "--downstream") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--downstream requires a value");
            }
            downstream = parse_nonnegative_window(argv[++idx], "downstream");
            continue;
        }
        if (arg == "--strand-aware") {
            strand_aware = true;
            continue;
        }
        throw std::runtime_error("unsupported resolve-gene argument: " + arg);
    }

    const auto gene = annotator.find_gene_window(gene_id, upstream, downstream, strand_aware);
    if (!gene) {
        throw std::runtime_error("gene ID not found in GFF3: " + gene_id);
    }

    std::cout << gene->seqid << ":" << gene->start << "-" << gene->end << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 2 && std::string(argv[1]) == "resolve-gene") {
            return run_resolve_gene(argc, argv);
        }
        if (argc >= 2 && std::string(argv[1]) == "view-json") {
            return run_view_json(argc, argv);
        }
        if (argc >= 2 && std::string(argv[1]) == "view-targets-json") {
            return run_view_targets_json(argc, argv);
        }
        if (argc >= 2 && std::string(argv[1]) == "view-targets-file-json") {
            return run_view_targets_json(argc, argv, true);
        }
        if (argc >= 2 && std::string(argv[1]) == "view-bed-jsonl") {
            return run_view_bed_jsonl(argc, argv);
        }
        return run_debug_fetch(argc, argv);
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return 2;
    }
}
