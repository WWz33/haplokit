#include "gff3.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>

using namespace gffsub;

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input.gff3> [options]\n"
        << "\n"
        << "Input/Region Options:\n"
        << "  -r, --region CHR:START-END\n"
        << "      Extract features overlapping the specified genomic region.\n"
        << "      Coordinates are 1-based and inclusive (GFF format).\n"
        << "      Example: -r chr1:1000000-2000000\n"
        << "\n"
        << "  -b, --bed FILE\n"
        << "      Extract features using genomic regions from a BED file.\n"
        << "      BED files use 0-based half-open coordinates, automatically\n"
        << "      converted to 1-based for internal processing.\n"
        << "\n"
        << "Feature Filter Options:\n"
        << "  -f, --feature TYPE\n"
        << "      Filter features by type (3rd column in GFF/GTF).\n"
        << "      Examples: gene, mRNA, exon, CDS, transcript\n"
        << "\n"
        << "  -L, --longest\n"
        << "      Keep only the longest transcript isoform for each gene.\n"
        << "      Per-gene comparison (AGAT logic): if gene has CDS isoforms,\n"
        << "      only compare by CDS length; otherwise compare by exon length.\n"
        << "\n"
        << "  -@, --threads N\n"
        << "      Number of threads for parallel processing (default: 1).\n"
        << "      Currently used with --longest for multi-chromosome parallelization.\n"
        << "\n"
        << "Output Options:\n"
        << "  -t, --output-format FMT\n"
        << "      Output format. Choices: gff3, gtf2, gtf3, bed\n"
        << "      gff3  - GFF3 format (##gff-version 3)\n"
        << "      gtf2  - GTF2 format (##gtf-version 2)\n"
        << "      gtf3  - GTF3/Ensembl format (##gtf-version 2.2.1)\n"
        << "      bed   - BED format (0-based half-open coordinates)\n"
        << "      Default: gff3\n"
        << "\n"
        << "  -o, --output FILE\n"
        << "      Output file path. If not specified, writes to stdout.\n"
        << "\n"
        << "  -h, --help\n"
        << "      Display this help message.\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 -f gene\n"
        << "  " << prog << " annotation.gff3 --bed regions.bed -f exon\n"
        << "  " << prog << " annotation.gff3 --longest\n"
        << "  " << prog << " annotation.gff3 --longest -@ 6\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 -t gtf3 -o out.gtf\n";
}

int main(int argc, char* argv[]) {
    std::string region_str;
    std::string bed_file;
    std::string feature;
    bool do_longest = false;
    size_t num_threads = 6;
    std::string output_format = "gff3";
    std::string output_file;

    static struct option long_options[] = {
        {"region",        required_argument, nullptr, 'r'},
        {"bed",           required_argument, nullptr, 'b'},
        {"feature",       required_argument, nullptr, 'f'},
        {"longest",       no_argument,       nullptr, 'L'},
        {"threads",       required_argument, nullptr, '@'},
        {"output-format", required_argument, nullptr, 't'},
        {"output",        required_argument, nullptr, 'o'},
        {"help",          no_argument,       nullptr, 'h'},
        {nullptr,        0,                 nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "r:b:f:L@:t:o:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'r': region_str = optarg; break;
            case 'b': bed_file = optarg; break;
            case 'f': feature = optarg; break;
            case 'L': do_longest = true; break;
            case '@': {
                size_t t = std::stoul(optarg);
                if (t == 0) t = 1;
                if (t > 256) t = 256; // cap to prevent over-subscription
                num_threads = t;
                break;
            }
            case 't': output_format = optarg; break;
            case 'o': output_file = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }

    if (optind >= argc) {
        usage(argv[0]);
        return 1;
    }

    // Validate output format
    OutputFormat fmt = OutputFormat::GFF3;
    if (output_format == "gtf") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf2") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf3") {
        fmt = OutputFormat::GTF3;
    } else if (output_format == "bed") {
        fmt = OutputFormat::BED;
    } else if (output_format == "gff3") {
        fmt = OutputFormat::GFF3;
    } else {
        std::cerr << "Error: unknown output format " << output_format << '\n';
        std::cerr << "Supported formats: gff3, gtf2, gtf3, bed\n";
        return 1;
    }

    std::string input_file = argv[optind];

    GffData data;
    IdIndex idx;

    // Parse input file
    if (parse_file(input_file, data, idx, InputFormat::GFF3) != 0) {
        std::cerr << "Error: cannot parse " << input_file << '\n';
        return 1;
    }

    Region region{"", 0, 0};
    std::optional<Region> parsed_region;

    // Apply region filters
    if (!region_str.empty()) {
        parsed_region = parse_region(region_str);
        if (!parsed_region) {
            std::cerr << "Error: invalid region format " << region_str << '\n';
            return 1;
        }
        region = *parsed_region;
        filter_by_region(data, region);
    }

    if (!bed_file.empty()) {
        filter_by_regions_from_file(data, bed_file);
    }

    // Apply feature filters
    if (do_longest) {
        filter_longest(data, idx, feature, num_threads);
    } else if (!feature.empty()) {
        filter_by_feature(data, feature);
    }

    // Output
    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!output_file.empty()) {
        out_file.open(output_file);
        if (!out_file.is_open()) {
            std::cerr << "Error: cannot open " << output_file << '\n';
            return 1;
        }
        out = &out_file;
    }

    switch (fmt) {
        case OutputFormat::GFF3: print_gff3(*out, data); break;
        case OutputFormat::GTF2:
        case OutputFormat::GTF3: print_gtf(*out, data, fmt); break;
        case OutputFormat::BED:  print_bed(*out, data); break;
    }

    return 0;
}
