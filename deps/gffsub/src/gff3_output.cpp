#include "gff3.hpp"
#include <iostream>
#include <sstream>
#include <unordered_set>

namespace gffsub {

void print_gff3(std::ostream& out, const GffData& data) {
    out << "##gff-version 3\n";
    for (const auto& rec : data) {
        if (!rec.kept) continue;
        std::string score_str = rec.score ? std::to_string(*rec.score) : ".";
        out << rec.seqid << '\t' << rec.source << '\t' << rec.type << '\t'
            << rec.start << '\t' << rec.end << '\t' << score_str << '\t'
            << rec.strand << '\t' << rec.phase << '\t' << rec.attr_raw << '\n';
    }
}

static std::string build_gtf_attrs(const std::string& gene_id_val, const std::string& transcript_id_val) {
    std::string result;
    result.reserve(64);
    result = "gene_id \"";
    result += gene_id_val;
    result += "\"; ";
    if (!transcript_id_val.empty()) {
        result += "transcript_id \"";
        result += transcript_id_val;
        result += "\"; ";
    }
    return result;
}

void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt) {
    // GTF header per AGAT spec
    if (fmt == OutputFormat::GTF3) {
        out << "##gtf-version 2.2.1\n";
    } else {
        out << "##gtf-version 2\n";
    }

    // Build mappings
    std::unordered_map<std::string, std::string> mRNA_to_gene;
    for (const auto& rec : data) {
        if (rec.kept && rec.type == "mRNA" && rec.parent_id && rec.id) {
            mRNA_to_gene[*rec.id] = *rec.parent_id;
        }
    }

    // GTF3 feature types (static to avoid recreation on each call)
    static const std::unordered_set<std::string> gtf3_types = {
        "gene", "transcript", "exon", "CDS", "start_codon", "stop_codon",
        "five_prime_utr", "three_prime_utr", "Selenocysteine", "mRNA"
    };

    for (const auto& rec : data) {
        if (!rec.kept) continue;

        // Filter by GTF3 feature types if using GTF3
        if (fmt == OutputFormat::GTF3) {
            if (gtf3_types.count(rec.type) == 0) continue;
        }

        std::string score_str = rec.score ? std::to_string(*rec.score) : ".";

        std::string gene_id_val;
        std::string transcript_id_val;

        if (rec.type == "gene") {
            gene_id_val = rec.id ? *rec.id : (rec.gene_id ? *rec.gene_id : "");
        } else if (rec.type == "mRNA" || rec.type == "transcript") {
            if (rec.parent_id) {
                gene_id_val = *rec.parent_id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
            transcript_id_val = rec.id ? *rec.id : (rec.transcript_id ? *rec.transcript_id : "");
        } else {
            // Child features
            if (rec.parent_id && mRNA_to_gene.count(*rec.parent_id)) {
                gene_id_val = mRNA_to_gene[*rec.parent_id];
                transcript_id_val = *rec.parent_id;
            } else if (rec.parent_id) {
                gene_id_val = *rec.parent_id;
                transcript_id_val = *rec.parent_id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
        }

        if (gene_id_val.empty()) continue;

        std::string gtf_type = rec.type;
        if (rec.type == "mRNA") {
            gtf_type = (fmt == OutputFormat::GTF3) ? "transcript" : "mRNA";
        }

        std::string attrs = build_gtf_attrs(gene_id_val, transcript_id_val);

        out << rec.seqid << '\t' << rec.source << '\t' << gtf_type << '\t'
            << rec.start << '\t' << rec.end << '\t' << score_str << '\t'
            << rec.strand << '\t' << rec.phase << '\t' << attrs << '\n';
    }
}

void print_bed(std::ostream& out, const GffData& data) {
    // BED is 0-based half-open, GFF is 1-based inclusive
    // BED start = GFF start - 1, BED end = GFF end
    for (const auto& rec : data) {
        if (!rec.kept) continue;

        std::string name = rec.id ? *rec.id : rec.type;
        std::string score_str = rec.score ? std::to_string(*rec.score) : "0";

        out << rec.seqid << '\t'
            << (rec.start - 1) << '\t'
            << rec.end << '\t'
            << name << '\t'
            << score_str << '\t'
            << rec.strand << '\n';
    }
}

}  // namespace gffsub
