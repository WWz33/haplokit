#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "gff3.hpp"

namespace haplokit {

struct CdsVariantAnnotation {
    std::string transcript_id;   // GFF mRNA ID, e.g. test1G0387.1
    std::string gene_id;         // parent gene ID when resolvable
    int cds_pos = 0;             // 1-based CDS coordinate of the SNP
    char ref_base = 'N';         // base on the transcript strand
    char alt_base = 'N';         // alt base on the transcript strand
    std::string codon_ref;       // reference codon (3 nt)
    std::string codon_alt;       // alternate codon (3 nt)
    int codon_index = 0;         // 1-based codon number
    int codon_offset = 0;        // 1..3 position within the codon
    std::string aa_ref;          // reference amino acid (one-letter)
    std::string aa_alt;          // alternate amino acid (one-letter)
    // SnpEff-style consequence: synonymous_variant, missense_variant,
    // stop_gained, stop_lost, start_lost.
    std::string effect;
    // HGVS-style protein change, e.g. p.Ala159Val; p.(=) when synonymous.
    std::string protein_change;
};

class CdsAnnotator {
public:
    // gff_path supplies the CDS interval structure; cds_fasta_path holds one
    // spliced CDS sequence per transcript (header = transcript ID). Both must
    // describe the same annotation.
    bool load(const std::string& gff_path, const std::string& cds_fasta_path);

    // Annotate a single SNV (ref/alt length 1) against every CDS transcript
    // overlapping (chrom, pos). Returns one entry per hit, ordered by
    // transcript ID. Non-SNP variants get an empty vector.
    std::vector<CdsVariantAnnotation> annotate_snp(
        const std::string& chrom, int64_t pos, const std::string& ref, const std::string& alt) const;

    const std::string& last_error() const { return last_error_; }

private:
    std::optional<gffsub::AnnotationIndex> index_;
    // transcript_id -> spliced CDS sequence (5'->3', already strand-oriented)
    std::unordered_map<std::string, std::string> cds_seqs_;
    std::string last_error_;
};

}  // namespace haplokit
