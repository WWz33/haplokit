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

// One variant of a haplotype, in VCF terms.
struct HaplotypeVariant {
    int64_t pos = 0;
    std::string ref;
    std::string alt;
};

// Consequence of a whole haplotype, obtained by applying all of its variants
// to the transcript at once (unlike CdsVariantAnnotation, which looks at one
// variant in isolation). This is what catches compound effects.
struct CdsHaplotypeAnnotation {
    std::string transcript_id;
    std::string gene_id;
    int nt_delta = 0;              // net nucleotide length change
    // synonymous_variant, missense_variant, start_lost, stop_gained,
    // stop_lost, frameshift_variant, inframe_insertion, inframe_deletion,
    // inframe_altering
    std::string effect;
    std::string protein_change;    // composite HGVS, e.g. p.[Ile5Lys;Tyr294Ter]
    std::string dna_change;        // applied changes joined with '+', e.g. c.14T>A+c.548T>A
    int first_diff_aa = 0;         // 1-based first differing amino acid, 0 = none
    int aa_ref_len = 0;
    int aa_alt_len = 0;
    int ptc_codon = 0;             // 1-based premature stop codon, 0 = none
    // NMD_sensitive when the PTC sits more than 50 nt upstream of the last
    // exon-exon junction; NMD_escape otherwise (including single-exon CDS).
    std::string nmd = "-";
    // Variants that could not be applied, with the reason.
    std::vector<std::string> skipped;
    // Set when the reference CDS itself looks wrong (no terminal stop, or an
    // internal stop before the last codon).
    std::string note;
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

    // Annotate a whole haplotype: apply every variant to each transcript's CDS
    // at once, translate, and compare with the reference protein. Returns one
    // entry per transcript that has at least one applicable variant.
    std::vector<CdsHaplotypeAnnotation> annotate_haplotype(
        const std::string& chrom, const std::vector<HaplotypeVariant>& variants) const;

    const std::string& last_error() const { return last_error_; }

private:
    struct TranscriptStructure {
        std::vector<gffsub::GffRecord> cds;     // CDS spans, genomic ascending
        std::vector<gffsub::GffRecord> exons;   // exon spans, genomic ascending
        bool minus = false;
        std::string gene_id;
        std::string protein;                    // reference translation
        std::vector<int64_t> junctions;         // exon-exon junctions inside the CDS
    };
    const TranscriptStructure* structure(const std::string& transcript_id) const;

    std::optional<gffsub::AnnotationIndex> index_;
    // transcript_id -> spliced CDS sequence (5'->3', already strand-oriented)
    std::unordered_map<std::string, std::string> cds_seqs_;
    mutable std::unordered_map<std::string, std::optional<TranscriptStructure>> structures_;
    std::string last_error_;
};

}  // namespace haplokit
