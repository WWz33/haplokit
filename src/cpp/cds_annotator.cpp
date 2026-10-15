#include "cds_annotator.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace haplokit {

namespace {

// NCBI translation table 1 (standard genetic code), ordered TTT..GGG with
// bases T,C,A,G. Stop codons are '*'.
const char kCodonTable[64] = {
    'F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W',
    'L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R',
    'I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R',
    'V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G',
};

int codon_index(char b1, char b2, char b3) {
    const auto base_code = [](char b) -> int {
        switch (b) {
            case 'T': case 'U': return 0;
            case 'C': return 1;
            case 'A': return 2;
            case 'G': return 3;
            default: return -1;
        }
    };
    const int c1 = base_code(b1), c2 = base_code(b2), c3 = base_code(b3);
    if (c1 < 0 || c2 < 0 || c3 < 0) return -1;
    return c1 * 16 + c2 * 4 + c3;
}

char translate_codon(char b1, char b2, char b3) {
    const int idx = codon_index(b1, b2, b3);
    return idx < 0 ? 'X' : kCodonTable[idx];
}

char complement(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'U': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default: return 'N';
    }
}

// Indexed by amino acid letter ('A'..'Z'); unknown letters map to Xaa.
const char* kAaNames[26] = {
    "Ala", // A
    "Xaa", // B
    "Cys", // C
    "Asp", // D
    "Glu", // E
    "Phe", // F
    "Gly", // G
    "His", // H
    "Ile", // I
    "Xaa", // J
    "Lys", // K
    "Leu", // L
    "Met", // M
    "Asn", // N
    "Xaa", // O
    "Pro", // P
    "Gln", // Q
    "Arg", // R
    "Ser", // S
    "Thr", // T
    "Xaa", // U
    "Val", // V
    "Trp", // W
    "Xaa", // X
    "Tyr", // Y
    "Xaa", // Z
};

std::string aa_name(char aa) {
    if (aa >= 'A' && aa <= 'Z') return kAaNames[aa - 'A'];
    return "Xaa";
}

// Header format: >ID<space or tab description...>. ID must match the GFF
// transcript (mRNA) ID or the CDS Parent value.
std::unordered_map<std::string, std::string> read_fasta(const std::string& path, std::string* error) {
    std::unordered_map<std::string, std::string> seqs;
    std::ifstream in(path);
    if (!in) {
        if (error) *error = "cannot open CDS FASTA: " + path;
        return seqs;
    }
    std::string line;
    std::string current_id;
    std::string current_seq;
    auto flush = [&]() {
        if (!current_id.empty() && !current_seq.empty()) {
            seqs.emplace(current_id, std::move(current_seq));
        }
        current_id.clear();
        current_seq.clear();
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>' || line[0] == ';') {
            flush();
            if (line[0] == '>') {
                std::string header = line.substr(1);
                const auto sep = header.find_first_of(" \t");
                current_id = sep == std::string::npos ? header : header.substr(0, sep);
            }
            continue;
        }
        current_seq += line;
    }
    flush();
    return seqs;
}

std::string to_upper_seq(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return out;
}

}  // namespace

bool CdsAnnotator::load(const std::string& gff_path, const std::string& cds_fasta_path) {
    cds_seqs_ = read_fasta(cds_fasta_path, &last_error_);
    if (cds_seqs_.empty()) {
        if (last_error_.empty()) last_error_ = "CDS FASTA has no sequences: " + cds_fasta_path;
        return false;
    }

    try {
        index_ = gffsub::AnnotationIndex::from_gff3(gff_path);
    } catch (const std::exception& e) {
        last_error_ = e.what();
        index_.reset();
        return false;
    }
    return true;
}

std::vector<CdsVariantAnnotation> CdsAnnotator::annotate_snp(
    const std::string& chrom, int64_t pos, const std::string& ref, const std::string& alt) const {
    std::vector<CdsVariantAnnotation> hits;
    if (!index_ || ref.size() != 1 || alt.size() != 1) return hits;

    for (const auto& rec : index_->overlap(chrom, pos, pos)) {
        if (rec.type != "CDS") continue;
        const std::string transcript_id(rec.parent_id ? *rec.parent_id : std::string());
        if (transcript_id.empty()) continue;
        const auto seq_it = cds_seqs_.find(transcript_id);
        if (seq_it == cds_seqs_.end()) continue;
        const std::string& cds = seq_it->second;

        // Collect the transcript's CDS records via children traversal.
        std::vector<gffsub::GffRecord> cds_records;
        for (const auto& child : index_->children_of(transcript_id)) {
            if (child.type == "CDS") cds_records.push_back(child);
        }
        if (cds_records.empty()) continue;
        std::sort(cds_records.begin(), cds_records.end(),
                  [](const gffsub::GffRecord& a, const gffsub::GffRecord& b) {
                      return a.start < b.start;
                  });

        // Map genomic pos to the CDS offset (1-based). Plus strand: walk spans
        // in ascending genomic order, offset counts up from each span's start.
        // Minus strand: the spliced CDS starts at the highest genomic span and
        // runs backwards, so walk spans in descending order and count down from
        // each span's end.
        const bool minus = rec.strand == '-';
        int64_t cds_pos = 0;
        bool in_cds = false;
        if (!minus) {
            for (const auto& span : cds_records) {
                if (pos >= span.start && pos <= span.end) {
                    cds_pos += (pos - span.start) + 1;
                    in_cds = true;
                    break;
                }
                cds_pos += span.end - span.start + 1;
            }
        } else {
            for (auto it = cds_records.rbegin(); it != cds_records.rend(); ++it) {
                if (pos >= it->start && pos <= it->end) {
                    cds_pos += (it->end - pos) + 1;
                    in_cds = true;
                    break;
                }
                cds_pos += it->end - it->start + 1;
            }
        }
        if (!in_cds) continue;

        CdsVariantAnnotation ann;
        ann.transcript_id = transcript_id;
        // gene id: the transcript's parent, when present.
        const auto tx = index_->find_by_id(transcript_id);
        if (tx && tx->parent_id) ann.gene_id = *tx->parent_id;

        const char ref_upper = static_cast<char>(std::toupper(ref[0]));
        const char alt_upper = static_cast<char>(std::toupper(alt[0]));
        ann.ref_base = minus ? complement(ref_upper) : ref_upper;
        ann.alt_base = minus ? complement(alt_upper) : alt_upper;
        ann.cds_pos = static_cast<int>(cds_pos);

        // Cross-check the FASTA base against the VCF REF; a mismatch means the
        // GFF and the FASTA describe different annotation versions.
        if (cds_pos <= 0 || static_cast<std::size_t>(cds_pos) > cds.size()) continue;
        const char fasta_base = cds[cds_pos - 1];
        if (fasta_base != ann.ref_base) continue;

        const int codon_index_1based = (static_cast<int>(cds_pos) + 2) / 3;
        const int codon_offset = (static_cast<int>(cds_pos) - 1) % 3 + 1;
        const std::size_t codon_start = static_cast<std::size_t>(codon_index_1based - 1) * 3;
        if (codon_start + 3 > cds.size()) continue;
        ann.codon_index = codon_index_1based;
        ann.codon_offset = codon_offset;
        ann.codon_ref = cds.substr(codon_start, 3);
        ann.codon_alt = ann.codon_ref;
        ann.codon_alt[codon_offset - 1] = ann.alt_base;

        ann.aa_ref = std::string(1, translate_codon(ann.codon_ref[0], ann.codon_ref[1], ann.codon_ref[2]));
        ann.aa_alt = std::string(1, translate_codon(ann.codon_alt[0], ann.codon_alt[1], ann.codon_alt[2]));

        // Effect names follow SnpEff's Sequence Ontology output; the protein
        // change strings follow its HGVS conventions (Ter for stop,
        // p.Met1? for a lost start, p.Ter5LysextTer? for a lost stop).
        if (ann.aa_ref == ann.aa_alt) {
            if (ann.aa_ref == "*") {
                ann.effect = "stop_retained_variant";
                ann.protein_change = "p.Ter" + std::to_string(codon_index_1based) + "=";
            } else {
                ann.effect = "synonymous_variant";
                ann.protein_change = "p.(=)";
            }
        } else if (ann.aa_ref == "*") {
            ann.effect = "stop_lost";
            ann.protein_change = "p.Ter" + std::to_string(codon_index_1based) + aa_name(ann.aa_alt[0]) + "extTer?";
        } else if (ann.aa_alt == "*") {
            ann.effect = "stop_gained";
            ann.protein_change = "p." + aa_name(ann.aa_ref[0]) + std::to_string(codon_index_1based) + "Ter";
        } else if (codon_index_1based == 1 && ann.codon_ref == "ATG") {
            ann.effect = "start_lost";
            ann.protein_change = "p." + aa_name(ann.aa_ref[0]) + "1?";
        } else {
            ann.effect = "missense_variant";
            ann.protein_change = "p." + aa_name(ann.aa_ref[0]) + std::to_string(codon_index_1based) + aa_name(ann.aa_alt[0]);
        }

        hits.push_back(std::move(ann));
    }

    std::sort(hits.begin(), hits.end(),
              [](const CdsVariantAnnotation& a, const CdsVariantAnnotation& b) {
                  return a.transcript_id < b.transcript_id;
              });
    return hits;
}

}  // namespace haplokit
