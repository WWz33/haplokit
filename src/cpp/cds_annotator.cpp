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

char complement(char b);

std::string reverse_complement(const std::string& seq) {
    std::string out;
    out.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        out.push_back(complement(*it));
    }
    return out;
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

std::string translate_seq(const std::string& cds) {
    std::string protein;
    protein.reserve(cds.size() / 3);
    for (std::size_t i = 0; i + 3 <= cds.size(); i += 3) {
        protein.push_back(translate_codon(cds[i], cds[i + 1], cds[i + 2]));
    }
    return protein;
}

// HGVS DNA change for one variant at a 1-based CDS coordinate. VCF alleles
// carry an anchor base, so insertions report alt[1:] and deletions report
// ref[1:] after the anchor.
std::string format_dna_change(int cds_pos, const std::string& ref, const std::string& alt) {
    const int ref_len = static_cast<int>(ref.size());
    const int alt_len = static_cast<int>(alt.size());
    if (ref_len == alt_len) {
        return "c." + std::to_string(cds_pos) + ref + ">" + alt;
    }
    if (alt_len > ref_len) {
        // c.14_15insA for an insertion anchored at c.14
        return "c." + std::to_string(cds_pos) + "_" + std::to_string(cds_pos + 1) + "ins"
            + alt.substr(1);
    }
    // c.15_17delAC for a deletion of AC after the anchor at c.14
    return "c." + std::to_string(cds_pos + 1) + "_" + std::to_string(cds_pos + ref_len - 1)
        + "del" + ref.substr(1);
}

std::string aa_three(char aa) { return aa == '*' ? std::string("Ter") : aa_name(aa); }

// Substitutions only (CDS length unchanged): list every changed codon.
// Bounded to kMaxListed entries so a haplotype with many variants stays
// readable; the full count is preserved by the trailing marker.
std::string format_substitution_list(const std::string& ref_prot, const std::string& alt_prot) {
    constexpr std::size_t kMaxListed = 12;
    std::vector<std::string> parts;
    std::size_t total = 0;
    for (std::size_t i = 0; i < ref_prot.size() && i < alt_prot.size(); ++i) {
        if (ref_prot[i] == alt_prot[i]) continue;
        ++total;
        if (parts.size() < kMaxListed) {
            parts.push_back(aa_three(ref_prot[i]) + std::to_string(i + 1) + aa_three(alt_prot[i]));
        }
    }
    if (parts.empty()) return "p.(=)";
    std::string joined;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) joined += ";";
        joined += parts[i];
    }
    if (total > parts.size()) joined += ";...";
    if (parts.size() == 1 && total == 1) return "p." + joined;
    return "p.[" + joined + "]";
}

std::string inserted_residues(const std::string& mid) {
    constexpr std::size_t kMaxListed = 6;
    if (mid.size() > kMaxListed) return std::to_string(mid.size()) + "aa";
    std::string out;
    for (char aa : mid) out += aa_three(aa);
    return out;
}

// Frame-preserving indel: report the changed reference span. The residue
// sequence is listed for short insertions and summarised by length otherwise,
// so the string stays bounded no matter how many variants the haplotype has.
std::string format_indel_change(const std::string& ref_prot, const std::string& alt_prot) {
    std::size_t prefix = 0;
    const std::size_t shorter = std::min(ref_prot.size(), alt_prot.size());
    while (prefix < shorter && ref_prot[prefix] == alt_prot[prefix]) ++prefix;
    std::size_t suffix = 0;
    while (suffix < shorter - prefix
           && ref_prot[ref_prot.size() - 1 - suffix] == alt_prot[alt_prot.size() - 1 - suffix]) {
        ++suffix;
    }
    const std::string ref_mid = ref_prot.substr(prefix, ref_prot.size() - prefix - suffix);
    const std::string alt_mid = alt_prot.substr(prefix, alt_prot.size() - prefix - suffix);
    const int ref_beg = static_cast<int>(prefix) + 1;
    const int ref_end = static_cast<int>(ref_prot.size() - suffix);

    if (ref_mid.empty()) {
        const int left = ref_beg - 1;
        if (left < 1) return "p.0_" + aa_three(ref_prot[0]) + "1ins" + inserted_residues(alt_mid);
        return "p." + aa_three(ref_prot[left - 1]) + std::to_string(left) + "_"
            + aa_three(ref_prot[left]) + std::to_string(left + 1) + "ins"
            + inserted_residues(alt_mid);
    }
    if (alt_mid.empty()) {
        return "p." + aa_three(ref_mid.front()) + std::to_string(ref_beg) + "_"
            + aa_three(ref_mid.back()) + std::to_string(ref_end) + "del";
    }
    // Trimming can leave the inserted/deleted residues attached to either end of
    // the changed reference span. Recognise that so a plain insertion or
    // deletion is not reported as a delins.
    if (alt_mid.size() > ref_mid.size()) {
        const std::size_t extra = alt_mid.size() - ref_mid.size();
        if (alt_mid.compare(0, ref_mid.size(), ref_mid) == 0) {
            const std::string inserted = alt_mid.substr(ref_mid.size());
            const int left = ref_end;
            if (left >= 1 && left < static_cast<int>(ref_prot.size())) {
                return "p." + aa_three(ref_prot[left - 1]) + std::to_string(left) + "_"
                    + aa_three(ref_prot[left]) + std::to_string(left + 1) + "ins"
                    + inserted_residues(inserted);
            }
        }
        if (alt_mid.compare(extra, ref_mid.size(), ref_mid) == 0) {
            const std::string inserted = alt_mid.substr(0, extra);
            const int right = ref_beg;
            if (right >= 2 && right <= static_cast<int>(ref_prot.size())) {
                return "p." + aa_three(ref_prot[right - 2]) + std::to_string(right - 1) + "_"
                    + aa_three(ref_prot[right - 1]) + std::to_string(right) + "ins"
                    + inserted_residues(inserted);
            }
        }
    } else if (ref_mid.size() > alt_mid.size()) {
        const std::size_t extra = ref_mid.size() - alt_mid.size();
        if (ref_mid.compare(0, alt_mid.size(), alt_mid) == 0) {
            const std::string deleted = ref_mid.substr(alt_mid.size());
            const int beg = ref_beg + static_cast<int>(alt_mid.size());
            return "p." + aa_three(deleted.front()) + std::to_string(beg) + "_"
                + aa_three(deleted.back()) + std::to_string(ref_end) + "del";
        }
        if (ref_mid.compare(extra, alt_mid.size(), alt_mid) == 0) {
            const std::string deleted = ref_mid.substr(0, extra);
            return "p." + aa_three(deleted.front()) + std::to_string(ref_beg) + "_"
                + aa_three(deleted.back()) + std::to_string(ref_beg + static_cast<int>(extra) - 1)
                + "del";
        }
    }
    return "p." + aa_three(ref_mid.front()) + std::to_string(ref_beg) + "_"
        + aa_three(ref_mid.back()) + std::to_string(ref_end) + "delins"
        + inserted_residues(alt_mid);
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

const CdsAnnotator::TranscriptStructure* CdsAnnotator::structure(
    const std::string& transcript_id) const {
    const auto cached = structures_.find(transcript_id);
    if (cached != structures_.end()) {
        return cached->second ? &*cached->second : nullptr;
    }
    std::optional<TranscriptStructure> built;
    const auto seq_it = cds_seqs_.find(transcript_id);
    if (index_ && seq_it != cds_seqs_.end()) {
        TranscriptStructure st;
        for (const auto& child : index_->children_of(transcript_id)) {
            if (child.type == "CDS") st.cds.push_back(child);
            else if (child.type == "exon") st.exons.push_back(child);
        }
        if (!st.cds.empty()) {
            std::sort(st.cds.begin(), st.cds.end(),
                      [](const gffsub::GffRecord& a, const gffsub::GffRecord& b) {
                          return a.start < b.start;
                      });
            std::sort(st.exons.begin(), st.exons.end(),
                      [](const gffsub::GffRecord& a, const gffsub::GffRecord& b) {
                          return a.start < b.start;
                      });
            st.minus = st.cds.front().strand == '-';
            const auto tx = index_->find_by_id(transcript_id);
            if (tx && tx->parent_id) st.gene_id = *tx->parent_id;
            st.protein = translate_seq(seq_it->second);
            // Exon-exon junctions that fall inside the CDS: the end of every
            // CDS span except the last (spans are genomic-ascending, and the
            // spliced CDS runs through them in that order on the plus strand
            // and in reverse on the minus strand, so the set of boundary
            // coordinates is the same either way).
            int64_t cum = 0;
            for (std::size_t i = 0; i + 1 < st.cds.size(); ++i) {
                cum += st.cds[i].end - st.cds[i].start + 1;
                st.junctions.push_back(cum);
            }
            built = std::move(st);
        }
    }
    auto& slot = structures_[transcript_id];
    slot = std::move(built);
    return slot ? &*slot : nullptr;
}

std::vector<CdsHaplotypeAnnotation> CdsAnnotator::annotate_haplotype(
    const std::string& chrom, const std::vector<HaplotypeVariant>& variants) const {
    std::vector<CdsHaplotypeAnnotation> out;
    if (!index_ || variants.empty()) return out;

    // Transcripts whose CDS is touched by at least one variant.
    std::vector<std::string> transcript_ids;
    for (const auto& variant : variants) {
        for (const auto& rec : index_->overlap(chrom, variant.pos, variant.pos + static_cast<int64_t>(variant.ref.size()) - 1)) {
            if (rec.type != "CDS" || !rec.parent_id) continue;
            if (std::find(transcript_ids.begin(), transcript_ids.end(), *rec.parent_id) == transcript_ids.end()) {
                transcript_ids.push_back(*rec.parent_id);
            }
        }
    }
    std::sort(transcript_ids.begin(), transcript_ids.end());

    for (const auto& transcript_id : transcript_ids) {
        const TranscriptStructure* st = structure(transcript_id);
        if (st == nullptr) continue;

        // Map every variant into CDS coordinates and apply it to a copy of the
        // reference CDS. Variants are applied right to left so that earlier
        // coordinates stay valid as indels shift the sequence.
        struct Applied {
            int cds_pos = 0;
            std::string ref;
            std::string alt;
            std::size_t order = 0;
        };
        std::vector<Applied> todo;
        CdsHaplotypeAnnotation ann;
        ann.transcript_id = transcript_id;
        ann.gene_id = st->gene_id;
        ann.aa_ref_len = static_cast<int>(st->protein.size());

        std::size_t order = 0;
        for (const auto& variant : variants) {
            ++order;
            int64_t cds_pos = 0;
            bool in_cds = false;
            int64_t span_upper = 0;   // genomic end of the span holding the variant
            if (!st->minus) {
                for (const auto& span : st->cds) {
                    if (variant.pos >= span.start && variant.pos <= span.end) {
                        cds_pos += (variant.pos - span.start) + 1;
                        span_upper = span.end;
                        in_cds = true;
                        break;
                    }
                    cds_pos += span.end - span.start + 1;
                }
            } else {
                for (auto it = st->cds.rbegin(); it != st->cds.rend(); ++it) {
                    if (variant.pos >= it->start && variant.pos <= it->end) {
                        cds_pos += (it->end - variant.pos) + 1;
                        span_upper = it->end;
                        in_cds = true;
                        break;
                    }
                    cds_pos += it->end - it->start + 1;
                }
            }
            if (!in_cds) continue;   // variant sits outside this transcript's CDS

            std::string ref = variant.ref;
            std::string alt = variant.alt;
            if (st->minus) {
                ref = reverse_complement(ref);
                alt = reverse_complement(alt);
            }
            // A REF that runs past the end of the containing CDS span crosses
            // an exon junction, where the genomic sequence is not contiguous in
            // the spliced CDS.
            const bool spans_junction =
                variant.pos + static_cast<int64_t>(variant.ref.size()) - 1 > span_upper;
            if (spans_junction) {
                ann.skipped.push_back("c." + std::to_string(cds_pos) + " spans exon junction");
                continue;
            }
            todo.push_back(Applied{static_cast<int>(cds_pos), ref, alt, order});
        }

        // Keep the entry when variants were rejected so the reason is visible;
        // drop it only when nothing about this transcript's CDS is affected.
        if (todo.empty() && ann.skipped.empty()) continue;

        const std::string ref_cds = cds_seqs_.at(transcript_id);
        std::string cds = ref_cds;
        std::sort(todo.begin(), todo.end(), [](const Applied& a, const Applied& b) {
            if (a.cds_pos != b.cds_pos) return a.cds_pos > b.cds_pos;
            return a.order > b.order;
        });

        int nt_delta = 0;
        int n_ins = 0;
        int n_del = 0;
        bool has_indel = false;
        std::vector<Applied> applied;
        for (const auto& item : todo) {
            const std::size_t beg = static_cast<std::size_t>(item.cds_pos) - 1;
            if (beg + item.ref.size() > cds.size() || cds.compare(beg, item.ref.size(), item.ref) != 0) {
                ann.skipped.push_back("c." + std::to_string(item.cds_pos)
                                      + " REF does not match the CDS sequence");
                continue;
            }
            cds.replace(beg, item.ref.size(), item.alt);
            nt_delta += static_cast<int>(item.alt.size()) - static_cast<int>(item.ref.size());
            if (item.alt.size() > item.ref.size()) ++n_ins;
            else if (item.alt.size() < item.ref.size()) ++n_del;
            if (item.alt.size() != item.ref.size()) has_indel = true;
            applied.push_back(item);
        }
        if (applied.empty()) {
            // Every variant was rejected: keep the entry so the reasons are
            // reported, but leave the consequence fields empty.
            ann.effect = "-";
            out.push_back(std::move(ann));
            continue;
        }

        // Report DNA changes in ascending CDS order.
        std::sort(applied.begin(), applied.end(), [](const Applied& a, const Applied& b) {
            return a.cds_pos < b.cds_pos;
        });
        for (const auto& item : applied) {
            if (!ann.dna_change.empty()) ann.dna_change += "+";
            ann.dna_change += format_dna_change(item.cds_pos, item.ref, item.alt);
        }

        // Only the region up to the first stop codon is protein; codons past it
        // are not translated, so compare the truncated sequences. This also
        // keeps a mis-annotated reference (internal stop) from producing
        // nonsense calls for variants that sit after it.
        const auto truncate_at_stop = [](const std::string& protein) {
            const auto stop = protein.find('*');
            return stop == std::string::npos ? protein : protein.substr(0, stop + 1);
        };
        const auto protein_length = [](const std::string& protein) {
            return static_cast<int>(protein.size()) - (protein.empty() || protein.back() != '*' ? 0 : 1);
        };
        const std::string ref_prot = truncate_at_stop(st->protein);
        const std::string alt_prot = truncate_at_stop(translate_seq(cds));
        ann.nt_delta = nt_delta;
        ann.aa_ref_len = protein_length(ref_prot);
        ann.aa_alt_len = protein_length(alt_prot);

        if (ref_prot.find('*') == std::string::npos) {
            ann.note = "reference CDS has no stop codon";
        } else if (ref_prot.size() != st->protein.size()) {
            ann.note = "reference CDS has an internal stop at codon "
                + std::to_string(ref_prot.size());
        }

        int first_diff = 0;
        const std::size_t comparable = std::min(ref_prot.size(), alt_prot.size());
        for (std::size_t i = 0; i < comparable; ++i) {
            if (ref_prot[i] != alt_prot[i]) {
                first_diff = static_cast<int>(i) + 1;
                break;
            }
        }
        if (first_diff == 0 && ref_prot.size() != alt_prot.size()) {
            first_diff = static_cast<int>(comparable) + 1;
        }
        ann.first_diff_aa = first_diff;

        const std::size_t ref_stop = ref_prot.find('*');
        const std::size_t alt_stop = alt_prot.find('*');
        const bool premature_stop = alt_stop != std::string::npos
            && (ref_stop == std::string::npos || alt_stop < ref_stop);
        ann.ptc_codon = premature_stop ? static_cast<int>(alt_stop) + 1 : 0;

        if (ref_prot == alt_prot) {
            ann.effect = "synonymous_variant";
        } else if (first_diff == 1 && ref_cds.compare(0, 3, "ATG") == 0) {
            // The reference start codon is gone (csq calls this start_lost
            // whether or not the change is also a frameshift).
            ann.effect = "start_lost";
        } else if (nt_delta % 3 != 0) {
            ann.effect = "frameshift_variant";
        } else if (has_indel) {
            // Frame preserved, but the haplotype carries indels: an insertion
            // and a deletion that cancel out still counts as inframe_altering.
            ann.effect = (n_ins > 0 && n_del > 0) ? "inframe_altering"
                : (nt_delta > 0 ? "inframe_insertion" : "inframe_deletion");
        } else if (premature_stop) {
            ann.effect = "stop_gained";
        } else if (ref_stop != std::string::npos
                   && (alt_stop == std::string::npos || alt_stop > ref_stop)) {
            ann.effect = "stop_lost";
        } else {
            ann.effect = "missense_variant";
        }

        // Protein change in SnpEff/HGVS style.
        if (ann.effect == "synonymous_variant") {
            ann.protein_change = "p.(=)";
        } else if (ann.effect == "start_lost") {
            ann.protein_change = "p." + aa_three(ref_prot[0]) + "1?";
        } else if (ann.effect == "stop_gained") {
            ann.protein_change = "p." + aa_three(ref_prot[ann.ptc_codon - 1])
                + std::to_string(ann.ptc_codon) + "Ter";
        } else if (ann.effect == "frameshift_variant") {
            const int aa_pos = first_diff > 0 ? first_diff : 1;
            const int distance = ann.ptc_codon > 0 ? ann.ptc_codon - aa_pos + 1 : 0;
            ann.protein_change = "p." + aa_three(ref_prot[aa_pos - 1]) + std::to_string(aa_pos)
                + "fs" + (ann.ptc_codon > 0 ? "Ter" + std::to_string(distance) : "?");
        } else if (ann.effect == "stop_lost") {
            const int stop_codon = static_cast<int>(ref_stop) + 1;
            const std::string new_aa = alt_stop != std::string::npos && alt_stop < alt_prot.size()
                ? aa_three(alt_prot[ref_stop]) : std::string("?");
            ann.protein_change = "p.Ter" + std::to_string(stop_codon) + new_aa + "extTer?";
        } else if (!has_indel) {
            ann.protein_change = format_substitution_list(ref_prot, alt_prot);
        } else {
            ann.protein_change = format_indel_change(ref_prot, alt_prot);
        }

        // NMD: a PTC more than 50 nt upstream of the last exon-exon junction
        // triggers nonsense-mediated decay.
        if (ann.ptc_codon > 0 && !st->junctions.empty()) {
            const int64_t ptc_cds = (static_cast<int64_t>(ann.ptc_codon) - 1) * 3 + 1;
            const int64_t last_junction = st->junctions.back();
            ann.nmd = (last_junction - ptc_cds) > 50 ? "NMD_sensitive" : "NMD_escape";
        } else if (ann.ptc_codon > 0) {
            ann.nmd = "NMD_escape";   // single-exon CDS: no junction to trigger NMD
        }

        out.push_back(std::move(ann));
    }
    return out;
}

}  // namespace haplokit
