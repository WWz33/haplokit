#include "gff3.hpp"
#include <algorithm>
#include <future>
#include <unordered_map>

namespace gffsub {

void filter_longest(GffData& data, IdIndex& /*idx*/, std::string_view feature_type, size_t num_threads) {
    // Build gene -> [mRNA indices] index once
    std::unordered_map<std::string, std::vector<int>> gene_to_mRNAs;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        if (rec.type == "mRNA" && rec.parent_id) {
            gene_to_mRNAs[*rec.parent_id].push_back(i);
        }
    }

    // Build mRNA -> [child indices] index once
    std::unordered_map<std::string, std::vector<int>> mrna_to_children;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        if (rec.parent_id) {
            mrna_to_children[*rec.parent_id].push_back(i);
        }
    }

    // Group genes by chromosome
    std::unordered_map<std::string, std::vector<int>> chrom_to_gene_idx;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        if (rec.type == "gene" && rec.id) {
            chrom_to_gene_idx[rec.seqid].push_back(i);
        }
    }

    auto process_chromosome = [&](const std::string& chrom, const std::vector<int>& gene_indices) {
        (void)chrom; // reserved for future debugging
        for (int gene_idx : gene_indices) {
            const auto& gene = data.records[gene_idx];
            if (!gene.id) continue;

            // Find mRNAs for this gene using index
            auto mRNA_it = gene_to_mRNAs.find(*gene.id);
            if (mRNA_it == gene_to_mRNAs.end()) {
                // Orphan gene with no mRNA children - mark as not kept (AGAT behavior)
                data.records[gene_idx].kept = false;
                continue;
            }
            if (mRNA_it->second.size() <= 1) continue;

            const auto& mRNA_indices = mRNA_it->second;

            // Per-gene check: does ANY mRNA have CDS?
            bool gene_has_cds = false;
            for (int mrna_idx : mRNA_indices) {
                const auto& mrna = data.records[mrna_idx];
                if (!mrna.id) continue;
                auto child_it = mrna_to_children.find(*mrna.id);
                if (child_it != mrna_to_children.end()) {
                    for (int child_idx : child_it->second) {
                        if (data.records[child_idx].type == "CDS") {
                            gene_has_cds = true;
                            break;
                        }
                    }
                }
                if (gene_has_cds) break;
            }

            // Find longest mRNA
            int longest_idx = -1;
            int64_t max_len = -1;

            for (int mrna_idx : mRNA_indices) {
                const auto& mrna = data.records[mrna_idx];
                if (!mrna.id) continue;

                auto child_it = mrna_to_children.find(*mrna.id);
                if (child_it == mrna_to_children.end()) continue;

                int64_t len = 0;
                bool found = false;

                if (gene_has_cds) {
                    for (int child_idx : child_it->second) {
                        const auto& child = data.records[child_idx];
                        if (child.type == "CDS") {
                            len += child.end - child.start + 1;
                            found = true;
                        }
                    }
                    if (!found) continue; // mRNA without CDS is skipped
                } else {
                    for (int child_idx : child_it->second) {
                        const auto& child = data.records[child_idx];
                        if (child.type == "exon") {
                            len += child.end - child.start + 1;
                            found = true;
                        }
                    }
                    if (!found) {
                        len = mrna.end - mrna.start + 1;
                    }
                }

                if (len > max_len) {
                    max_len = len;
                    longest_idx = mrna_idx;
                }
            }

            // Mark longest as kept, others as not kept
            if (longest_idx >= 0) {
                for (int mrna_idx : mRNA_indices) {
                    data.records[mrna_idx].kept = (mrna_idx == longest_idx);
                }
                // Mark children of longest mRNA as kept, others as not kept
                for (int mrna_idx : mRNA_indices) {
                    const auto& mrna = data.records[mrna_idx];
                    if (!mrna.id) continue;
                    auto child_it = mrna_to_children.find(*mrna.id);
                    if (child_it != mrna_to_children.end()) {
                        bool is_longest = (mrna_idx == longest_idx);
                        for (int child_idx : child_it->second) {
                            data.records[child_idx].kept = is_longest;
                        }
                    }
                }
            }
        }
    };

    if (num_threads <= 1) {
        for (auto& [chrom, gene_indices] : chrom_to_gene_idx) {
            process_chromosome(chrom, gene_indices);
        }
    } else {
        std::vector<std::future<void>> futures;
        for (auto& [chrom, gene_indices] : chrom_to_gene_idx) {
            auto f = std::async(std::launch::async, [&, chrom, gene_indices]() {
                process_chromosome(chrom, gene_indices);
            });
            futures.push_back(std::move(f));
        }
        for (auto& f : futures) {
            f.get();
        }
    }

    if (!feature_type.empty()) {
        filter_by_feature(data, feature_type);
    }
}

}  // namespace gffsub
