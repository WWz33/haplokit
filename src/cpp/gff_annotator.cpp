#include "gff_annotator.h"
#include "gff3.hpp"

#include <algorithm>
#include <limits>

namespace haptools {

bool GffAnnotator::load(const std::string& gff_path) {
    gffsub::GffData data;
    gffsub::IdIndex idx;
    if (gffsub::parse_file(gff_path, data, idx, gffsub::InputFormat::GFF3) != 0) {
        return false;
    }

    for (const auto& rec : data) {
        if (rec.type == "gene" && rec.id) {
            genes_[rec.seqid].push_back(GeneRecord{*rec.id, rec.start, rec.end, rec.strand});
        }
    }

    for (auto& [chrom, gene_list] : genes_) {
        std::sort(gene_list.begin(), gene_list.end(),
                  [](const GeneRecord& a, const GeneRecord& b) { return a.start < b.start; });
    }

    return true;
}

GeneAnnotation GffAnnotator::annotate(const std::string& chrom, int64_t start, int64_t end) const {
    GeneAnnotation ann;

    auto it = genes_.find(chrom);
    if (it == genes_.end() || it->second.empty()) {
        ann.mode = "none";
        return ann;
    }

    const auto& gene_list = it->second;

    for (const auto& gene : gene_list) {
        if (gene.start > end) break;
        if (gene.end >= start) {
            ann.mode = "overlap";
            ann.id = gene.id;
            ann.seqid = chrom;
            ann.start = gene.start;
            ann.end = gene.end;
            ann.strand = gene.strand;
            return ann;
        }
    }

    const GeneRecord* nearest = nullptr;
    int64_t nearest_dist = std::numeric_limits<int64_t>::max();
    for (const auto& gene : gene_list) {
        int64_t dist;
        if (end < gene.start) {
            dist = gene.start - end;
        } else if (start > gene.end) {
            dist = start - gene.end;
        } else {
            dist = 0;
        }
        if (dist < nearest_dist) {
            nearest = &gene;
            nearest_dist = dist;
        }
    }

    if (nearest) {
        ann.mode = "nearest";
        ann.id = nearest->id;
        ann.seqid = chrom;
        ann.start = nearest->start;
        ann.end = nearest->end;
        ann.strand = nearest->strand;
        ann.distance = nearest_dist;
        return ann;
    }

    ann.mode = "none";
    return ann;
}

}  // namespace haptools
