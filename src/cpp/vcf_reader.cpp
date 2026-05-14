#include "vcf_reader.h"

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
}

namespace haplokit {
namespace {

std::vector<std::string> resolve_samples(
    bcf_hdr_t* hdr,
    const std::vector<std::string>& requested,
    std::vector<int>* sample_indexes) {
    std::vector<std::string> resolved;
    if (sample_indexes == nullptr) {
        throw std::runtime_error("sample index output is required");
    }

    if (requested.empty()) {
        const int sample_count = bcf_hdr_nsamples(hdr);
        resolved.reserve(sample_count);
        sample_indexes->reserve(sample_count);
        for (int idx = 0; idx < sample_count; ++idx) {
            resolved.emplace_back(hdr->samples[idx]);
            sample_indexes->push_back(idx);
        }
        return resolved;
    }

    resolved.reserve(requested.size());
    sample_indexes->reserve(requested.size());
    for (const auto& sample : requested) {
        const int idx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, sample.c_str());
        if (idx < 0) {
            throw std::runtime_error("sample not found in VCF header: " + sample);
        }
        resolved.push_back(sample);
        sample_indexes->push_back(idx);
    }
    return resolved;
}

VariantRecord build_variant_record(
    bcf_hdr_t* hdr,
    bcf1_t* rec,
    const std::vector<int>& sample_indexes,
    int32_t** gt_buffer,
    int* gt_buffer_len) {
    if (hdr == nullptr || rec == nullptr || gt_buffer == nullptr || gt_buffer_len == nullptr) {
        throw std::runtime_error("cannot build variant record from null input");
    }

    bcf_unpack(rec, BCF_UN_ALL);

    VariantRecord row;
    row.chrom = bcf_hdr_id2name(hdr, rec->rid);
    row.pos = rec->pos + 1;
    row.ref = rec->d.allele[0];
    row.alt.reserve(rec->n_allele > 0 ? rec->n_allele - 1 : 0);
    for (int allele_index = 1; allele_index < rec->n_allele; ++allele_index) {
        row.alt.emplace_back(rec->d.allele[allele_index]);
    }

    const int ngt = bcf_get_genotypes(hdr, rec, gt_buffer, gt_buffer_len);
    const int total_samples = bcf_hdr_nsamples(hdr);
    if (total_samples < 0) {
        throw std::runtime_error("invalid sample count in VCF header");
    }
    if (total_samples > 0 && ngt < 0) {
        throw std::runtime_error("failed to decode GT field");
    }

    const int ploidy = total_samples == 0 ? 0 : ngt / total_samples;
    row.genotypes.reserve(sample_indexes.size());
    for (const int sample_index : sample_indexes) {
        GenotypeCall call;
        if (sample_index < 0 || sample_index >= total_samples || ploidy < 2) {
            row.genotypes.push_back(call);
            continue;
        }

        const int base = sample_index * ploidy;
        if (base + 1 >= ngt) {
            row.genotypes.push_back(call);
            continue;
        }

        const int32_t allele1 = (*gt_buffer)[base];
        const int32_t allele2 = (*gt_buffer)[base + 1];
        if (bcf_gt_is_missing(allele1) || bcf_gt_is_missing(allele2)) {
            row.genotypes.push_back(call);
            continue;
        }

        call.allele1 = bcf_gt_allele(allele1);
        call.allele2 = bcf_gt_allele(allele2);
        call.phased = bcf_gt_is_phased(allele1) || bcf_gt_is_phased(allele2);
        row.genotypes.push_back(call);
    }

    return row;
}

}  // namespace

VcfReader::VcfReader(std::string path) : path_(std::move(path)) {}

RegionData VcfReader::fetch(const Region& region, const std::vector<std::string>& samples) const {
    RegionData data;
    bcf_srs_t* sr = bcf_sr_init();
    int32_t* gt_buffer = nullptr;
    int gt_buffer_len = 0;

    if (sr == nullptr) {
        throw std::runtime_error("failed to initialize synced bcf reader");
    }

    try {
        const std::string query =
            region.chrom + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
        if (bcf_sr_set_regions(sr, query.c_str(), 0) != 0) {
            throw std::runtime_error("failed to set region query: " + query);
        }

        bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
        if (!bcf_sr_add_reader(sr, path_.c_str())) {
            throw std::runtime_error("failed to add indexed VCF/BCF reader: " + std::string(bcf_sr_strerror(sr->errnum)));
        }

        bcf_hdr_t* hdr = bcf_sr_get_header(sr, 0);
        if (hdr == nullptr) {
            throw std::runtime_error("failed to read VCF header");
        }

        std::vector<int> sample_indexes;
        data.samples = resolve_samples(hdr, samples, &sample_indexes);

        // Enable sample subset at the reader level for BCF performance
        if (!samples.empty()) {
            std::string sample_list;
            for (std::size_t i = 0; i < samples.size(); ++i) {
                if (i > 0) sample_list += ",";
                sample_list += samples[i];
            }
            bcf_sr_set_samples(sr, sample_list.c_str(), 0);
        }

        while (bcf_sr_next_line(sr) > 0) {
            bcf1_t* rec = bcf_sr_get_line(sr, 0);
            if (rec == nullptr) {
                continue;
            }
            data.variants.push_back(build_variant_record(hdr, rec, sample_indexes, &gt_buffer, &gt_buffer_len));
        }
    } catch (...) {
        std::free(gt_buffer);
        bcf_sr_destroy(sr);
        throw;
    }

    std::free(gt_buffer);
    bcf_sr_destroy(sr);
    return data;
}

}  // namespace haplokit
