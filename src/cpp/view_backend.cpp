#include "view_backend.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <nlohmann/json.hpp>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace haptools {
namespace {

std::optional<std::string> normalize_call(const GenotypeCall& call, bool impute_ref) {
    if (call.allele1 < 0 || call.allele2 < 0) {
        if (impute_ref) {
            return std::string("0/0");
        }
        return std::nullopt;
    }
    if (call.allele1 != call.allele2) {
        return std::nullopt;
    }
    if (call.allele1 == 0) {
        return std::string("0/0");
    }
    return std::to_string(call.allele1) + "/" + std::to_string(call.allele2);
}

std::string join_states(const std::vector<std::string>& states) {
    std::ostringstream buffer;
    for (std::size_t idx = 0; idx < states.size(); ++idx) {
        if (idx > 0) {
            buffer << "|";
        }
        buffer << states[idx];
    }
    return buffer.str();
}

std::vector<HaplotypeDetailRow> build_sample_profiles(const RegionData& data, bool impute_ref) {
    std::vector<HaplotypeDetailRow> rows;
    for (std::size_t sample_idx = 0; sample_idx < data.samples.size(); ++sample_idx) {
        std::vector<std::string> states;
        states.reserve(data.variants.size());
        bool skip = false;
        for (const auto& variant : data.variants) {
            const auto normalized = normalize_call(variant.genotypes[sample_idx], impute_ref);
            if (!normalized.has_value()) {
                skip = true;
                break;
            }
            states.push_back(*normalized);
        }
        if (skip) {
            continue;
        }
        rows.push_back(HaplotypeDetailRow{data.samples[sample_idx], join_states(states)});
    }
    return rows;
}

std::vector<HaplotypeDetailRow> build_approx_accessions(
    const RegionData& data,
    bool impute_ref,
    double max_diff) {
    auto raw_profiles = build_sample_profiles(data, impute_ref);

    std::vector<HaplotypeDetailRow> rows;
    std::vector<std::vector<std::string>> representatives;
    std::vector<std::string> labels;
    for (const auto& profile : raw_profiles) {
        const auto& sample_name = profile.sample;
        const auto& states_str = profile.hap;
        // Parse the joined state string back for distance calc
        std::vector<std::string> states;
        std::istringstream iss(states_str);
        std::string token;
        while (std::getline(iss, token, '|')) {
            states.push_back(token);
        }

        std::optional<std::string> assigned;
        for (std::size_t rep_idx = 0; rep_idx < representatives.size(); ++rep_idx) {
            const auto& representative = representatives[rep_idx];
            const auto comparable = std::max(representative.size(), states.size());
            if (comparable == 0) {
                continue;
            }
            int diffs = 0;
            for (std::size_t state_idx = 0; state_idx < representative.size() && state_idx < states.size(); ++state_idx) {
                if (representative[state_idx] != states[state_idx]) {
                    ++diffs;
                }
            }
            if (static_cast<double>(diffs) / static_cast<double>(comparable) <= max_diff) {
                assigned = labels[rep_idx];
                break;
            }
        }
        if (!assigned.has_value()) {
            std::ostringstream label;
            label << "A" << std::setw(3) << std::setfill('0') << (labels.size() + 1);
            assigned = label.str();
            representatives.push_back(states);
            labels.push_back(*assigned);
        }
        rows.push_back(HaplotypeDetailRow{sample_name, *assigned});
    }
    return rows;
}

std::vector<HaplotypeSummaryRow> summarize_accessions(const std::vector<HaplotypeDetailRow>& accessions) {
    std::map<std::string, int> counts;
    for (const auto& accession : accessions) {
        counts[accession.hap] += 1;
    }

    std::vector<HaplotypeSummaryRow> rows;
    rows.reserve(counts.size());
    for (const auto& entry : counts) {
        rows.push_back(HaplotypeSummaryRow{entry.first, entry.second});
    }
    std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
        if (left.count != right.count) {
            return left.count > right.count;
        }
        return left.hap < right.hap;
    });
    return rows;
}

}  // namespace

ViewResult build_view_result(const RegionData& data, const ViewOptions& options) {
    ViewResult result;
    result.grouping_mode =
        options.max_diff.has_value()
            ? (options.by == GroupBy::Region ? "approx-region" : "approx-site")
            : (options.by == GroupBy::Region ? "strict-region" : "strict-site");
    result.grouping_method = options.max_diff.has_value() ? "max-diff" : "exact";
    result.output_mode = options.output_mode == OutputMode::Both ? "both"
        : options.output_mode == OutputMode::Summary ? "summary" : "detail";
    result.imputed_ref = options.impute;
    result.max_diff = options.max_diff;
    result.variant_count = static_cast<int>(data.variants.size());
    result.sample_count = static_cast<int>(data.samples.size());
    result.sites.reserve(data.variants.size());
    for (const auto& variant : data.variants) {
        std::ostringstream allele;
        allele << variant.ref << "/";
        for (std::size_t idx = 0; idx < variant.alt.size(); ++idx) {
            if (idx > 0) {
                allele << ",";
            }
            allele << variant.alt[idx];
        }
        result.sites.push_back(SiteRow{variant.chrom, variant.pos, allele.str()});
    }

    std::vector<HaplotypeDetailRow> accessions;
    if (!data.variants.empty()) {
        if (options.max_diff.has_value()) {
            accessions = build_approx_accessions(data, options.impute, *options.max_diff);
        } else {
            accessions = build_sample_profiles(data, options.impute);
        }
    }

    const auto summary = summarize_accessions(accessions);
    result.haplotype_count = static_cast<int>(summary.size());

    // Always populate both fields for Both mode; otherwise per output_mode
    if (options.output_mode == OutputMode::Both) {
        result.haplotypes = summary;
        result.accessions = std::move(accessions);
    } else if (options.output_mode == OutputMode::Summary) {
        result.haplotypes = summary;
    } else {
        result.accessions = std::move(accessions);
    }
    return result;
}

std::string serialize_view_result_json(const ViewResult& result) {
    using json = nlohmann::json;
    json j;
    j["grouping_method"] = result.grouping_method;
    j["grouping_mode"] = result.grouping_mode;
    j["haplotype_count"] = result.haplotype_count;
    j["imputed_ref"] = result.imputed_ref;
    j["max_diff"] = result.max_diff.has_value() ? json(*result.max_diff) : json(nullptr);
    j["output_mode"] = result.output_mode;
    j["sample_count"] = result.sample_count;
    j["variant_count"] = result.variant_count;

    json sites = json::array();
    for (const auto& site : result.sites) {
        sites.push_back({{"allele", site.allele}, {"chrom", site.chrom}, {"pos", site.pos}});
    }
    j["sites"] = sites;

    if (result.output_mode == "summary" || result.output_mode == "both") {
        json haps = json::array();
        for (const auto& h : result.haplotypes) {
            haps.push_back({{"count", h.count}, {"hap", h.hap}});
        }
        j["haplotypes"] = haps;
    }

    if (result.output_mode == "detail" || result.output_mode == "both") {
        json accs = json::array();
        for (const auto& a : result.accessions) {
            accs.push_back({{"hap", a.hap}, {"sample", a.sample}});
        }
        j["accessions"] = accs;
    }

    if (result.annotation.has_value()) {
        const auto& ann = *result.annotation;
        json ann_json;
        ann_json["mode"] = ann.mode;
        if (ann.mode != "none") {
            ann_json["id"] = ann.id;
            ann_json["seqid"] = ann.seqid;
            ann_json["start"] = ann.start;
            ann_json["end"] = ann.end;
            ann_json["strand"] = std::string(1, ann.strand);
        }
        if (ann.mode == "nearest" && ann.distance.has_value()) {
            ann_json["distance"] = *ann.distance;
        }
        j["annotation"] = ann_json;
    }

    return j.dump();
}

std::vector<std::string> load_sample_list(const std::string& path) {
    std::ifstream handle(path);
    if (!handle) {
        throw std::runtime_error("failed to open sample list: " + path);
    }

    std::vector<std::string> samples;
    std::string line;
    while (std::getline(handle, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        samples.push_back(line);
    }
    return samples;
}

}  // namespace haptools
