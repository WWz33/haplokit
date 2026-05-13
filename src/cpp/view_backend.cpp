#include "view_backend.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
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

std::vector<HaplotypeDetailRow> build_exact_accessions(const RegionData& data, bool impute_ref) {
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

int state_distance(const std::string& left, const std::string& right) {
    return left == right ? 0 : 1;
}

std::vector<HaplotypeDetailRow> build_approx_accessions(
    const RegionData& data,
    bool impute_ref,
    double max_diff) {
    std::vector<std::pair<std::string, std::vector<std::string>>> raw_profiles;
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
        raw_profiles.emplace_back(data.samples[sample_idx], std::move(states));
    }

    std::vector<HaplotypeDetailRow> rows;
    std::vector<std::vector<std::string>> representatives;
    std::vector<std::string> labels;
    for (const auto& profile : raw_profiles) {
        const auto& sample_name = profile.first;
        const auto& states = profile.second;
        std::optional<std::string> assigned;
        for (std::size_t rep_idx = 0; rep_idx < representatives.size(); ++rep_idx) {
            const auto& representative = representatives[rep_idx];
            const auto comparable = std::max(representative.size(), states.size());
            if (comparable == 0) {
                continue;
            }
            int diffs = 0;
            for (std::size_t state_idx = 0; state_idx < representative.size() && state_idx < states.size(); ++state_idx) {
                diffs += state_distance(representative[state_idx], states[state_idx]);
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

std::string escape_json(const std::string& value) {
    std::ostringstream escaped;
    for (const char ch : value) {
        switch (ch) {
            case '\\':
                escaped << "\\\\";
                break;
            case '"':
                escaped << "\\\"";
                break;
            case '\n':
                escaped << "\\n";
                break;
            case '\r':
                escaped << "\\r";
                break;
            case '\t':
                escaped << "\\t";
                break;
            default:
                escaped << ch;
                break;
        }
    }
    return escaped.str();
}

std::string json_string(const std::string& value) {
    return "\"" + escape_json(value) + "\"";
}

std::string json_bool(bool value) {
    return value ? "true" : "false";
}

std::string json_number(double value) {
    std::ostringstream buffer;
    buffer << std::setprecision(15) << value;
    return buffer.str();
}

}  // namespace

ViewResult build_view_result(const RegionData& data, const ViewOptions& options) {
    ViewResult result;
    result.grouping_mode =
        options.max_diff.has_value()
            ? (options.by == GroupBy::Region ? "approx-region" : "approx-site")
            : (options.by == GroupBy::Region ? "strict-region" : "strict-site");
    result.grouping_method = options.max_diff.has_value() ? "max-diff" : "exact";
    result.output_mode = options.output_mode == OutputMode::Summary ? "summary" : "detail";
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
            accessions = build_exact_accessions(data, options.impute);
        }
    }

    const auto summary = summarize_accessions(accessions);
    result.haplotype_count = static_cast<int>(summary.size());
    if (options.output_mode == OutputMode::Summary) {
        result.haplotypes = summary;
    } else {
        result.accessions = std::move(accessions);
    }
    return result;
}

std::string serialize_view_result_json(const ViewResult& result) {
    std::ostringstream json;
    json << "{";
    json << "\"grouping_method\":" << json_string(result.grouping_method) << ",";
    json << "\"grouping_mode\":" << json_string(result.grouping_mode) << ",";
    json << "\"haplotype_count\":" << result.haplotype_count << ",";
    json << "\"imputed_ref\":" << json_bool(result.imputed_ref) << ",";
    json << "\"max_diff\":";
    if (result.max_diff.has_value()) {
        json << json_number(*result.max_diff);
    } else {
        json << "null";
    }
    json << ",";
    json << "\"output_mode\":" << json_string(result.output_mode) << ",";
    json << "\"sample_count\":" << result.sample_count << ",";
    json << "\"variant_count\":" << result.variant_count << ",";
    json << "\"sites\":[";
    for (std::size_t idx = 0; idx < result.sites.size(); ++idx) {
        if (idx > 0) {
            json << ",";
        }
        json << "{"
             << "\"allele\":" << json_string(result.sites[idx].allele) << ","
             << "\"chrom\":" << json_string(result.sites[idx].chrom) << ","
             << "\"pos\":" << result.sites[idx].pos
             << "}";
    }
    json << "]";

    if (result.output_mode == "summary") {
        json << ",\"haplotypes\":[";
        for (std::size_t idx = 0; idx < result.haplotypes.size(); ++idx) {
            if (idx > 0) {
                json << ",";
            }
            json << "{"
                 << "\"count\":" << result.haplotypes[idx].count << ","
                 << "\"hap\":" << json_string(result.haplotypes[idx].hap)
                 << "}";
        }
        json << "]";
    } else {
        json << ",\"accessions\":[";
        for (std::size_t idx = 0; idx < result.accessions.size(); ++idx) {
            if (idx > 0) {
                json << ",";
            }
            json << "{"
                 << "\"hap\":" << json_string(result.accessions[idx].hap) << ","
                 << "\"sample\":" << json_string(result.accessions[idx].sample)
                 << "}";
        }
        json << "]";
    }

    if (result.annotation.has_value()) {
        const auto& ann = *result.annotation;
        json << ",\"annotation\":{";
        json << "\"mode\":" << json_string(ann.mode);
        if (ann.mode != "none") {
            json << ",\"id\":" << json_string(ann.id);
            json << ",\"seqid\":" << json_string(ann.seqid);
            json << ",\"start\":" << ann.start;
            json << ",\"end\":" << ann.end;
            json << ",\"strand\":" << json_string(std::string(1, ann.strand));
        }
        if (ann.mode == "nearest" && ann.distance.has_value()) {
            json << ",\"distance\":" << *ann.distance;
        }
        json << "}";
    }

    json << "}";
    return json.str();
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
