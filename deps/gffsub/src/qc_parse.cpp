#include "qc_parse.hpp"

#include "qc_validate.hpp"

#include <fstream>
#include <optional>
#include <sstream>

namespace gffsub {

DirectiveParseResult parse_directives(const std::string& path) {
    DirectiveParseResult result;
    std::ifstream in{path};
    std::string line;
    int line_num = 0;
    int gff_version_count = 0;
    while (std::getline(in, line)) {
        ++line_num;
        if (is_fasta_boundary(line)) {
            break;
        }
        if (line.rfind("##gff-version", 0) == 0) {
            ++gff_version_count;
            std::istringstream fields{line};
            std::string directive;
            std::string version;
            std::string extra;
            fields >> directive >> version;
            if (directive != "##gff-version" || !is_gff3_version(version) || (fields >> extra)) {
                result.issues.push_back({line_num, "invalid_gff_version",
                                         "##gff-version must declare exactly one version beginning with 3"});
            }
            if (line_num != 1) {
                result.issues.push_back({line_num, "invalid_gff_version",
                                         "##gff-version must be the topmost line"});
            }
            continue;
        }
        if (line.rfind("##sequence-region ", 0) != 0) {
            if (!line.empty() && line[0] != '#') {
                if (count_tab_delimited_columns(line) != 9) {
                    result.issues.push_back({line_num, "invalid_column_count",
                                             "feature lines must contain exactly 9 tab-delimited columns"});
                    continue;
                }
                const auto last_tab = line.rfind('\t');
                const auto attrs = std::string_view{line}.substr(last_tab + 1);
                if (const auto error = attribute_syntax_error(attrs)) {
                    result.issues.push_back({line_num, "invalid_attribute_syntax", *error});
                } else if (const auto empty_value_tag = empty_attribute_value_tag(attrs)) {
                    result.issues.push_back({line_num, "invalid_attribute_value",
                                             "attribute tag " + *empty_value_tag + " must have a non-empty value"});
                } else if (const auto duplicate_tag = duplicate_attribute_tag(attrs)) {
                    result.issues.push_back({line_num, "duplicate_attribute_tag",
                                             "attribute tag " + *duplicate_tag + " appears more than once"});
                } else if (const auto tag = invalid_multi_value_attribute_tag(attrs)) {
                    result.issues.push_back({line_num, "invalid_attribute_multivalue",
                                             "attribute tag " + *tag + " must not contain comma-separated values"});
                }
            }
            continue;
        }

        std::istringstream fields{line};
        std::string directive;
        std::string seqid;
        std::string start_text;
        std::string end_text;
        std::string extra;
        if (!(fields >> directive >> seqid >> start_text >> end_text) || (fields >> extra)) {
            result.issues.push_back({line_num, "invalid_sequence_region",
                                     "malformed ##sequence-region directive"});
            continue;
        }
        if (const auto error = seqid_syntax_error(seqid)) {
            result.issues.push_back({line_num, "invalid_sequence_region",
                                     "invalid ##sequence-region seqid " + seqid + ": " + *error});
            continue;
        }
        int64_t start = 0;
        int64_t end = 0;
        if (!parse_positive_int64(start_text, start) || !parse_positive_int64(end_text, end)) {
            result.issues.push_back({line_num, "invalid_sequence_region",
                                     "invalid ##sequence-region coordinates for " + seqid});
            continue;
        }
        if (start < 1 || end < 1 || start > end) {
            result.issues.push_back({line_num, "invalid_sequence_region",
                                     "invalid ##sequence-region coordinates for " + seqid});
            continue;
        }
        if (result.sequence_regions.find(seqid) != result.sequence_regions.end()) {
            result.issues.push_back({line_num, "duplicate_sequence_region",
                                     "##sequence-region appears more than once for " + seqid});
            continue;
        }
        result.sequence_regions[seqid] = Region{seqid, start, end};
    }
    if (gff_version_count == 0) {
        result.issues.push_back({-1, "invalid_gff_version", "missing ##gff-version directive"});
    } else if (gff_version_count > 1) {
        result.issues.push_back({-1, "invalid_gff_version", "##gff-version appears more than once"});
    }
    return result;
}

QcParseResult parse_qc_records(const std::string& path) {
    QcParseResult result;
    std::ifstream in{path};
    if (!in.is_open()) {
        return result;
    }
    result.opened = true;
    std::string line;
    int line_num = 0;
    bool in_fasta = false;
    while (std::getline(in, line)) {
        ++line_num;
        if (in_fasta) {
            continue;
        }
        if (is_fasta_boundary(line)) {
            in_fasta = true;
            continue;
        }
        if (line.empty() || line[0] == '#') {
            continue;
        }
        const auto cols = split_tab_fields(line);
        if (cols.size() < 9) {
            continue;
        }
        GffRecord rec;
        rec.seqid = cols[0];
        rec.source = cols[1];
        rec.type = cols[2];
        rec.line_idx = static_cast<int>(result.data.size());
        rec.kept = true;

        if (cols[1].empty()) {
            result.issues.push_back({line_num, "invalid_source",
                                     "source field must not be empty; use . when source is unknown"});
        }

        if (!parse_qc_int64(cols[3], rec.start) || !parse_qc_int64(cols[4], rec.end)) {
            result.issues.push_back({line_num, "invalid_coordinate",
                                     "start and end must be integer 1-based coordinates"});
            continue;
        }

        if (!parse_qc_score(cols[5], rec.score)) {
            result.issues.push_back({line_num, "invalid_score",
                                     "score must be a finite floating point number or ."});
            rec.score = std::nullopt;
        }

        if (cols[6].empty()) {
            result.issues.push_back({line_num, "invalid_strand",
                                     "strand field must be +, -, ., or ?"});
        }
        if (cols[7].empty()) {
            result.issues.push_back({line_num, "invalid_phase",
                                     "phase field must be ., 0, 1, or 2"});
        }
        rec.strand = cols[6].empty() ? '.' : cols[6][0];
        rec.phase = cols[7].empty() ? '.' : cols[7][0];
        rec.attr_raw = cols[8];
        rec.id = raw_attr_value(cols[8], "ID");
        rec.parent_id = raw_attr_value(cols[8], "Parent");
        rec.gene_id = raw_attr_value(cols[8], "gene_id");
        rec.transcript_id = raw_attr_value(cols[8], "transcript_id");

        result.data.append(rec);
    }
    return result;
}

}  // namespace gffsub
