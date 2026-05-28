#ifndef QC_VALIDATE_HPP
#define QC_VALIDATE_HPP

#include "annotation.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gffsub {

bool is_gff3_version(std::string_view version);
int count_tab_delimited_columns(const std::string& line);
bool is_fasta_boundary(const std::string& line);
std::vector<std::string> split_tab_fields(const std::string& line);
std::optional<std::string> raw_attr_value(std::string_view attrs, std::string_view key);
bool parse_qc_int64(std::string_view value, int64_t& out);
bool parse_qc_score(std::string_view value, std::optional<double>& out);
bool parse_positive_int64(std::string_view value, int64_t& out);
std::optional<std::string> attribute_syntax_error(std::string_view attrs);
std::optional<std::string> duplicate_attribute_tag(std::string_view attrs);
std::optional<std::string> empty_attribute_value_tag(std::string_view attrs);
std::optional<std::string> invalid_multi_value_attribute_tag(std::string_view attrs);
std::optional<std::string> percent_encoding_error(std::string_view attrs);
std::optional<std::string> attribute_escape_error(std::string_view attrs);
std::optional<std::string> target_attribute_error(std::string_view target);
std::optional<std::string> gap_attribute_error(std::string_view gap);
std::optional<std::string> database_accession_error(std::string_view label, std::string_view value);
bool has_circular_true(const std::unordered_map<std::string, std::vector<std::string>>& attrs);
bool has_invalid_is_circular(const std::unordered_map<std::string, std::vector<std::string>>& attrs);
bool allowed_discontinuous_id(const std::vector<const GffRecord*>& records);
std::unordered_set<std::string> find_parent_cycle_ids(
    const std::unordered_map<std::string, std::vector<std::string>>& parents_by_id);
std::optional<std::string> seqid_syntax_error(std::string_view seqid);
std::optional<std::string> feature_type_syntax_error(std::string_view type);

}  // namespace gffsub

#endif  // QC_VALIDATE_HPP
