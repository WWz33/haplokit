#ifndef QUERY_SUMMARY_HPP
#define QUERY_SUMMARY_HPP

#include "annotation.hpp"

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

namespace gffsub {

struct SummaryRow {
    std::string query_id;
    std::string matched_id;
    std::string matched_by;
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
    char strand = '.';
    std::string type;
    std::string parent_id;
    size_t child_count = 0;
    size_t transcript_count = 0;
    size_t exon_count = 0;
    int64_t cds_length = 0;
    std::string status;
    std::vector<std::string> attrs;
};

std::string record_id(const GffRecord& rec);
SummaryRow make_summary_row(const AnnotationIndex& index,
                            const std::string& query_id,
                            const std::string& matched_by,
                            const GffRecord& rec);
SummaryRow make_not_found_row(const std::string& query_id, const std::string& matched_by);
std::string infer_gene_match_key(const AnnotationIndex& index, const std::string& query, const GffRecord& rec);
std::vector<std::string> extract_output_attrs(const std::string& attrs, const std::vector<std::string>& keys);
void print_summary_tsv(std::ostream& out,
                       const std::vector<SummaryRow>& rows,
                       const std::vector<std::string>& output_attrs);
void print_summary_json(std::ostream& out,
                        const std::vector<SummaryRow>& rows,
                        const std::vector<std::string>& output_attrs);

}  // namespace gffsub

#endif  // QUERY_SUMMARY_HPP
