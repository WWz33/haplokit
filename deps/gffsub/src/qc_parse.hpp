#ifndef QC_PARSE_HPP
#define QC_PARSE_HPP

#include "gff3.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace gffsub {

struct DirectiveIssue {
    int line_idx;
    std::string code;
    std::string message;
};

struct DirectiveParseResult {
    std::unordered_map<std::string, Region> sequence_regions;
    std::vector<DirectiveIssue> issues;
};

struct QcParseResult {
    GffData data;
    std::vector<DirectiveIssue> issues;
    bool opened = false;
};

DirectiveParseResult parse_directives(const std::string& path);
QcParseResult parse_qc_records(const std::string& path);

}  // namespace gffsub

#endif  // QC_PARSE_HPP
