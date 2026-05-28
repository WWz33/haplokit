#include "gtf_parser.hpp"

#include <optional>
#include <string>

namespace gffsub {
namespace {

std::optional<std::string> extract_quoted_value(const std::string& attrs, const char* key) {
    size_t pos = attrs.find(key);
    if (pos == std::string::npos) return std::nullopt;
    size_t q1 = attrs.find('"', pos);
    if (q1 == std::string::npos) return std::nullopt;
    size_t q2 = attrs.find('"', q1 + 1);
    if (q2 == std::string::npos) return std::nullopt;
    return attrs.substr(q1 + 1, q2 - q1 - 1);
}

}  // namespace

void apply_gtf_attributes(GffRecord& rec) {
    if (!rec.gene_id) {
        rec.gene_id = extract_quoted_value(rec.attr_raw, "gene_id");
    }
    if (!rec.transcript_id) {
        rec.transcript_id = extract_quoted_value(rec.attr_raw, "transcript_id");
    }
}

}  // namespace gffsub
