#pragma once

#include <optional>
#include <string>

namespace haplokit {

struct Region {
    std::string chrom;
    int start = 0;
    int end = 0;
};

Region parse_region(const std::string& value);
bool is_site_region(const Region& region);

}  // namespace haplokit
