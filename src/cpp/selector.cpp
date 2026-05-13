#include "selector.h"

#include <stdexcept>

namespace haptools {

Region parse_region(const std::string& value) {
    const auto colon = value.find(':');
    if (colon == std::string::npos) {
        throw std::invalid_argument("region must contain ':'");
    }

    Region region;
    region.chrom = value.substr(0, colon);
    const auto coords = value.substr(colon + 1);
    const auto dash = coords.find('-');
    if (dash == std::string::npos) {
        region.start = std::stoi(coords);
        region.end = region.start;
        return region;
    }

    region.start = std::stoi(coords.substr(0, dash));
    region.end = std::stoi(coords.substr(dash + 1));
    return region;
}

bool is_site_region(const Region& region) {
    return region.start == region.end;
}

}  // namespace haptools
