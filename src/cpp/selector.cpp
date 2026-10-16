#include "selector.h"

#include <stdexcept>
#include <string>

namespace haplokit {
namespace {

int parse_coordinate(const std::string& text, const std::string& value) {
    std::size_t consumed = 0;
    int parsed = 0;
    try {
        parsed = std::stoi(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("region coordinate is not a number in '" + value + "'");
    }
    if (consumed != text.size()) {
        throw std::invalid_argument("region coordinate is not a number in '" + value + "'");
    }
    return parsed;
}

}  // namespace

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
        region.start = parse_coordinate(coords, value);
        region.end = region.start;
    } else {
        region.start = parse_coordinate(coords.substr(0, dash), value);
        region.end = parse_coordinate(coords.substr(dash + 1), value);
    }

    if (region.chrom.empty()) {
        throw std::invalid_argument("region must name a sequence: '" + value + "'");
    }
    if (region.start < 1) {
        throw std::invalid_argument("region start must be >= 1 in '" + value + "'");
    }
    if (region.start > region.end) {
        throw std::invalid_argument("region start must be <= end in '" + value + "'");
    }
    return region;
}

bool is_site_region(const Region& region) {
    return region.start == region.end;
}

}  // namespace haplokit
