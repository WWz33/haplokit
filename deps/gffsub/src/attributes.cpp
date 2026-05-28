#include "gff3.hpp"

namespace gffsub {

std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs) {
    std::unordered_map<std::string, std::vector<std::string>> parsed;
    size_t pos = 0;
    while (pos < attrs.size()) {
        const size_t key_end = attrs.find('=', pos);
        if (key_end == std::string_view::npos) {
            break;
        }

        const std::string key{attrs.substr(pos, key_end - pos)};
        const size_t value_start = key_end + 1;
        size_t value_end = attrs.find(';', value_start);
        if (value_end == std::string_view::npos) {
            value_end = attrs.size();
        }

        const std::string value{attrs.substr(value_start, value_end - value_start)};
        if (!key.empty() && !value.empty()) {
            size_t part_start = 0;
            while (part_start <= value.size()) {
                size_t part_end = value.find(',', part_start);
                if (part_end == std::string::npos) {
                    part_end = value.size();
                }
                const std::string part = value.substr(part_start, part_end - part_start);
                if (!part.empty()) {
                    parsed[key].push_back(part);
                }
                if (part_end == value.size()) {
                    break;
                }
                part_start = part_end + 1;
            }
        }

        pos = (value_end < attrs.size()) ? value_end + 1 : attrs.size();
    }
    return parsed;
}

}  // namespace gffsub
