#include "qc_validate.hpp"

#include "gff3.hpp"

#include <algorithm>
#include <cctype>
#include <exception>
#include <sstream>

namespace gffsub {
namespace {

bool is_digits(std::string_view value) {
    return !value.empty() &&
           std::all_of(value.begin(), value.end(), [](const unsigned char c) { return std::isdigit(c); });
}

bool is_hex_digit(char c) {
    return std::isxdigit(static_cast<unsigned char>(c));
}

bool allows_multiple_attribute_values(std::string_view tag) {
    return tag == "Parent" || tag == "Alias" || tag == "Note" || tag == "Dbxref" || tag == "Ontology_term";
}

bool is_gap_operation(char op) {
    return op == 'M' || op == 'I' || op == 'D' || op == 'F' || op == 'R';
}

std::vector<std::string> sorted_attribute_values(
    const std::unordered_map<std::string, std::vector<std::string>>& attrs,
    const char* key) {
    const auto it = attrs.find(key);
    if (it == attrs.end()) {
        return {};
    }
    auto values = it->second;
    values.erase(std::remove(values.begin(), values.end(), ""), values.end());
    std::sort(values.begin(), values.end());
    return values;
}

bool share_discontinuous_relation(
    const std::unordered_map<std::string, std::vector<std::string>>& first_attrs,
    const std::unordered_map<std::string, std::vector<std::string>>& next_attrs) {
    const auto first_parents = sorted_attribute_values(first_attrs, "Parent");
    if (!first_parents.empty() && first_parents == sorted_attribute_values(next_attrs, "Parent")) {
        return true;
    }

    const auto first_derives = sorted_attribute_values(first_attrs, "Derives_from");
    return !first_derives.empty() && first_derives == sorted_attribute_values(next_attrs, "Derives_from");
}

}  // namespace

bool is_gff3_version(std::string_view version) {
    return version == "3" || version.rfind("3.", 0) == 0;
}

int count_tab_delimited_columns(const std::string& line) {
    return static_cast<int>(std::count(line.begin(), line.end(), '\t')) + 1;
}

bool is_fasta_boundary(const std::string& line) {
    return line.rfind("##FASTA", 0) == 0 || line.rfind(">", 0) == 0;
}

std::vector<std::string> split_tab_fields(const std::string& line) {
    std::vector<std::string> cols;
    cols.reserve(9);
    size_t start = 0;
    while (true) {
        const auto pos = line.find('\t', start);
        if (pos == std::string::npos) {
            cols.emplace_back(line.substr(start));
            break;
        }
        cols.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return cols;
}

std::optional<std::string> raw_attr_value(std::string_view attrs, std::string_view key) {
    size_t pos = 0;
    while (pos < attrs.size()) {
        const auto key_end = attrs.find('=', pos);
        if (key_end == std::string_view::npos) {
            const auto next = attrs.find(';', pos);
            pos = (next == std::string_view::npos) ? attrs.size() : next + 1;
            continue;
        }

        const auto found_key = attrs.substr(pos, key_end - pos);
        const size_t value_start = key_end + 1;
        auto value_end = attrs.find(';', value_start);
        if (value_end == std::string_view::npos) {
            value_end = attrs.size();
        }

        if (found_key == key) {
            return std::string{attrs.substr(value_start, value_end - value_start)};
        }
        pos = (value_end < attrs.size()) ? value_end + 1 : attrs.size();
    }
    return std::nullopt;
}

bool parse_qc_int64(std::string_view value, int64_t& out) {
    if (!is_digits(value)) {
        return false;
    }
    try {
        const std::string copy{value};
        size_t parsed = 0;
        const auto parsed_value = std::stoll(copy, &parsed);
        if (parsed != copy.size()) {
            return false;
        }
        out = parsed_value;
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool parse_qc_score(std::string_view value, std::optional<double>& out) {
    if (value == ".") {
        out = std::nullopt;
        return true;
    }
    if (value.empty()) {
        return false;
    }
    try {
        const std::string copy{value};
        size_t parsed = 0;
        const double parsed_value = std::stod(copy, &parsed);
        if (parsed != copy.size()) {
            return false;
        }
        out = parsed_value;
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool parse_positive_int64(std::string_view value, int64_t& out) {
    if (!is_digits(value)) {
        return false;
    }
    try {
        const std::string copy{value};
        size_t parsed = 0;
        const auto parsed_value = std::stoll(copy, &parsed);
        if (parsed != copy.size() || parsed_value < 1) {
            return false;
        }
        out = parsed_value;
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

std::optional<std::string> attribute_syntax_error(std::string_view attrs) {
    if (attrs == ".") {
        return std::nullopt;
    }
    if (attrs.empty()) {
        return "attributes field must be . or semicolon-separated tag=value fields";
    }

    size_t pos = 0;
    while (pos <= attrs.size()) {
        size_t end = attrs.find(';', pos);
        if (end == std::string_view::npos) {
            end = attrs.size();
        }
        const auto field = attrs.substr(pos, end - pos);
        if (!field.empty()) {
            const auto equals = field.find('=');
            if (equals == std::string_view::npos || equals == 0 || field.find('=', equals + 1) != std::string_view::npos) {
                return "attributes must be semicolon-separated tag=value fields";
            }
        }
        if (end == attrs.size()) {
            break;
        }
        pos = end + 1;
    }
    return std::nullopt;
}

std::optional<std::string> duplicate_attribute_tag(std::string_view attrs) {
    if (attrs == ".") {
        return std::nullopt;
    }

    std::unordered_set<std::string> seen;
    size_t pos = 0;
    while (pos <= attrs.size()) {
        size_t end = attrs.find(';', pos);
        if (end == std::string_view::npos) {
            end = attrs.size();
        }
        const auto field = attrs.substr(pos, end - pos);
        if (!field.empty()) {
            const auto equals = field.find('=');
            const std::string tag{field.substr(0, equals)};
            if (!seen.insert(tag).second) {
                return tag;
            }
        }
        if (end == attrs.size()) {
            break;
        }
        pos = end + 1;
    }
    return std::nullopt;
}

std::optional<std::string> empty_attribute_value_tag(std::string_view attrs) {
    if (attrs == ".") {
        return std::nullopt;
    }

    size_t pos = 0;
    while (pos <= attrs.size()) {
        size_t end = attrs.find(';', pos);
        if (end == std::string_view::npos) {
            end = attrs.size();
        }
        const auto field = attrs.substr(pos, end - pos);
        if (!field.empty()) {
            const auto equals = field.find('=');
            if (equals != std::string_view::npos && equals + 1 == field.size()) {
                return std::string{field.substr(0, equals)};
            }
            if (equals != std::string_view::npos) {
                const auto tag = field.substr(0, equals);
                const auto value = field.substr(equals + 1);
                if (!value.empty() && allows_multiple_attribute_values(tag) &&
                    (value.front() == ',' || value.back() == ',' ||
                     value.find(",,") != std::string_view::npos)) {
                    return std::string{tag};
                }
            }
        }
        if (end == attrs.size()) {
            break;
        }
        pos = end + 1;
    }
    return std::nullopt;
}

std::optional<std::string> invalid_multi_value_attribute_tag(std::string_view attrs) {
    if (attrs == ".") {
        return std::nullopt;
    }

    size_t pos = 0;
    while (pos <= attrs.size()) {
        size_t end = attrs.find(';', pos);
        if (end == std::string_view::npos) {
            end = attrs.size();
        }
        const auto field = attrs.substr(pos, end - pos);
        if (!field.empty()) {
            const auto equals = field.find('=');
            const auto tag = field.substr(0, equals);
            const auto value = field.substr(equals + 1);
            if (value.find(',') != std::string_view::npos && !allows_multiple_attribute_values(tag)) {
                return std::string{tag};
            }
        }
        if (end == attrs.size()) {
            break;
        }
        pos = end + 1;
    }
    return std::nullopt;
}

std::optional<std::string> percent_encoding_error(std::string_view attrs) {
    for (size_t i = 0; i < attrs.size(); ++i) {
        if (attrs[i] != '%') {
            continue;
        }
        if (i + 2 >= attrs.size() || !is_hex_digit(attrs[i + 1]) || !is_hex_digit(attrs[i + 2])) {
            return "percent escapes must be % followed by two hex digits";
        }
    }
    return std::nullopt;
}

std::optional<std::string> attribute_escape_error(std::string_view attrs) {
    if (attrs.find('&') != std::string_view::npos) {
        return "ampersand in attributes must be percent-escaped as %26";
    }
    if (attrs.find('"') != std::string_view::npos) {
        return "double quote in attributes must be percent-escaped as %22";
    }
    return std::nullopt;
}

std::optional<std::string> target_attribute_error(std::string_view target) {
    std::istringstream fields{std::string{target}};
    std::vector<std::string> tokens;
    std::string token;
    while (fields >> token) {
        tokens.push_back(token);
    }

    if (tokens.size() != 3 && tokens.size() != 4) {
        return "Target must have target_id start end [strand]";
    }

    int64_t start = 0;
    int64_t end = 0;
    if (!parse_positive_int64(tokens[1], start) || !parse_positive_int64(tokens[2], end)) {
        return "Target start and end must be positive integers";
    }
    if (start > end) {
        return "Target start is greater than end";
    }
    if (tokens.size() == 4 && tokens[3] != "+" && tokens[3] != "-") {
        return "Target strand must be + or -";
    }
    return std::nullopt;
}

std::optional<std::string> gap_attribute_error(std::string_view gap) {
    std::istringstream fields{std::string{gap}};
    std::string token;
    bool saw_token = false;
    while (fields >> token) {
        saw_token = true;
        if (token.size() < 2 || !is_gap_operation(token[0]) || !is_digits(std::string_view{token}.substr(1))) {
            return "Gap must contain operation-length pairs such as M8 D3 M6";
        }
        int64_t length = 0;
        if (!parse_positive_int64(std::string_view{token}.substr(1), length)) {
            return "Gap operation lengths must be positive integers";
        }
    }
    if (!saw_token) {
        return "Gap must contain at least one operation-length pair";
    }
    return std::nullopt;
}

std::optional<std::string> database_accession_error(std::string_view label, std::string_view value) {
    const auto colon = value.find(':');
    if (colon == std::string_view::npos || colon == 0 || colon + 1 == value.size()) {
        return std::string{label} + " must have database:accession";
    }
    return std::nullopt;
}

bool has_circular_true(const std::unordered_map<std::string, std::vector<std::string>>& attrs) {
    const auto circular_it = attrs.find("Is_circular");
    return circular_it != attrs.end() &&
           std::find(circular_it->second.begin(), circular_it->second.end(), "true") != circular_it->second.end();
}

bool has_invalid_is_circular(const std::unordered_map<std::string, std::vector<std::string>>& attrs) {
    const auto circular_it = attrs.find("Is_circular");
    return circular_it != attrs.end() &&
           std::find_if(circular_it->second.begin(), circular_it->second.end(),
                        [](const std::string& value) { return value != "true"; }) != circular_it->second.end();
}

bool allowed_discontinuous_id(const std::vector<const GffRecord*>& records) {
    if (records.size() < 2) {
        return true;
    }

    const auto* first = records.front();
    const auto first_attrs = parse_attributes(first->attr_raw);
    for (size_t i = 1; i < records.size(); ++i) {
        const auto* next = records[i];
        if (next->seqid != first->seqid || next->source != first->source || next->type != first->type ||
            next->strand != first->strand) {
            return false;
        }

        const auto next_attrs = parse_attributes(next->attr_raw);
        if (!share_discontinuous_relation(first_attrs, next_attrs)) {
            return false;
        }
    }
    return true;
}

std::unordered_set<std::string> find_parent_cycle_ids(
    const std::unordered_map<std::string, std::vector<std::string>>& parents_by_id) {
    std::unordered_map<std::string, int> state;
    std::unordered_map<std::string, size_t> stack_pos;
    std::unordered_set<std::string> cycle_ids;
    std::vector<std::string> stack;

    struct Frame {
        std::string id;
        size_t next_parent = 0;
    };

    for (const auto& [id, parents] : parents_by_id) {
        (void)parents;
        if (state.find(id) != state.end()) {
            continue;
        }

        std::vector<Frame> frames{{id, 0}};
        state[id] = 1;
        stack_pos[id] = stack.size();
        stack.push_back(id);

        while (!frames.empty()) {
            auto& frame = frames.back();
            const auto parents_it = parents_by_id.find(frame.id);
            if (parents_it == parents_by_id.end() || frame.next_parent == parents_it->second.size()) {
                stack.pop_back();
                stack_pos.erase(frame.id);
                state[frame.id] = 2;
                frames.pop_back();
                continue;
            }

            const auto& parent_id = parents_it->second[frame.next_parent++];
            if (parents_by_id.find(parent_id) == parents_by_id.end()) {
                continue;
            }

            const auto parent_state_it = state.find(parent_id);
            if (parent_state_it == state.end()) {
                state[parent_id] = 1;
                stack_pos[parent_id] = stack.size();
                stack.push_back(parent_id);
                frames.push_back({parent_id, 0});
                continue;
            }
            if (parent_state_it->second == 1) {
                const auto cycle_start_it = stack_pos.find(parent_id);
                if (cycle_start_it != stack_pos.end()) {
                    for (size_t i = cycle_start_it->second; i < stack.size(); ++i) {
                        cycle_ids.insert(stack[i]);
                    }
                }
            }
        }
    }

    return cycle_ids;
}

std::optional<std::string> seqid_syntax_error(std::string_view seqid) {
    if (seqid.empty()) {
        return "seqid must not be empty";
    }
    if (seqid.front() == '>') {
        return "seqid must not begin with unescaped >";
    }
    for (size_t i = 0; i < seqid.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(seqid[i]);
        if (std::isspace(c)) {
            return "seqid must not contain unescaped whitespace";
        }
        if (seqid[i] == '%') {
            if (i + 2 >= seqid.size() || !is_hex_digit(seqid[i + 1]) || !is_hex_digit(seqid[i + 2])) {
                return "seqid percent escapes must be % followed by two hex digits";
            }
            i += 2;
            continue;
        }
        if (std::isalnum(c) || seqid[i] == '.' || seqid[i] == ':' || seqid[i] == '^' ||
            seqid[i] == '*' || seqid[i] == '$' || seqid[i] == '@' || seqid[i] == '!' ||
            seqid[i] == '+' || seqid[i] == '_' || seqid[i] == '?' || seqid[i] == '-' ||
            seqid[i] == '|') {
            continue;
        }
        return std::string{"seqid contains unescaped character "} + seqid[i];
    }
    return std::nullopt;
}

std::optional<std::string> feature_type_syntax_error(std::string_view type) {
    if (type.empty() || type == ".") {
        return "feature type must be a Sequence Ontology term or accession";
    }
    if (type.rfind("SO:", 0) == 0) {
        const auto accession = type.substr(3);
        if (accession.size() != 7 ||
            !std::all_of(accession.begin(), accession.end(), [](const unsigned char c) { return std::isdigit(c); })) {
            return "Sequence Ontology accession must be SO: followed by seven digits";
        }
    }
    return std::nullopt;
}

}  // namespace gffsub
