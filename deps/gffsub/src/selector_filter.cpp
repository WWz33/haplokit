#include "selector_filter.hpp"

#include "gff3.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <sstream>

namespace gffsub {

namespace {

bool is_supported_filter_field(std::string_view field) {
    if (field == "seqid" || field == "source" || field == "type" || field == "start" ||
        field == "end" || field == "length" || field == "score" || field == "strand" ||
        field == "phase" || field == "attrs" || field == "attributes") {
        return true;
    }
    if (field.rfind("attr.", 0) == 0 && field.size() > 5) {
        return true;
    }
    return field == "ID" || field == "id" || field == "Name" || field == "name" ||
           field == "Parent" || field == "parent" || field == "Alias" || field == "alias" ||
           field == "Dbxref" || field == "dbxref" || field == "Note" || field == "note" ||
           field == "biotype" ||
           field == "gene_id" || field == "transcript_id" || field == "locus_tag";
}

std::string join_filter_values(const std::vector<std::string>& values) {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << values[i];
    }
    return out.str();
}

std::optional<std::string> record_field_value(const GffRecord& rec, std::string_view field) {
    if (field == "seqid") return rec.seqid;
    if (field == "source") return rec.source;
    if (field == "type") return rec.type;
    if (field == "start") return std::to_string(rec.start);
    if (field == "end") return std::to_string(rec.end);
    if (field == "length") return std::to_string(rec.end - rec.start + 1);
    if (field == "score") return rec.score ? std::to_string(*rec.score) : ".";
    if (field == "strand") return std::string(1, rec.strand);
    if (field == "phase") return std::string(1, rec.phase);
    if (field == "attrs") return rec.attr_raw;
    if (field == "attributes") return rec.attr_raw;

    std::string attr_key;
    if (field.rfind("attr.", 0) == 0) {
        attr_key = std::string{field.substr(5)};
    } else if (field == "ID" || field == "id") {
        attr_key = "ID";
    } else if (field == "Name" || field == "name") {
        attr_key = "Name";
    } else if (field == "Parent" || field == "parent") {
        attr_key = "Parent";
    } else if (field == "Alias" || field == "alias") {
        attr_key = "Alias";
    } else if (field == "Dbxref" || field == "dbxref") {
        attr_key = "Dbxref";
    } else if (field == "Note" || field == "note") {
        attr_key = "Note";
    } else if (field == "biotype" ||
               field == "gene_id" || field == "transcript_id" || field == "locus_tag") {
        attr_key = std::string{field};
    } else {
        return std::nullopt;
    }

    const auto attrs = parse_attributes(rec.attr_raw);
    const auto it = attrs.find(attr_key);
    if (it == attrs.end()) {
        return ".";
    }
    return join_filter_values(it->second);
}

bool parse_number(std::string_view value, double& out) {
    if (value == "." || value.empty()) {
        return false;
    }
    try {
        size_t parsed = 0;
        out = std::stod(std::string{value}, &parsed);
        return parsed == value.size() && std::isfinite(out);
    } catch (const std::exception&) {
        return false;
    }
}

bool is_numeric_filter_field(std::string_view field) {
    return field == "start" || field == "end" || field == "length" || field == "score";
}

std::string lowercase_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool field_matches_grep(const GffRecord& rec, const GrepFilter& filter) {
    const auto value = record_field_value(rec, filter.field);
    if (!value) {
        return false;
    }
    if (filter.use_regex) {
        return filter.compiled && std::regex_search(*value, *filter.compiled);
    }
    if (filter.ignore_case) {
        return lowercase_copy(*value).find(lowercase_copy(filter.pattern)) != std::string::npos;
    }
    return value->find(filter.pattern) != std::string::npos;
}

bool field_matches_expr(const GffRecord& rec, const ExprFilter& filter) {
    const auto value = record_field_value(rec, filter.field);
    if (!value) {
        return false;
    }
    const auto comparable_value = filter.ignore_case ? lowercase_copy(*value) : *value;
    const auto comparable_filter = filter.ignore_case ? lowercase_copy(filter.value) : filter.value;
    switch (filter.op) {
        case ExprOp::Equal:
        case ExprOp::NotEqual: {
            if (is_numeric_filter_field(filter.field)) {
                double lhs = 0.0;
                double rhs = 0.0;
                const bool numeric_match = parse_number(*value, lhs) && parse_number(filter.value, rhs) && lhs == rhs;
                return filter.op == ExprOp::Equal ? numeric_match : !numeric_match;
            }
            return filter.op == ExprOp::Equal ? comparable_value == comparable_filter : comparable_value != comparable_filter;
        }
        case ExprOp::Regex: return filter.compiled && std::regex_search(*value, *filter.compiled);
        case ExprOp::NotRegex: return filter.compiled && !std::regex_search(*value, *filter.compiled);
        case ExprOp::Less:
        case ExprOp::LessEqual:
        case ExprOp::Greater:
        case ExprOp::GreaterEqual: {
            double lhs = 0.0;
            double rhs = 0.0;
            if (!parse_number(*value, lhs) || !parse_number(filter.value, rhs)) {
                return false;
            }
            if (filter.op == ExprOp::Less) return lhs < rhs;
            if (filter.op == ExprOp::LessEqual) return lhs <= rhs;
            if (filter.op == ExprOp::Greater) return lhs > rhs;
            return lhs >= rhs;
        }
    }
    return false;
}

bool expr_node_matches(const GffRecord& rec, const ExprNode& node) {
    switch (node.kind) {
        case ExprNode::Kind::Predicate:
            return field_matches_expr(rec, node.predicate);
        case ExprNode::Kind::Not:
            return node.children.empty() || !expr_node_matches(rec, node.children.front());
        case ExprNode::Kind::And:
            for (const auto& child : node.children) {
                if (!expr_node_matches(rec, child)) {
                    return false;
                }
            }
            return true;
        case ExprNode::Kind::Or:
            for (const auto& child : node.children) {
                if (expr_node_matches(rec, child)) {
                    return true;
                }
            }
            return false;
    }
    return false;
}

std::string trim_copy(std::string_view value) {
    size_t start = 0;
    while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start]))) {
        ++start;
    }

    size_t end = value.size();
    while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }

    return std::string{value.substr(start, end - start)};
}

std::string unquote_expr_value(std::string value) {
    value = trim_copy(value);
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
        std::string out;
        out.reserve(value.size() - 2);
        for (size_t i = 1; i + 1 < value.size(); ++i) {
            if (value[i] == '\\' && i + 2 < value.size()) {
                const char next = value[i + 1];
                if (next == '"' || next == '\\') {
                    out.push_back(next);
                    ++i;
                    continue;
                }
            }
            out.push_back(value[i]);
        }
        return out;
    }
    return value;
}

bool parse_expr_condition(std::string_view condition, ExprFilter& out, std::string& error) {
    const std::pair<std::string_view, ExprOp> ops[] = {
        {"!~", ExprOp::NotRegex},
        {"==", ExprOp::Equal},
        {"!=", ExprOp::NotEqual},
        {"<=", ExprOp::LessEqual},
        {">=", ExprOp::GreaterEqual},
        {"~", ExprOp::Regex},
        {"<", ExprOp::Less},
        {">", ExprOp::Greater}
    };

    bool in_quotes = false;
    bool escaped = false;
    for (size_t i = 0; i < condition.size(); ++i) {
        const char ch = condition[i];
        if (escaped) {
            escaped = false;
            continue;
        }
        if (ch == '\\') {
            escaped = true;
            continue;
        }
        if (ch == '"') {
            in_quotes = !in_quotes;
            continue;
        }
        if (in_quotes) {
            continue;
        }
        for (const auto& [token, op] : ops) {
            if (condition.substr(i, token.size()) == token) {
                out.field = trim_copy(condition.substr(0, i));
                out.op = op;
                out.value = unquote_expr_value(std::string{condition.substr(i + token.size())});
                if (out.field.empty() || out.value.empty()) {
                    error = "empty field or value in expression condition " + std::string{condition};
                    return false;
                }
                return true;
            }
        }
    }
    error = "expected operator ==, !=, ~, !~, <, <=, >, or >= in expression condition " + std::string{condition};
    return false;
}

class ExprParser {
public:
    explicit ExprParser(std::string_view expr) : expr_(expr) {}

    bool parse(ExprNode& out, std::string& error) {
        skip_space();
        if (!parse_or(out, error)) {
            return false;
        }
        skip_space();
        if (pos_ != expr_.size()) {
            error = "unexpected token near " + std::string{expr_.substr(pos_)};
            return false;
        }
        return true;
    }

private:
    std::string_view expr_;
    size_t pos_ = 0;

    void skip_space() {
        while (pos_ < expr_.size() && std::isspace(static_cast<unsigned char>(expr_[pos_]))) {
            ++pos_;
        }
    }

    bool consume(std::string_view token) {
        skip_space();
        if (expr_.substr(pos_, token.size()) == token) {
            pos_ += token.size();
            return true;
        }
        return false;
    }

    bool parse_or(ExprNode& out, std::string& error) {
        ExprNode first;
        if (!parse_and(first, error)) {
            return false;
        }
        std::vector<ExprNode> nodes;
        nodes.push_back(std::move(first));
        while (consume("||")) {
            ExprNode next;
            if (!parse_and(next, error)) {
                return false;
            }
            nodes.push_back(std::move(next));
        }
        if (nodes.size() == 1) {
            out = std::move(nodes.front());
        } else {
            out.kind = ExprNode::Kind::Or;
            out.children = std::move(nodes);
        }
        return true;
    }

    bool parse_and(ExprNode& out, std::string& error) {
        ExprNode first;
        if (!parse_unary(first, error)) {
            return false;
        }
        std::vector<ExprNode> nodes;
        nodes.push_back(std::move(first));
        while (consume("&&")) {
            ExprNode next;
            if (!parse_unary(next, error)) {
                return false;
            }
            nodes.push_back(std::move(next));
        }
        if (nodes.size() == 1) {
            out = std::move(nodes.front());
        } else {
            out.kind = ExprNode::Kind::And;
            out.children = std::move(nodes);
        }
        return true;
    }

    bool parse_unary(ExprNode& out, std::string& error) {
        skip_space();
        if (pos_ < expr_.size() && expr_[pos_] == '!' && expr_.substr(pos_, 2) != "!=" && expr_.substr(pos_, 2) != "!~") {
            ++pos_;
            ExprNode child;
            if (!parse_unary(child, error)) {
                return false;
            }
            out.kind = ExprNode::Kind::Not;
            out.children.push_back(std::move(child));
            return true;
        }
        return parse_primary(out, error);
    }

    bool parse_primary(ExprNode& out, std::string& error) {
        skip_space();
        if (pos_ >= expr_.size()) {
            error = "unexpected end of expression";
            return false;
        }
        if (expr_[pos_] == '(') {
            ++pos_;
            if (!parse_or(out, error)) {
                return false;
            }
            if (!consume(")")) {
                error = "missing closing )";
                return false;
            }
            return true;
        }

        const size_t start = pos_;
        bool in_quotes = false;
        bool escaped = false;
        while (pos_ < expr_.size()) {
            const char ch = expr_[pos_];
            if (escaped) {
                escaped = false;
                ++pos_;
                continue;
            }
            if (ch == '\\') {
                escaped = true;
                ++pos_;
                continue;
            }
            if (ch == '"') {
                in_quotes = !in_quotes;
                ++pos_;
                continue;
            }
            if (!in_quotes && ch == ')') {
                break;
            }
            if (!in_quotes && pos_ + 1 < expr_.size() &&
                ((expr_[pos_] == '&' && expr_[pos_ + 1] == '&') ||
                 (expr_[pos_] == '|' && expr_[pos_ + 1] == '|'))) {
                break;
            }
            ++pos_;
        }

        const auto condition = trim_copy(expr_.substr(start, pos_ - start));
        if (condition.empty()) {
            error = "empty condition in expression";
            return false;
        }
        ExprFilter filter;
        if (!parse_expr_condition(condition, filter, error)) {
            return false;
        }
        out.kind = ExprNode::Kind::Predicate;
        out.predicate = std::move(filter);
        return true;
    }
};

}  // namespace

std::optional<std::pair<std::string, std::string>> parse_field_pattern(std::string_view value) {
    const auto colon = value.find(':');
    if (colon == std::string_view::npos || colon == 0 || colon + 1 == value.size()) {
        return std::nullopt;
    }
    return std::pair<std::string, std::string>{std::string{value.substr(0, colon)}, std::string{value.substr(colon + 1)}};
}

void filter_by_grep(GffData& data, const std::vector<GrepFilter>& filters, bool invert) {
    for (auto& rec : data) {
        if (!rec.kept) {
            continue;
        }
        bool matched = filters.empty();
        for (const auto& filter : filters) {
            if (field_matches_grep(rec, filter)) {
                matched = true;
                break;
            }
        }
        if (invert ? matched : !matched) {
            rec.kept = false;
        }
    }
}

void filter_by_expr(GffData& data, const std::vector<ExprNode>& filters, bool include) {
    for (auto& rec : data) {
        if (!rec.kept) {
            continue;
        }
        for (const auto& filter : filters) {
            const bool matched = expr_node_matches(rec, filter);
            if ((include && !matched) || (!include && matched)) {
                rec.kept = false;
                break;
            }
        }
    }
}

bool parse_expr_filters(std::string_view expr, std::vector<ExprNode>& out, std::string& error) {
    ExprNode node;
    ExprParser parser{expr};
    if (!parser.parse(node, error)) {
        return false;
    }
    out.push_back(std::move(node));
    return true;
}

bool compile_filter_regexes(std::vector<GrepFilter>& grep_filters,
                            std::vector<ExprNode>& include_filters,
                            std::vector<ExprNode>& exclude_filters,
                            bool ignore_case,
                            std::string& error) {
    const auto flags = ignore_case ? (std::regex::ECMAScript | std::regex::icase) : std::regex::ECMAScript;
    for (auto& filter : grep_filters) {
        if (!is_supported_filter_field(filter.field)) {
            error = "unknown filter field " + filter.field;
            return false;
        }
        filter.ignore_case = ignore_case;
        if (!filter.use_regex) {
            continue;
        }
        try {
            filter.compiled.emplace(filter.pattern, flags);
        } catch (const std::regex_error&) {
            error = "invalid regex in --grep-regex for " + filter.field + ":" + filter.pattern;
            return false;
        }
    }

    auto compile_expr_predicate = [&](ExprFilter& filter) {
        if (!is_supported_filter_field(filter.field)) {
            error = "unknown expression field " + filter.field;
            return false;
        }
        filter.ignore_case = ignore_case;
        if (filter.op != ExprOp::Regex && filter.op != ExprOp::NotRegex) {
            return true;
        }
        try {
            filter.compiled.emplace(filter.value, flags);
        } catch (const std::regex_error&) {
            error = "invalid regex in expression for " + filter.field + "~" + filter.value;
            return false;
        }
        return true;
    };

    auto compile_expr_node = [&](auto&& self, ExprNode& node) -> bool {
        if (node.kind == ExprNode::Kind::Predicate) {
            return compile_expr_predicate(node.predicate);
        }
        for (auto& child : node.children) {
            if (!self(self, child)) {
                return false;
            }
        }
        return true;
    };

    for (auto& filter : include_filters) {
        if (!compile_expr_node(compile_expr_node, filter)) {
            return false;
        }
    }
    for (auto& filter : exclude_filters) {
        if (!compile_expr_node(compile_expr_node, filter)) {
            return false;
        }
    }
    return true;
}

}  // namespace gffsub
