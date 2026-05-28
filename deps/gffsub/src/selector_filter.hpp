#ifndef SELECTOR_FILTER_HPP
#define SELECTOR_FILTER_HPP

#include "annotation.hpp"

#include <optional>
#include <regex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace gffsub {

struct GrepFilter {
    std::string field;
    std::string pattern;
    bool use_regex = false;
    bool ignore_case = false;
    std::optional<std::regex> compiled;
};

enum class ExprOp {
    Equal,
    NotEqual,
    Regex,
    NotRegex,
    Less,
    LessEqual,
    Greater,
    GreaterEqual
};

struct ExprFilter {
    std::string field;
    ExprOp op;
    std::string value;
    bool ignore_case = false;
    std::optional<std::regex> compiled;
};

struct ExprNode {
    enum class Kind {
        Predicate,
        Not,
        And,
        Or
    };

    Kind kind = Kind::Predicate;
    ExprFilter predicate;
    std::vector<ExprNode> children;
};

std::optional<std::pair<std::string, std::string>> parse_field_pattern(std::string_view value);
bool parse_expr_filters(std::string_view expr, std::vector<ExprNode>& out, std::string& error);
bool compile_filter_regexes(std::vector<GrepFilter>& grep_filters,
                            std::vector<ExprNode>& include_filters,
                            std::vector<ExprNode>& exclude_filters,
                            bool ignore_case,
                            std::string& error);
void filter_by_grep(GffData& data, const std::vector<GrepFilter>& filters, bool invert);
void filter_by_expr(GffData& data, const std::vector<ExprNode>& filters, bool include);

}  // namespace gffsub

#endif  // SELECTOR_FILTER_HPP
