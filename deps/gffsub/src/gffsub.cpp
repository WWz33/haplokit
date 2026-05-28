#include "gff3.hpp"
#include "qc_parse.hpp"
#include "qc_validate.hpp"
#include "query_summary.hpp"
#include "selector_filter.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <getopt.h>

using namespace gffsub;

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input.gff3> [options]\n"
        << "\n"
        << "Selector Options:\n"
        << "  --id ID\n"
        << "      Extract a feature by exact GFF3 ID. May be repeated.\n"
        << "      Top-level selector for exact-ID extraction.\n"
        << "  --ids FILE, --id-list FILE\n"
        << "      Extract exact feature IDs listed one per non-empty line.\n"
        << "      Top-level selector for batch exact-ID extraction.\n"
        << "  --name NAME\n"
        << "      Extract a gene by ID/Name/gene_id/locus_tag/Alias/Dbxref.\n"
        << "      Top-level selector for common gene naming keys.\n"
        << "  --where KEY=VALUE, --attr KEY=VALUE\n"
        << "      Extract features by an exact GFF3 attribute value. May be repeated.\n"
        << "      Top-level selector for exact column-9 KEY=VALUE matches. --attr is a compatibility alias.\n"
        << "  --grep FIELD:PATTERN\n"
        << "      Keep records whose GFF field or attribute contains PATTERN. May be repeated.\n"
        << "  --grep-regex FIELD:REGEX\n"
        << "      Keep records whose GFF field or attribute matches REGEX. May be repeated.\n"
        << "  --grep-file FILE --grep-field FIELD\n"
        << "      Keep records matching one non-empty pattern per file line against FIELD.\n"
        << "  -I, --include-expr EXPR\n"
        << "      Keep records matching a GFF expression, such as type==\"gene\" && attr.biotype==\"protein_coding\".\n"
        << "  -E, --exclude-expr EXPR\n"
        << "      Drop records matching a GFF expression, such as attr.Note~\"transposon|retroelement\".\n"
        << "  -v, --invert-match\n"
        << "      Invert --grep/--grep-regex/--grep-file matches.\n"
        << "  --ignore-case\n"
        << "      Match grep and expression strings case-insensitively.\n"
        << "  -C, --children, --include-children\n"
        << "      Include descendants of records matched by --id, --ids, --name, --where, or --nearest.\n"
        << "  --parents, --include-parents\n"
        << "      Include ancestors of records matched by --id, --ids, --name, --where, or --nearest.\n"
        << "  --model, --gene-model\n"
        << "      Include the full gene model for records matched by --id, --ids, --name, --where, or --nearest.\n"
        << "  --nearest REGION, --nearest-gene REGION\n"
        << "      Keep the nearest gene to a 1-based inclusive CHR:START-END region.\n"
        << "  --out-attrs KEYS, --output-attrs KEYS\n"
        << "      Output selected attributes as TSV/JSON fields. Implies query summary output.\n"
        << "      --attrs remains as a deprecated alias.\n"
        << "  --summary FMT, --summary-format FMT\n"
        << "      Output summary instead of GFF3. Choices: tsv, json.\n"
        << "  --up N, --down N, --upstream N, --downstream N, --strand-aware\n"
        << "      Extract records overlapping the expanded window around --id.\n"
        << "  --qc\n"
        << "      Run annotation QC.\n"
        << "\n"
        << "Input/Region Options:\n"
        << "  --seqid SEQID\n"
        << "      Extract features whose first column exactly matches SEQID.\n"
        << "\n"
        << "  --source SOURCE\n"
        << "      Extract features whose second column exactly matches SOURCE.\n"
        << "\n"
        << "  --score SCORE\n"
        << "      Filter features by score (6th column). Use . for missing score.\n"
        << "\n"
        << "  -r, --region CHR:START-END\n"
        << "      Extract features overlapping the specified genomic region.\n"
        << "      Coordinates are 1-based and inclusive (GFF format).\n"
        << "      Example: -r chr1:1000000-2000000\n"
        << "\n"
        << "  -b, --bed FILE\n"
        << "      Extract features using genomic regions from a BED file.\n"
        << "      BED files use 0-based half-open coordinates, automatically\n"
        << "      converted to 1-based for internal processing.\n"
        << "\n"
        << "Feature Filter Options:\n"
        << "  -f, --feature TYPE, --type TYPE\n"
        << "      Filter features by type (3rd column in GFF/GTF).\n"
        << "      Examples: gene, mRNA, exon, CDS, transcript\n"
        << "\n"
        << "  --strand STRAND\n"
        << "      Filter features by strand (7th column). Choices: +, -, ., ?\n"
        << "\n"
        << "  --phase PHASE\n"
        << "      Filter features by phase (8th column). Choices: 0, 1, 2, .\n"
        << "\n"
        << "  -L, --longest\n"
        << "      Keep only the longest transcript isoform for each gene.\n"
        << "      Per-gene comparison (AGAT logic): if gene has CDS isoforms,\n"
        << "      only compare by CDS length; otherwise compare by exon length.\n"
        << "\n"
        << "  -@, --threads N\n"
        << "      Number of threads for parallel processing (default: 1).\n"
        << "      Currently used with --longest for multi-chromosome parallelization.\n"
        << "\n"
        << "Output Options:\n"
        << "  -t, --format FMT, --output-format FMT\n"
        << "      Output format. Choices: gff3, gtf2, gtf3, bed\n"
        << "      --output-format remains as a verbose alias.\n"
        << "      gff3  - GFF3 format (##gff-version 3)\n"
        << "      gtf2  - GTF2 format (##gtf-version 2)\n"
        << "      gtf3  - GTF3/Ensembl format (##gtf-version 2.2.1)\n"
        << "      bed   - BED format (0-based half-open coordinates)\n"
        << "      Default: gff3\n"
        << "\n"
        << "  -o, --output FILE\n"
        << "      Output file path. If not specified, writes to stdout.\n"
        << "\n"
        << "  -h, --help\n"
        << "      Display this help message.\n"
        << "\n"
        << "Attribute column example:\n"
        << "  ID=gene0001;Name=ABC1;Alias=ABC-1;Dbxref=GeneID:123\n"
        << "  --id searches ID exactly; --name searches common gene naming keys.\n"
        << "  --where KEY=VALUE filters records; --out-attrs prints selected attributes.\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " annotation.gff3 --id GeneA\n"
        << "  " << prog << " annotation.gff3 --ids genes.txt\n"
        << "  " << prog << " annotation.gff3 --id GeneA -C\n"
        << "  " << prog << " annotation.gff3 --id GeneA --summary tsv\n"
        << "  " << prog << " annotation.gff3 --id GeneA --out-attrs ID,Name,Parent\n"
        << "  " << prog << " annotation.gff3 --id GeneA --up 2000 --down 500\n"
        << "  " << prog << " annotation.gff3 --qc\n"
        << "  " << prog << " annotation.gff3 --name ABC1\n"
        << "  " << prog << " annotation.gff3 --where biotype=protein_coding\n"
        << "  " << prog << " annotation.gff3 --grep-file genes.txt --grep-field ID\n"
        << "  " << prog << " annotation.gff3 --grep-regex 'ID:^Glyma\\\\.01G'\n"
        << "  " << prog << " annotation.gff3 -I 'type==\"gene\" && attr.biotype==\"protein_coding\"'\n"
        << "  " << prog << " annotation.gff3 -E 'attr.Note~\"transposon|retroelement\"'\n"
        << "  " << prog << " annotation.gff3 --seqid chr1 -f gene\n"
        << "  " << prog << " annotation.gff3 --source Gnomon -f mRNA\n"
        << "  " << prog << " annotation.gff3 --score 42.5\n"
        << "  " << prog << " annotation.gff3 --strand - -f gene\n"
        << "  " << prog << " annotation.gff3 --phase 0 -f CDS\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 -f gene\n"
        << "  " << prog << " annotation.gff3 --bed regions.bed -f exon\n"
        << "  " << prog << " annotation.gff3 --longest\n"
        << "  " << prog << " annotation.gff3 --longest -@ 6\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 --format gtf3 -o out.gtf\n";
}

static void query_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " query <input.gff3> [options]\n"
        << "\n"
        << "Most workflows can use the top-level form: " << prog << " <input.gff3> [selectors]\n"
        << "\n"
        << "Query Options:\n"
        << "  --id ID                 Query a feature by ID.\n"
        << "  --name NAME             Query a gene by ID/Name/gene_id/locus_tag/Alias/Dbxref.\n"
        << "  --ids FILE              Query one feature ID per line.\n"
        << "  --id-list FILE          Verbose alias for --ids.\n"
        << "  --region CHR:START-END  Query features overlapping a 1-based inclusive region.\n"
        << "  --type TYPE             Restrict query output by feature type.\n"
        << "  --where KEY=VALUE       Query features by an exact GFF3 attribute value.\n"
        << "  --attr KEY=VALUE        Compatibility alias for --where.\n"
        << "  --out-attrs KEYS        Output selected attributes as extra TSV/JSON fields.\n"
        << "  --output-attrs KEYS     Verbose alias for --out-attrs.\n"
        << "  --attrs KEYS            Deprecated alias for --out-attrs.\n"
        << "  -C, --children           Include descendants of matched IDs.\n"
        << "  --include-children       Verbose alias for --children.\n"
        << "  --parents                Include ancestors of matched IDs.\n"
        << "  --include-parents        Verbose alias for --parents.\n"
        << "  --model                  Include the full gene model for matched records.\n"
        << "  --gene-model             Verbose alias for --model.\n"
        << "  --nearest REGION         Query the nearest gene to a 1-based inclusive CHR:START-END region.\n"
        << "  --nearest-gene REGION    Verbose alias for --nearest.\n"
        << "  --summary FMT           Output query summary instead of GFF3. Choices: tsv, json.\n"
        << "  --summary-format FMT    Verbose alias for --summary.\n"
        << "  -h, --help              Display this help message.\n"
        << "\n"
        << "Attribute column example:\n"
        << "  ID=gene0001;Name=ABC1;Alias=ABC-1;Dbxref=GeneID:123\n"
        << "  --where Alias=ABC-1 filters records; --out-attrs ID,Name prints fields.\n";
}

static void window_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " window <input.gff3> [options]\n"
        << "\n"
        << "Top-level equivalent: " << prog << " <input.gff3> --id ID --up N --down N\n"
        << "\n"
        << "Window Options:\n"
        << "  --id ID                 Target feature ID or gene lookup key.\n"
        << "  --up N, --upstream N    Bases to add upstream of the target (default: 0).\n"
        << "  --down N, --downstream N\n"
        << "                          Bases to add downstream of the target (default: 0).\n"
        << "  --strand-aware          Interpret upstream/downstream by feature strand.\n"
        << "  -h, --help              Display this help message.\n";
}

static void qc_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " qc <input.gff3>\n"
        << "\n"
        << "Top-level equivalent: " << prog << " <input.gff3> --qc\n"
        << "\n"
        << "QC checks:\n"
        << "  duplicate_id      Repeated ID attributes.\n"
        << "  invalid_attribute_value  Empty tag=value attribute value.\n"
        << "  invalid_attribute_escape  Unescaped ampersand or double quote in attributes.\n"
        << "  invalid_dbxref    Malformed Dbxref attribute.\n"
        << "  invalid_source    Empty source column; use . when source is unknown.\n"
        << "  missing_derives_from  Derives_from points to an absent ID.\n"
        << "  invalid_range     start greater than end.\n"
        << "  invalid_gap       Malformed Gap attribute.\n"
        << "  invalid_is_circular  Is_circular value is not true.\n"
        << "  invalid_ontology_term  Malformed Ontology_term attribute.\n"
        << "  invalid_percent_encoding  Malformed percent encoding in attributes.\n"
        << "  invalid_target    Malformed Target attribute.\n"
        << "  missing_parent    Parent points to an absent ID.\n"
        << "  parent_cycle      Parent relationships contain a cycle.\n"
        << "  child_outside_parent  Child coordinates outside parent coordinates.\n";
}

static bool append_unique(GffData& out, std::unordered_set<int>& seen, const GffRecord& rec) {
    if (!seen.insert(rec.line_idx).second) {
        return false;
    }
    out.append(rec);
    return true;
}

static std::string trim_copy(std::string_view value) {
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

static std::vector<std::string> split_attr_keys(std::string_view keys) {
    std::vector<std::string> result;
    size_t pos = 0;
    while (pos <= keys.size()) {
        const size_t comma = keys.find(',', pos);
        const size_t end = (comma == std::string_view::npos) ? keys.size() : comma;
        const auto key = trim_copy(keys.substr(pos, end - pos));
        if (!key.empty()) {
            result.push_back(key);
        }
        if (comma == std::string_view::npos) {
            break;
        }
        pos = comma + 1;
    }
    return result;
}

static std::optional<char> parse_strand_filter(std::string_view value) {
    if (value.size() != 1) {
        return std::nullopt;
    }
    const char strand = value[0];
    if (strand == '+' || strand == '-' || strand == '.' || strand == '?') {
        return strand;
    }
    return std::nullopt;
}

static std::optional<char> parse_phase_filter(std::string_view value) {
    if (value.size() != 1) {
        return std::nullopt;
    }
    const char phase = value[0];
    if (phase == '0' || phase == '1' || phase == '2' || phase == '.') {
        return phase;
    }
    return std::nullopt;
}

static std::optional<std::optional<double>> parse_score_filter(std::string_view value) {
    if (value == ".") {
        return std::optional<double>{};
    }
    if (value.empty()) {
        return std::nullopt;
    }
    try {
        size_t parsed = 0;
        const double score = std::stod(std::string{value}, &parsed);
        if (parsed == value.size() && std::isfinite(score)) {
            return std::optional<double>{score};
        }
    } catch (const std::exception&) {
    }
    return std::nullopt;
}

static std::string join_values(const std::vector<std::string>& values) {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << values[i];
    }
    return out.str();
}

static std::vector<char*> argv_from(std::vector<std::string>& args) {
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (auto& arg : args) {
        argv.push_back(arg.data());
    }
    return argv;
}

static void print_qc_row(std::ostream& out,
                         const char* severity,
                         const char* code,
                         int line_idx,
                         const std::string& id,
                         const std::string& message) {
    out << severity << '\t'
        << code << '\t'
        << line_idx << '\t'
        << id << '\t'
        << message << '\n';
}

static void print_percent_encoding_qc(std::ostream& out,
                                      std::string_view value,
                                      std::string_view field,
                                      int line_idx,
                                      const std::string& id) {
    if (const auto error = percent_encoding_error(value)) {
        const auto message = field.empty() ? *error : std::string{field} + " " + *error;
        print_qc_row(out, "error", "invalid_percent_encoding", line_idx, id, message);
    }
}

static int run_query(int argc, char* argv[], const char* prog) {
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        query_usage(prog);
        return 0;
    }

    if (argc < 2) {
        query_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    std::vector<std::string> ids;
    std::string name;
    std::string id_list_file;
    std::string region_str;
    std::string nearest_region_str;
    std::string feature_type;
    std::string summary_format;
    std::vector<std::string> output_attrs;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    bool include_children = false;
    bool include_parents = false;
    bool include_model = false;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            ids.push_back(*value);
        } else if (arg == "--name") {
            auto value = require_value("--name");
            if (!value) return 1;
            name = *value;
        } else if (arg == "--ids" || arg == "--id-list") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            id_list_file = *value;
        } else if (arg == "--region") {
            auto value = require_value("--region");
            if (!value) return 1;
            region_str = *value;
        } else if (arg == "--nearest" || arg == "--nearest-gene") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            nearest_region_str = *value;
        } else if (arg == "--type") {
            auto value = require_value("--type");
            if (!value) return 1;
            feature_type = *value;
        } else if (arg == "--where" || arg == "--attr") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            const auto equal_pos = value->find('=');
            if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value->size()) {
                std::cerr << "Error: " << arg << " expects KEY=VALUE\n";
                return 1;
            }
            attr_filters.emplace_back(value->substr(0, equal_pos), value->substr(equal_pos + 1));
        } else if (arg == "--output-attrs" || arg == "--out-attrs" || arg == "--attrs") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            const auto keys = split_attr_keys(*value);
            if (keys.empty()) {
                std::cerr << "Error: " << arg << " expects a comma-separated list of keys\n";
                return 1;
            }
            output_attrs.insert(output_attrs.end(), keys.begin(), keys.end());
        } else if (arg == "--summary" || arg == "--summary-format") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            summary_format = *value;
            if (summary_format != "tsv" && summary_format != "json") {
                std::cerr << "Error: " << arg << " expects tsv or json\n";
                return 1;
            }
        } else if (arg == "-C" || arg == "--children" || arg == "--include-children") {
            include_children = true;
        } else if (arg == "--parents" || arg == "--include-parents") {
            include_parents = true;
        } else if (arg == "--model" || arg == "--gene-model") {
            include_model = true;
        } else if (arg == "-h" || arg == "--help") {
            query_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown query option " << arg << '\n';
            query_usage(prog);
            return 1;
        }
    }

    if ((include_children || include_parents || include_model) && ids.empty() && id_list_file.empty() && name.empty() && attr_filters.empty() && nearest_region_str.empty()) {
        std::cerr << "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest\n";
        return 1;
    }

    if (!id_list_file.empty()) {
        std::ifstream in{id_list_file};
        if (!in.is_open()) {
            std::cerr << "Error: cannot open " << id_list_file << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty()) {
                ids.push_back(line);
            }
        }
    }

    gffsub::AnnotationIndex index = gffsub::AnnotationIndex::from_gff3(input_file);
    GffData result;
    std::unordered_set<int> seen;
    std::vector<SummaryRow> summary_rows;
    const bool emit_summary = !summary_format.empty() || !output_attrs.empty();

    auto add_summary = [&](const std::string& query_id, const std::string& matched_by, const GffRecord& rec) {
        if (emit_summary) {
            auto row = make_summary_row(index, query_id, matched_by, rec);
            if (!output_attrs.empty()) {
                row.attrs = extract_output_attrs(rec.attr_raw, output_attrs);
            }
            summary_rows.push_back(std::move(row));
        }
    };

    auto add_match = [&](const GffRecord& rec, const std::string& query_id, const std::string& matched_by) {
        if (include_model) {
            const auto model = index.gene_model(record_id(rec));
            if (model) {
                for (const auto& model_rec : model->records) {
                    if (feature_type.empty() || model_rec.type == feature_type) {
                        if (append_unique(result, seen, model_rec)) {
                            add_summary(query_id, "model", model_rec);
                        }
                    }
                }
                return;
            }
        }
        if (feature_type.empty() || rec.type == feature_type) {
            if (append_unique(result, seen, rec)) {
                add_summary(query_id, matched_by, rec);
            }
        }
        if (include_parents && rec.id) {
            std::vector<GffRecord> stack = index.parents_of(*rec.id);
            std::unordered_set<int> visited_parents;
            while (!stack.empty()) {
                const auto parent = stack.back();
                stack.pop_back();
                if (!visited_parents.insert(parent.line_idx).second) {
                    continue;
                }
                if (feature_type.empty() || parent.type == feature_type) {
                    if (append_unique(result, seen, parent)) {
                        add_summary(query_id, "parent", parent);
                    }
                }
                if (parent.id) {
                    for (const auto& grandparent : index.parents_of(*parent.id)) {
                        stack.push_back(grandparent);
                    }
                }
            }
        }
        if (include_children && rec.id) {
            for (const auto& child : index.descendants_of(*rec.id)) {
                if (feature_type.empty() || child.type == feature_type) {
                    if (append_unique(result, seen, child)) {
                        add_summary(query_id, "child", child);
                    }
                }
            }
        }
    };

    for (const auto& id : ids) {
        const auto rec = index.find_by_id(id);
        if (rec) {
            add_match(*rec, id, "ID");
        } else if (emit_summary) {
            auto row = make_not_found_row(id, "ID");
            row.attrs.assign(output_attrs.size(), "");
            summary_rows.push_back(std::move(row));
        }
    }

    if (!name.empty()) {
        const auto rec = index.find_gene(name);
        if (rec) {
            add_match(*rec, name, infer_gene_match_key(index, name, *rec));
        } else if (emit_summary) {
            auto row = make_not_found_row(name, "name");
            row.attrs.assign(output_attrs.size(), "");
            summary_rows.push_back(std::move(row));
        }
    }

    if (!region_str.empty()) {
        const auto region = parse_region(region_str);
        if (!region) {
            std::cerr << "Error: invalid region format " << region_str << '\n';
            return 1;
        }
        for (const auto& rec : index.overlap(region->seqid, region->start, region->end)) {
            add_match(rec, region_str, "region");
        }
    }
    if (!nearest_region_str.empty()) {
        const auto region = parse_region(nearest_region_str);
        if (!region) {
            std::cerr << "Error: invalid nearest region format " << nearest_region_str << '\n';
            return 1;
        }
        const auto rec = index.nearest_gene(region->seqid, region->start, region->end);
        if (rec) {
            add_match(*rec, nearest_region_str, "nearest");
        } else if (emit_summary) {
            auto row = make_not_found_row(nearest_region_str, "nearest");
            row.attrs.assign(output_attrs.size(), "");
            summary_rows.push_back(std::move(row));
        }
    }

    for (const auto& [key, value] : attr_filters) {
        bool matched = false;
        for (const auto& rec : index.with_attribute(key, value)) {
            matched = true;
            add_match(rec, key + "=" + value, key);
        }
        if (!matched && emit_summary) {
            auto row = make_not_found_row(key + "=" + value, key);
            row.attrs.assign(output_attrs.size(), "");
            summary_rows.push_back(std::move(row));
        }
    }

    std::sort(result.records.begin(), result.records.end(),
              [](const GffRecord& lhs, const GffRecord& rhs) {
                  return lhs.line_idx < rhs.line_idx;
              });
    if (summary_format == "json") {
        print_summary_json(std::cout, summary_rows, output_attrs);
    } else if (emit_summary) {
        print_summary_tsv(std::cout, summary_rows, output_attrs);
    } else {
        print_gff3(std::cout, result);
    }
    return 0;
}

static int run_window(int argc, char* argv[], const char* prog) {
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        window_usage(prog);
        return 0;
    }

    if (argc < 2) {
        window_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    std::string id;
    int64_t upstream = 0;
    int64_t downstream = 0;
    bool strand_aware = false;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };
        auto parse_non_negative = [&](const char* option, const std::string& text, int64_t& value) -> bool {
            try {
                value = std::stoll(text);
            } catch (const std::exception&) {
                std::cerr << "Error: " << option << " must be a non-negative integer\n";
                return false;
            }
            if (value < 0) {
                std::cerr << "Error: " << option << " must be non-negative\n";
                return false;
            }
            return true;
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            if (!id.empty()) {
                std::cerr << "Error: window requires exactly one --id\n";
                return 1;
            }
            id = *value;
        } else if (arg == "--up" || arg == "--upstream") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            if (!parse_non_negative(arg.c_str(), *value, upstream)) {
                return 1;
            }
        } else if (arg == "--down" || arg == "--downstream") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            if (!parse_non_negative(arg.c_str(), *value, downstream)) {
                return 1;
            }
        } else if (arg == "--strand-aware") {
            strand_aware = true;
        } else if (arg == "-h" || arg == "--help") {
            window_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown window option " << arg << '\n';
            window_usage(prog);
            return 1;
        }
    }

    if (id.empty()) {
        std::cerr << "Error: window requires --id\n";
        return 1;
    }

    const auto index = gffsub::AnnotationIndex::from_gff3(input_file);
    auto target = index.find_by_id(id);
    if (!target) {
        target = index.find_gene(id);
    }
    if (!target) {
        std::cerr << "Error: cannot find " << id << '\n';
        return 1;
    }

    const auto region = window_region(*target, upstream, downstream, strand_aware);
    GffData result;
    for (const auto& rec : index.overlap(region.seqid, region.start, region.end)) {
        result.append(rec);
    }
    print_gff3(std::cout, result);
    return 0;
}

static int run_qc(int argc, char* argv[], const char* prog) {
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        qc_usage(prog);
        return 0;
    }

    if (argc != 2) {
        qc_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    const auto qc_parse = parse_qc_records(input_file);
    if (!qc_parse.opened) {
        std::cerr << "Error: cannot parse " << input_file << '\n';
        return 1;
    }
    const auto& data = qc_parse.data;
    const auto directive_result = parse_directives(input_file);

    std::unordered_map<std::string, const GffRecord*> by_id;
    std::unordered_map<std::string, std::vector<const GffRecord*>> records_by_id;
    std::unordered_map<std::string, std::vector<std::string>> parents_by_id;
    std::unordered_set<std::string> circular_seqids;
    for (const auto& rec : data.records) {
        const auto attrs = parse_attributes(rec.attr_raw);
        if (rec.type == "region" && has_circular_true(attrs)) {
            circular_seqids.insert(rec.seqid);
        }
        if (rec.id) {
            records_by_id[*rec.id].push_back(&rec);
            by_id.emplace(*rec.id, &rec);

            const auto parent_it = attrs.find("Parent");
            if (parent_it != attrs.end()) {
                parents_by_id.emplace(*rec.id, parent_it->second);
            }
        }
    }
    const auto parent_cycle_ids = find_parent_cycle_ids(parents_by_id);

    std::cout << "severity\tcode\tline_idx\tid\tmessage\n";

    for (const auto& [id, records] : records_by_id) {
        if (records.size() > 1 && !allowed_discontinuous_id(records)) {
            print_qc_row(std::cout, "error", "duplicate_id", -1, id, "ID appears more than once");
        }
    }
    for (const auto& issue : directive_result.issues) {
        print_qc_row(std::cout, "error", issue.code.c_str(), issue.line_idx, ".", issue.message);
    }
    for (const auto& issue : qc_parse.issues) {
        print_qc_row(std::cout, "error", issue.code.c_str(), issue.line_idx, ".", issue.message);
    }

    for (const auto& rec : data.records) {
        const std::string id = record_id(rec);
        const bool has_valid_coordinates = rec.start >= 1 && rec.end >= 1 && rec.start <= rec.end;
        const auto attrs = parse_attributes(rec.attr_raw);
        print_percent_encoding_qc(std::cout, rec.attr_raw, "", rec.line_idx, id);
        print_percent_encoding_qc(std::cout, rec.source, "source", rec.line_idx, id);
        print_percent_encoding_qc(std::cout, rec.type, "feature type", rec.line_idx, id);
        if (const auto error = attribute_escape_error(rec.attr_raw)) {
            print_qc_row(std::cout, "error", "invalid_attribute_escape", rec.line_idx, id, *error);
        }
        if (const auto error = seqid_syntax_error(rec.seqid)) {
            print_qc_row(std::cout, "error", "invalid_seqid", rec.line_idx, id, *error);
        }
        if (const auto error = feature_type_syntax_error(rec.type)) {
            print_qc_row(std::cout, "error", "invalid_feature_type", rec.line_idx, id, *error);
        }
        if (rec.start < 1 || rec.end < 1) {
            print_qc_row(std::cout, "error", "invalid_coordinate", rec.line_idx, id,
                         "start and end must be positive 1-based coordinates");
        }
        if (rec.start > rec.end) {
            print_qc_row(std::cout, "error", "invalid_range", rec.line_idx, id, "start is greater than end");
        }
        if (rec.score && !std::isfinite(*rec.score)) {
            print_qc_row(std::cout, "error", "invalid_score", rec.line_idx, id,
                         "score must be a finite floating point number or .");
        }
        const auto sequence_region_it = directive_result.sequence_regions.find(rec.seqid);
        if (has_valid_coordinates && sequence_region_it != directive_result.sequence_regions.end()) {
            const auto& sequence_region = sequence_region_it->second;
            const bool outside_end = rec.end > sequence_region.end && circular_seqids.find(rec.seqid) == circular_seqids.end();
            if (rec.start < sequence_region.start || outside_end) {
                print_qc_row(std::cout, "error", "outside_sequence_region", rec.line_idx, id,
                             "feature is outside ##sequence-region " + rec.seqid);
            }
        }
        if (rec.strand != '+' && rec.strand != '-' && rec.strand != '.' && rec.strand != '?') {
            print_qc_row(std::cout, "error", "invalid_strand", rec.line_idx, id,
                         std::string{"strand "} + rec.strand + " must be +, -, ., or ?");
        }
        if (rec.type != "CDS" && rec.phase != '.') {
            print_qc_row(std::cout, "error", "invalid_phase", rec.line_idx, id,
                         std::string{"non-CDS phase "} + rec.phase + " must be .");
        }
        if (rec.type == "CDS" && rec.phase != '0' && rec.phase != '1' && rec.phase != '2') {
            print_qc_row(std::cout, "error", "invalid_cds_phase", rec.line_idx, id,
                         std::string{"CDS phase "} + rec.phase + " must be 0, 1, or 2");
        }
        if (has_invalid_is_circular(attrs)) {
            print_qc_row(std::cout, "error", "invalid_is_circular", rec.line_idx, id,
                         "Is_circular value must be true");
        }

        const auto target_it = attrs.find("Target");
        if (target_it != attrs.end()) {
            for (const auto& target : target_it->second) {
                if (const auto error = target_attribute_error(target)) {
                    print_qc_row(std::cout, "error", "invalid_target", rec.line_idx, id, *error);
                }
            }
        }
        const auto gap_it = attrs.find("Gap");
        if (gap_it != attrs.end()) {
            for (const auto& gap : gap_it->second) {
                if (const auto error = gap_attribute_error(gap)) {
                    print_qc_row(std::cout, "error", "invalid_gap", rec.line_idx, id, *error);
                }
            }
        }
        const auto dbxref_it = attrs.find("Dbxref");
        if (dbxref_it != attrs.end()) {
            for (const auto& dbxref : dbxref_it->second) {
                if (const auto error = database_accession_error("Dbxref", dbxref)) {
                    print_qc_row(std::cout, "error", "invalid_dbxref", rec.line_idx, id, *error);
                }
            }
        }
        const auto ontology_term_it = attrs.find("Ontology_term");
        if (ontology_term_it != attrs.end()) {
            for (const auto& ontology_term : ontology_term_it->second) {
                if (const auto error = database_accession_error("Ontology_term", ontology_term)) {
                    print_qc_row(std::cout, "error", "invalid_ontology_term", rec.line_idx, id, *error);
                }
            }
        }
        const auto derives_it = attrs.find("Derives_from");
        if (derives_it != attrs.end()) {
            for (const auto& source_id : derives_it->second) {
                if (by_id.find(source_id) == by_id.end()) {
                    print_qc_row(std::cout, "error", "missing_derives_from", rec.line_idx, id,
                                 "Derives_from " + source_id + " was not found");
                }
            }
        }

        const auto parent_it = attrs.find("Parent");
        if (rec.id && parent_cycle_ids.find(*rec.id) != parent_cycle_ids.end()) {
            print_qc_row(std::cout, "error", "parent_cycle", rec.line_idx, id,
                         "Parent relationships contain a cycle");
        }
        if (parent_it == attrs.end()) {
            continue;
        }
        std::unordered_set<std::string> seen_parents;
        for (const auto& parent_id : parent_it->second) {
            if (!seen_parents.insert(parent_id).second) {
                print_qc_row(std::cout, "error", "duplicate_parent", rec.line_idx, id,
                             "Parent " + parent_id + " appears more than once");
                continue;
            }
            const auto parent_record_it = by_id.find(parent_id);
            if (parent_record_it == by_id.end()) {
                print_qc_row(std::cout, "error", "missing_parent", rec.line_idx, id,
                             "Parent " + parent_id + " was not found");
                continue;
            }

            const auto& parent = *parent_record_it->second;
            if (rec.seqid != parent.seqid || rec.start < parent.start || rec.end > parent.end) {
                print_qc_row(std::cout, "warning", "child_outside_parent", rec.line_idx, id,
                             "child is outside Parent " + parent_id);
            }
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc > 1 && std::string(argv[1]) == "query") {
        return run_query(argc - 1, argv + 1, argv[0]);
    }
    if (argc > 1 && std::string(argv[1]) == "window") {
        return run_window(argc - 1, argv + 1, argv[0]);
    }
    if (argc > 1 && std::string(argv[1]) == "qc") {
        return run_qc(argc - 1, argv + 1, argv[0]);
    }

    std::vector<std::string> ids;
    std::string id_list_file;
    std::string name;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    std::vector<GrepFilter> grep_filters;
    std::string grep_file;
    std::string grep_field;
    bool grep_file_regex = false;
    std::vector<ExprNode> include_expr_filters;
    std::vector<ExprNode> exclude_expr_filters;
    bool invert_grep = false;
    bool ignore_case = false;
    bool include_children = false;
    bool include_parents = false;
    bool include_model = false;
    std::vector<std::string> output_attrs;
    std::string summary_format;
    std::string upstream_arg;
    std::string downstream_arg;
    bool strand_aware = false;
    bool do_qc = false;
    std::string region_str;
    std::string nearest_region_str;
    std::string seqid_filter;
    std::string source_filter;
    std::optional<std::optional<double>> score_filter;
    std::optional<char> strand_filter;
    std::optional<char> phase_filter;
    std::string bed_file;
    std::string feature;
    bool do_longest = false;
    size_t num_threads = 6;
    bool threads_set = false;
    std::string output_format = "gff3";
    std::string output_file;

    enum {
        OPT_ID = 1000,
        OPT_ID_LIST,
        OPT_NAME,
        OPT_ATTR,
        OPT_OUTPUT_ATTRS,
        OPT_SUMMARY_FORMAT,
        OPT_PARENTS,
        OPT_MODEL,
        OPT_NEAREST,
        OPT_UPSTREAM,
        OPT_DOWNSTREAM,
        OPT_STRAND_AWARE,
        OPT_QC,
        OPT_SEQID,
        OPT_SOURCE,
        OPT_SCORE,
        OPT_STRAND_FILTER,
        OPT_PHASE,
        OPT_GREP,
        OPT_GREP_REGEX,
        OPT_GREP_FILE,
        OPT_GREP_FIELD,
        OPT_GREP_FILE_REGEX,
        OPT_INCLUDE_EXPR,
        OPT_EXCLUDE_EXPR,
        OPT_INVERT_MATCH,
        OPT_IGNORE_CASE
    };
    static struct option long_options[] = {
        {"id",            required_argument, nullptr, OPT_ID},
        {"ids",           required_argument, nullptr, OPT_ID_LIST},
        {"id-list",       required_argument, nullptr, OPT_ID_LIST},
        {"name",          required_argument, nullptr, OPT_NAME},
        {"where",         required_argument, nullptr, OPT_ATTR},
        {"attr",          required_argument, nullptr, OPT_ATTR},
        {"grep",          required_argument, nullptr, OPT_GREP},
        {"grep-regex",    required_argument, nullptr, OPT_GREP_REGEX},
        {"grep-file",     required_argument, nullptr, OPT_GREP_FILE},
        {"grep-field",    required_argument, nullptr, OPT_GREP_FIELD},
        {"grep-file-regex", no_argument,     nullptr, OPT_GREP_FILE_REGEX},
        {"include-expr",  required_argument, nullptr, OPT_INCLUDE_EXPR},
        {"exclude-expr",  required_argument, nullptr, OPT_EXCLUDE_EXPR},
        {"invert-match",  no_argument,       nullptr, OPT_INVERT_MATCH},
        {"ignore-case",   no_argument,       nullptr, OPT_IGNORE_CASE},
        {"output-attrs",  required_argument, nullptr, OPT_OUTPUT_ATTRS},
        {"out-attrs",     required_argument, nullptr, OPT_OUTPUT_ATTRS},
        {"attrs",         required_argument, nullptr, OPT_OUTPUT_ATTRS},
        {"summary",       required_argument, nullptr, OPT_SUMMARY_FORMAT},
        {"summary-format", required_argument, nullptr, OPT_SUMMARY_FORMAT},
        {"parents",       no_argument,       nullptr, OPT_PARENTS},
        {"include-parents", no_argument,      nullptr, OPT_PARENTS},
        {"model",         no_argument,       nullptr, OPT_MODEL},
        {"gene-model",    no_argument,       nullptr, OPT_MODEL},
        {"nearest",       required_argument, nullptr, OPT_NEAREST},
        {"nearest-gene",  required_argument, nullptr, OPT_NEAREST},
        {"up",            required_argument, nullptr, OPT_UPSTREAM},
        {"upstream",      required_argument, nullptr, OPT_UPSTREAM},
        {"down",          required_argument, nullptr, OPT_DOWNSTREAM},
        {"downstream",    required_argument, nullptr, OPT_DOWNSTREAM},
        {"strand-aware",  no_argument,       nullptr, OPT_STRAND_AWARE},
        {"qc",            no_argument,       nullptr, OPT_QC},
        {"seqid",         required_argument, nullptr, OPT_SEQID},
        {"source",        required_argument, nullptr, OPT_SOURCE},
        {"score",         required_argument, nullptr, OPT_SCORE},
        {"strand",        required_argument, nullptr, OPT_STRAND_FILTER},
        {"phase",         required_argument, nullptr, OPT_PHASE},
        {"children",      no_argument,       nullptr, 'C'},
        {"include-children", no_argument,     nullptr, 'C'},
        {"region",        required_argument, nullptr, 'r'},
        {"bed",           required_argument, nullptr, 'b'},
        {"feature",       required_argument, nullptr, 'f'},
        {"type",          required_argument, nullptr, 'f'},
        {"longest",       no_argument,       nullptr, 'L'},
        {"threads",       required_argument, nullptr, '@'},
        {"format",        required_argument, nullptr, 't'},
        {"output-format", required_argument, nullptr, 't'},
        {"output",        required_argument, nullptr, 'o'},
        {"help",          no_argument,       nullptr, 'h'},
        {nullptr,        0,                 nullptr, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "r:b:f:CL@:t:o:hI:E:v", long_options, &option_index)) != -1) {
        switch (opt) {
            case OPT_ID: ids.emplace_back(optarg); break;
            case OPT_ID_LIST: id_list_file = optarg; break;
            case OPT_NAME: name = optarg; break;
            case OPT_ATTR: {
                const std::string value{optarg};
                const auto equal_pos = value.find('=');
                if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value.size()) {
                    std::cerr << "Error: --" << long_options[option_index].name << " expects KEY=VALUE\n";
                    return 1;
                }
                attr_filters.emplace_back(value.substr(0, equal_pos), value.substr(equal_pos + 1));
                break;
            }
            case OPT_GREP:
            case OPT_GREP_REGEX: {
                const auto parsed = parse_field_pattern(optarg);
                if (!parsed) {
                    std::cerr << "Error: --" << long_options[option_index].name << " expects FIELD:PATTERN\n";
                    return 1;
                }
                GrepFilter filter;
                filter.field = parsed->first;
                filter.pattern = parsed->second;
                filter.use_regex = opt == OPT_GREP_REGEX;
                grep_filters.push_back(std::move(filter));
                break;
            }
            case OPT_GREP_FILE:
                grep_file = optarg;
                break;
            case OPT_GREP_FIELD:
                grep_field = optarg;
                break;
            case OPT_GREP_FILE_REGEX:
                grep_file_regex = true;
                break;
            case OPT_INCLUDE_EXPR:
            case 'I': {
                std::string error;
                if (!parse_expr_filters(optarg, include_expr_filters, error)) {
                    std::cerr << "Error: invalid include expression: " << error << '\n';
                    return 1;
                }
                break;
            }
            case OPT_EXCLUDE_EXPR:
            case 'E': {
                std::string error;
                if (!parse_expr_filters(optarg, exclude_expr_filters, error)) {
                    std::cerr << "Error: invalid exclude expression: " << error << '\n';
                    return 1;
                }
                break;
            }
            case OPT_INVERT_MATCH:
            case 'v':
                invert_grep = true;
                break;
            case OPT_IGNORE_CASE:
                ignore_case = true;
                break;
            case OPT_OUTPUT_ATTRS: {
                const auto keys = split_attr_keys(optarg);
                if (keys.empty()) {
                    std::cerr << "Error: --out-attrs expects a comma-separated list of keys\n";
                    return 1;
                }
                output_attrs.insert(output_attrs.end(), keys.begin(), keys.end());
                break;
            }
            case OPT_SUMMARY_FORMAT:
                summary_format = optarg;
                if (summary_format != "tsv" && summary_format != "json") {
                    std::cerr << "Error: --" << long_options[option_index].name << " expects tsv or json\n";
                    return 1;
                }
                break;
            case OPT_PARENTS: include_parents = true; break;
            case OPT_MODEL: include_model = true; break;
            case OPT_NEAREST: nearest_region_str = optarg; break;
            case 'C': include_children = true; break;
            case OPT_UPSTREAM: upstream_arg = optarg; break;
            case OPT_DOWNSTREAM: downstream_arg = optarg; break;
            case OPT_STRAND_AWARE: strand_aware = true; break;
            case OPT_QC: do_qc = true; break;
            case OPT_SEQID: seqid_filter = optarg; break;
            case OPT_SOURCE: source_filter = optarg; break;
            case OPT_SCORE: {
                score_filter = parse_score_filter(optarg);
                if (!score_filter) {
                    std::cerr << "Error: --score expects a finite floating point number or .\n";
                    return 1;
                }
                break;
            }
            case OPT_STRAND_FILTER: {
                strand_filter = parse_strand_filter(optarg);
                if (!strand_filter) {
                    std::cerr << "Error: --strand expects one of +, -, ., ?\n";
                    return 1;
                }
                break;
            }
            case OPT_PHASE: {
                phase_filter = parse_phase_filter(optarg);
                if (!phase_filter) {
                    std::cerr << "Error: --phase expects one of 0, 1, 2, .\n";
                    return 1;
                }
                break;
            }
            case 'r': region_str = optarg; break;
            case 'b': bed_file = optarg; break;
            case 'f': feature = optarg; break;
            case 'L': do_longest = true; break;
            case '@': {
                threads_set = true;
                size_t t = 0;
                try {
                    t = std::stoul(optarg);
                } catch (const std::exception&) {
                    std::cerr << "Error: --threads must be a non-negative integer\n";
                    return 1;
                }
                if (t == 0) t = 1;
                if (t > 256) t = 256; // cap to prevent over-subscription
                num_threads = t;
                break;
            }
            case 't': output_format = optarg; break;
            case 'o': output_file = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }

    if (optind >= argc) {
        usage(argv[0]);
        return 1;
    }

    if (!grep_file.empty()) {
        if (grep_field.empty()) {
            std::cerr << "Error: --grep-file requires --grep-field\n";
            return 1;
        }
        std::ifstream in{grep_file};
        if (!in.is_open()) {
            std::cerr << "Error: cannot open " << grep_file << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(in, line)) {
            const auto pattern = trim_copy(line);
            if (pattern.empty()) {
                continue;
            }
            GrepFilter filter;
            filter.field = grep_field;
            filter.pattern = pattern;
            filter.use_regex = grep_file_regex;
            grep_filters.push_back(std::move(filter));
        }
    } else if (!grep_field.empty() || grep_file_regex) {
        std::cerr << "Error: --grep-field and --grep-file-regex require --grep-file\n";
        return 1;
    }
    if (invert_grep && grep_filters.empty()) {
        std::cerr << "Error: --invert-match requires --grep, --grep-regex, or --grep-file\n";
        return 1;
    }
    {
        std::string error;
        if (!compile_filter_regexes(grep_filters, include_expr_filters, exclude_expr_filters, ignore_case, error)) {
            std::cerr << "Error: " << error << '\n';
            return 1;
        }
    }

    const bool has_query_style_selector = !ids.empty() || !id_list_file.empty() || !name.empty() || !attr_filters.empty() || !nearest_region_str.empty();
    if ((include_children || include_parents || include_model) && !has_query_style_selector) {
        std::cerr << "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest\n";
        return 1;
    }

    // Validate output format
    OutputFormat fmt = OutputFormat::GFF3;
    if (output_format == "gtf") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf2") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf3") {
        fmt = OutputFormat::GTF3;
    } else if (output_format == "bed") {
        fmt = OutputFormat::BED;
    } else if (output_format == "gff3") {
        fmt = OutputFormat::GFF3;
    } else {
        std::cerr << "Error: unknown output format " << output_format << '\n';
        std::cerr << "Supported formats: gff3, gtf2, gtf3, bed\n";
        return 1;
    }

    std::string input_file = argv[optind];

    if (do_qc) {
        if (!ids.empty() || !id_list_file.empty() || !name.empty() || !attr_filters.empty() || !nearest_region_str.empty() || include_children || include_parents || include_model ||
            !output_attrs.empty() || !summary_format.empty() || !upstream_arg.empty() || !downstream_arg.empty() ||
            strand_aware || score_filter || strand_filter || phase_filter || !region_str.empty() || !seqid_filter.empty() || !source_filter.empty() || !bed_file.empty() || !feature.empty() || do_longest ||
            !grep_filters.empty() || !grep_file.empty() || !grep_field.empty() || grep_file_regex || !include_expr_filters.empty() || !exclude_expr_filters.empty() || invert_grep || ignore_case ||
            output_format != "gff3" || !output_file.empty()) {
            std::cerr << "Error: --qc only supports the input file\n";
            return 1;
        }

        std::vector<std::string> qc_args{"qc", input_file};
        auto qc_argv = argv_from(qc_args);
        return run_qc(static_cast<int>(qc_argv.size()), qc_argv.data(), argv[0]);
    }

    if (!upstream_arg.empty() || !downstream_arg.empty() || strand_aware) {
        if (ids.size() != 1) {
            std::cerr << "Error: window shortcut requires exactly one --id\n";
            return 1;
        }
        if (!id_list_file.empty() || !name.empty() || !attr_filters.empty() || !nearest_region_str.empty() || include_children || include_parents || include_model ||
            !output_attrs.empty() || !summary_format.empty() || !region_str.empty() || !bed_file.empty() ||
            !seqid_filter.empty() || !source_filter.empty() || score_filter || strand_filter || phase_filter || !feature.empty() || do_longest ||
            !grep_filters.empty() || !grep_file.empty() || !grep_field.empty() || grep_file_regex || !include_expr_filters.empty() || !exclude_expr_filters.empty() || invert_grep || ignore_case ||
            output_format != "gff3" || !output_file.empty()) {
            std::cerr << "Error: window shortcut only supports --id, --up/--upstream, --down/--downstream, and --strand-aware\n";
            return 1;
        }

        std::vector<std::string> window_args{"window", input_file};
        for (const auto& id : ids) {
            window_args.push_back("--id");
            window_args.push_back(id);
        }
        if (!upstream_arg.empty()) {
            window_args.push_back("--up");
            window_args.push_back(upstream_arg);
        }
        if (!downstream_arg.empty()) {
            window_args.push_back("--down");
            window_args.push_back(downstream_arg);
        }
        if (strand_aware) {
            window_args.push_back("--strand-aware");
        }

        auto window_argv = argv_from(window_args);
        return run_window(static_cast<int>(window_argv.size()), window_argv.data(), argv[0]);
    }

    const bool can_dispatch_summary_to_query = seqid_filter.empty() && source_filter.empty() && !score_filter && !strand_filter && !phase_filter && bed_file.empty() && !do_longest && !threads_set &&
                                               grep_filters.empty() && grep_file.empty() && grep_field.empty() && !grep_file_regex && include_expr_filters.empty() && exclude_expr_filters.empty() && !invert_grep && !ignore_case &&
                                               output_format == "gff3" && output_file.empty();
    const bool can_dispatch_default_selector_to_query = can_dispatch_summary_to_query && region_str.empty();

    if (!summary_format.empty() || !output_attrs.empty()) {
        if (!can_dispatch_summary_to_query) {
            std::cerr << "Error: --summary/--summary-format/--out-attrs only supports query-style selectors; "
                      << "do not combine with --seqid, --source, --score, --strand, --phase, --bed, --longest, --threads, --format/--output-format, or --output\n";
            return 1;
        }
    }

    if (!summary_format.empty() || !output_attrs.empty() ||
        (has_query_style_selector && can_dispatch_default_selector_to_query)) {
        std::vector<std::string> query_args{"query", input_file};
        for (const auto& id : ids) {
            query_args.push_back("--id");
            query_args.push_back(id);
        }
        if (!id_list_file.empty()) {
            query_args.push_back("--ids");
            query_args.push_back(id_list_file);
        }
        if (!name.empty()) {
            query_args.push_back("--name");
            query_args.push_back(name);
        }
        if (!region_str.empty()) {
            query_args.push_back("--region");
            query_args.push_back(region_str);
        }
        if (!nearest_region_str.empty()) {
            query_args.push_back("--nearest");
            query_args.push_back(nearest_region_str);
        }
        if (!feature.empty()) {
            query_args.push_back("--type");
            query_args.push_back(feature);
        }
        for (const auto& [key, value] : attr_filters) {
            query_args.push_back("--where");
            query_args.push_back(key + "=" + value);
        }
        if (include_children) {
            query_args.push_back("--children");
        }
        if (include_parents) {
            query_args.push_back("--parents");
        }
        if (include_model) {
            query_args.push_back("--model");
        }
        if (!summary_format.empty()) {
            query_args.push_back("--summary");
            query_args.push_back(summary_format);
        }
        if (!output_attrs.empty()) {
            query_args.push_back("--out-attrs");
            query_args.push_back(join_values(output_attrs));
        }

        auto query_argv = argv_from(query_args);
        return run_query(static_cast<int>(query_argv.size()), query_argv.data(), argv[0]);
    }

    if (!id_list_file.empty()) {
        std::ifstream in{id_list_file};
        if (!in.is_open()) {
            std::cerr << "Error: cannot open " << id_list_file << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty()) {
                ids.push_back(line);
            }
        }
    }

    GffData data;
    IdIndex idx;

    // Parse input file
    if (parse_file(input_file, data, idx, InputFormat::GFF3) != 0) {
        std::cerr << "Error: cannot parse " << input_file << '\n';
        return 1;
    }

    if (!ids.empty() || !name.empty() || !attr_filters.empty() || !nearest_region_str.empty()) {
        const auto index = gffsub::AnnotationIndex::from_gff3(input_file);
        std::unordered_set<int> selected_lines;
        auto add_selected = [&](const GffRecord& rec) {
            const bool newly_selected = selected_lines.insert(rec.line_idx).second;
            if (!newly_selected) {
                return;
            }
            if (include_model) {
                const auto model = index.gene_model(record_id(rec));
                if (model) {
                    for (const auto& model_rec : model->records) {
                        selected_lines.insert(model_rec.line_idx);
                    }
                    return;
                }
            }
            if (include_parents && rec.id) {
                std::vector<GffRecord> stack = index.parents_of(*rec.id);
                std::unordered_set<int> visited_parents;
                while (!stack.empty()) {
                    const auto parent = stack.back();
                    stack.pop_back();
                    if (!visited_parents.insert(parent.line_idx).second) {
                        continue;
                    }
                    selected_lines.insert(parent.line_idx);
                    if (parent.id) {
                        for (const auto& grandparent : index.parents_of(*parent.id)) {
                            stack.push_back(grandparent);
                        }
                    }
                }
            }
            if (include_children && rec.id) {
                for (const auto& child : index.descendants_of(*rec.id)) {
                    selected_lines.insert(child.line_idx);
                }
            }
        };
        for (const auto& id : ids) {
            const auto rec = index.find_by_id(id);
            if (rec) {
                add_selected(*rec);
            }
        }
        if (!name.empty()) {
            const auto rec = index.find_gene(name);
            if (rec) {
                add_selected(*rec);
            }
        }
        for (const auto& [key, value] : attr_filters) {
            for (const auto& rec : index.with_attribute(key, value)) {
                add_selected(rec);
            }
        }
        if (!nearest_region_str.empty()) {
            const auto nearest_region = parse_region(nearest_region_str);
            if (!nearest_region) {
                std::cerr << "Error: invalid nearest region format " << nearest_region_str << '\n';
                return 1;
            }
            const auto rec = index.nearest_gene(nearest_region->seqid, nearest_region->start, nearest_region->end);
            if (rec) {
                add_selected(*rec);
            }
        }
        for (auto& rec : data.records) {
            rec.kept = selected_lines.count(rec.line_idx) > 0;
        }
    }

    Region region{"", 0, 0};
    std::optional<Region> parsed_region;

    // Apply region filters
    if (!region_str.empty()) {
        parsed_region = parse_region(region_str);
        if (!parsed_region) {
            std::cerr << "Error: invalid region format " << region_str << '\n';
            return 1;
        }
        region = *parsed_region;
        filter_by_region(data, region);
    }

    if (!bed_file.empty()) {
        filter_by_regions_from_file(data, bed_file);
    }

    if (!seqid_filter.empty()) {
        filter_by_seqid(data, seqid_filter);
    }

    if (!source_filter.empty()) {
        filter_by_source(data, source_filter);
    }

    if (score_filter) {
        filter_by_score(data, *score_filter);
    }

    if (strand_filter) {
        filter_by_strand(data, *strand_filter);
    }

    if (phase_filter) {
        filter_by_phase(data, *phase_filter);
    }

    if (!grep_filters.empty()) {
        filter_by_grep(data, grep_filters, invert_grep);
    }

    if (!include_expr_filters.empty()) {
        filter_by_expr(data, include_expr_filters, true);
    }

    if (!exclude_expr_filters.empty()) {
        filter_by_expr(data, exclude_expr_filters, false);
    }

    // Apply feature filters
    if (do_longest) {
        filter_longest_isoform(data, idx, feature, num_threads);
    } else if (!feature.empty()) {
        filter_by_feature(data, feature);
    }

    // Output
    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!output_file.empty()) {
        out_file.open(output_file);
        if (!out_file.is_open()) {
            std::cerr << "Error: cannot open " << output_file << '\n';
            return 1;
        }
        out = &out_file;
    }

    switch (fmt) {
        case OutputFormat::GFF3: print_gff3(*out, data); break;
        case OutputFormat::GTF2: print_gtf(*out, data, fmt); break;
        case OutputFormat::GTF3: print_gtf3(*out, data); break;
        case OutputFormat::BED:  print_bed(*out, data); break;
    }

    return 0;
}
