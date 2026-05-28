#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

static bool write_test_annotation(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) {
        return false;
    }

    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=gene0001;Name=ABC1;gene_id=G1;locus_tag=Locus1;Alias=ABC-1,LegacyABC;Dbxref=GeneID:123\n"
        << "chr1\tsrc\tmRNA\t100\t400\t42.5\t+\t.\tID=tx1;Parent=gene0001;Name=ABC1.1\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tID=exon1;Parent=tx1\n"
        << "chr1\tsrc\tCDS\t150\t170\t.\t+\t0\tID=cds1;Parent=tx1\n"
        << "chr1\tsrc\tCDS\t200\t220\t.\t+\t1\tID=cds2;Parent=tx1\n"
        << "chr1\tsrc\tCDS\t230\t250\t.\t+\t2\tID=cds3;Parent=tx1\n"
        << "chr1\tsrc\tgene\t600\t700\t.\t-\t.\tID=gene0002;Name=XYZ1;biotype=protein_coding;Note=transposon-like\n"
        << "chr2\tother\tgene\t100\t200\t.\t+\t.\tID=gene0003;Name=CHR2\n"
        << "chr2\tother\texon\t250\t280\t.\t.\t.\tID=exon2\n"
        << "chr2\tother\tmRNA\t300\t380\t.\t+\t.\tID=orphan_tx\n"
        << "chr2\tother\texon\t320\t360\t.\t+\t.\tID=orphan_exon;Parent=orphan_tx\n";
    return true;
}

static bool write_qc_annotation(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) {
        return false;
    }

    out << "##gff-version 3\n"
        << "##sequence-region chr1 1 1000\n"
        << "##sequence-region chr_bad 100 1\n"
        << "##sequence-region chr_signed +1 500\n"
        << "##sequence-region chr1 1 1200\n"
        << "##sequence-region chr2 1 500\n"
        << "##sequence-region chr#bad 1 500\n"
        << "##sequence-region chr_extra 1 500 extra\n"
        << "##gff-version 3 extra\n"
        << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=dup_gene\n"
        << "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=dup_gene\n"
        << "chr1\tsrc\tmRNA\t410\t480\t.\t+\t.\tID=disc_tx;Parent=disc_gene\n"
        << "chr1\tsrc\tCDS\t410\t430\t.\t+\t0\tID=disc_cds;Parent=disc_tx\n"
        << "chr1\tsrc\tCDS\t450\t480\t.\t+\t2\tID=disc_cds;Parent=disc_tx\n"
        << "chr1\tsrc\tpolypeptide\t410\t430\t.\t+\t.\tID=disc_poly;Derives_from=disc_tx\n"
        << "chr1\tsrc\tpolypeptide\t450\t480\t.\t+\t.\tID=disc_poly;Derives_from=disc_tx\n"
        << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\n"
        << "chr1\tsrc\tgene\t150\t180\tnan\t+\t.\tID=bad_score\n"
        << "chr1\t\tgene\t181\t190\t.\t+\t.\tID=bad_empty_source\n"
        << "chr1\tsrc\tgene\tabc\t180\t.\t+\t.\tID=bad_coordinate_text\n"
        << "chr1\tsrc\tgene\t+180\t190\t.\t+\t.\tID=bad_coordinate_sign\n"
        << "chr1\tsrc\tgene\t180\t190\t.\t+\t.\tID\n"
        << "chr1\tsrc\tgene\t180\t190\t.\t+\t.\t\n"
        << "chr1\tsrc\tgene\t181\t190\t.\t+\t.\tID=bad_attr_equals;Note=A=B\n"
        << "chr1\tsrc\tgene\t190\t195\t.\t+\t.\tID=dup_attr;Name=A;Name=B\n"
        << "chr1\tsrc\tgene\t196\t198\t.\t+\t.\tID=bad_multivalue;Name=A,B\n"
        << "chr1\tsrc\tgene\t199\t199\t.\t+\t.\tID=bad_empty_attr;Name=\n"
        << "chr1\tsrc\tgene\t201\t205\t.\t+\t.\tID=bad_empty_list_item;Alias=A,\n"
        << "bad seq\tsrc\tgene\t200\t210\t.\t+\t.\tID=bad_seqid\n"
        << "chr#bad\tsrc\tgene\t211\t215\t.\t+\t.\tID=bad_seqid_char\n"
        << "chr1\tsrc\t.\t220\t230\t.\t+\t.\tID=bad_type\n"
        << "chr1\tsrc\tSO:abc\t240\t250\t.\t+\t.\tID=bad_so_type\n"
        << "chr1\tsrc\tSO:123\t251\t255\t.\t+\t.\tID=bad_so_width\n"
        << "chr1\tsrc\tgene\t0\t50\t.\t+\t.\tID=bad_coordinate\n"
        << "chr1\tsrc\tmRNA\t500\t600\t.\t+\t.\tID=orphan_tx;Parent=missing_gene\n"
        << "chr1\tsrc\tmRNA\t610\t620\t.\t+\t.\tID=extra_column_child;Parent=missing_extra_parent\textra\n"
        << "chr1\tsrc\tgene\t700\t800\t.\tx\t.\tID=bad_strand\n"
        << "chr1\tsrc\tgene\t720\t730\t.\t\t.\tID=bad_empty_strand\n"
        << "chr1\tsrc\tgene\t731\t740\t.\t+\t\tID=bad_empty_phase\n"
        << "chr1\tsrc\tgene\t900\t950\t.\t+\tx\tID=bad_phase\n"
        << "chr1\tsrc\texon\t700\t710\t.\t+\t1\tID=bad_non_cds_phase\n"
        << "chr2\tsrc\tgene\t490\t510\t.\t+\t.\tID=outside_region\n"
        << "chr1\tsrc\tregion\t1\t1200\t.\t+\t.\tID=circular_region;Is_circular=true\n"
        << "chr1\tsrc\tgene\t1190\t1210\t.\t+\t.\tID=circular_gene\n"
        << "chr1\tsrc\tregion\t1\t1200\t.\t+\t.\tID=bad_circular;Is_circular=false\n"
        << "chr1\tsrc\tmRNA\t260\t270\t.\t+\t.\tID=cycle_a;Parent=cycle_b\n"
        << "chr1\tsrc\tmRNA\t260\t270\t.\t+\t.\tID=cycle_b;Parent=cycle_a\n"
        << "chr1\tsrc\tmatch_part\t280\t290\t.\t+\t.\tID=bad_target;Target=read1 20 10\n"
        << "chr1\tsrc\tmatch_part\t300\t320\t.\t+\t.\tID=bad_gap;Gap=M8 X3\n"
        << "chr1\tsrc\tpolypeptide\t330\t360\t.\t+\t.\tID=bad_derives;Derives_from=missing_transcript\n"
        << "chr1\tsrc\tgene\t370\t380\t.\t+\t.\tID=bad_dbxref;Dbxref=GeneID\n"
        << "chr1\tsrc\tgene\t390\t395\t.\t+\t.\tID=bad_ontology;Ontology_term=SO\n"
        << "chr1\tsrc\tgene\t396\t398\t.\t+\t.\tID=bad_percent;Note=bad%escape\n"
        << "chr1\tsrc%zz\tgene\t406\t407\t.\t+\t.\tID=bad_source_percent\n"
        << "chr1\tsrc\tgene%zz\t408\t409\t.\t+\t.\tID=bad_type_percent\n"
        << "chr1\tsrc\tgene\t399\t400\t.\t+\t.\tID=bad_escape;Note=A&B\n"
        << "chr1\tsrc\tgene\t401\t405\t.\t+\t.\tID=bad_quote;Note=A\"B\n"
        << "chr1\tsrc\texon\t120\t130\t.\t+\t.\tID=dup_parent;Parent=dup_gene,dup_gene\n"
        << "chr1\tsrc\tCDS\t520\t540\t.\t+\t.\tID=bad_cds_phase;Parent=orphan_tx\n"
        << ">chr1\n"
        << "chr1\tsrc\tmRNA\t800\t810\t.\t+\t.\tID=fasta_feature;Parent=fasta_parent\n";
    return true;
}

static bool write_qc_digit_coordinate_annotation(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) {
        return false;
    }

    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t+180\t190\t.\t+\t.\tID=bad_coordinate_sign\n";
    return true;
}

static std::string read_file(const std::string& path) {
    std::ifstream in{path};
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

static int run_command(const std::string& command) {
    const int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "command failed: " << command << '\n';
    }
    return status;
}

static bool contains(const std::string& text, const std::string& needle) {
    return text.find(needle) != std::string::npos;
}

static int require_contains(const std::string& path, const std::string& needle) {
    const auto text = read_file(path);
    if (!contains(text, needle)) {
        std::cerr << "missing '" << needle << "' in " << path << '\n';
        return 1;
    }
    return 0;
}

static int require_not_contains(const std::string& path, const std::string& needle) {
    const auto text = read_file(path);
    if (contains(text, needle)) {
        std::cerr << "unexpected '" << needle << "' in " << path << '\n';
        return 1;
    }
    return 0;
}

static int compare_files(const std::string& lhs_path, const std::string& rhs_path) {
    const auto lhs = read_file(lhs_path);
    const auto rhs = read_file(rhs_path);
    if (lhs != rhs) {
        std::cerr << "output mismatch: " << lhs_path << " vs " << rhs_path << '\n';
        return 1;
    }
    return 0;
}

static void cleanup_outputs() {
    std::remove("cli_selector_smoke.gff3");
    std::remove("cli_selector_qc.gff3");
    std::remove("cli_selector_ids.txt");
    std::remove("cli_selector_patterns.txt");
    std::remove("selector_id.gff3");
    std::remove("selector_query_id.gff3");
    std::remove("selector_attr_id.gff3");
    std::remove("selector_where_id.gff3");
    std::remove("selector_grep_id.gff3");
    std::remove("selector_grep_regex_id.gff3");
    std::remove("selector_grep_regex_seqid.gff3");
    std::remove("selector_grep_file.gff3");
    std::remove("selector_grep_invert.gff3");
    std::remove("selector_grep_ignore_case.gff3");
    std::remove("selector_include_expr_biotype.gff3");
    std::remove("selector_include_expr_logic.gff3");
    std::remove("selector_include_expr_quoted_regex.gff3");
    std::remove("selector_include_expr_numeric.gff3");
    std::remove("selector_include_expr_score.gff3");
    std::remove("selector_exclude_expr_note.gff3");
    std::remove("selector_expr_bad.out");
    std::remove("selector_expr_bad.err");
    std::remove("selector_grep_bad.out");
    std::remove("selector_grep_bad.err");
    std::remove("selector_invert_bad.out");
    std::remove("selector_invert_bad.err");
    std::remove("selector_id_list.gff3");
    std::remove("selector_id_list_verbose.gff3");
    std::remove("selector_query_id_list.gff3");
    std::remove("selector_query_id_list_verbose.gff3");
    std::remove("selector_id_list_children.gff3");
    std::remove("selector_id_list_children_alias.gff3");
    std::remove("selector_id_list_children_verbose.gff3");
    std::remove("selector_query_id_list_children.gff3");
    std::remove("selector_query_id_list_children_alias.gff3");
    std::remove("selector_query_id_list_children_verbose.gff3");
    std::remove("selector_name.gff3");
    std::remove("selector_query_name.gff3");
    std::remove("selector_alias.gff3");
    std::remove("selector_dbxref.gff3");
    std::remove("selector_name_summary.tsv");
    std::remove("selector_name_summary_verbose.tsv");
    std::remove("selector_gene_id_summary.tsv");
    std::remove("selector_locus_tag_summary.tsv");
    std::remove("selector_alias_summary.tsv");
    std::remove("selector_dbxref_summary.tsv");
    std::remove("selector_parent.gff3");
    std::remove("selector_parent_attr.gff3");
    std::remove("selector_query_parent.gff3");
    std::remove("selector_query_parent_attr.gff3");
    std::remove("selector_children.gff3");
    std::remove("selector_children_alias.gff3");
    std::remove("selector_query_children.gff3");
    std::remove("selector_query_children_alias.gff3");
    std::remove("selector_query_children_short.gff3");
    std::remove("selector_children_short.gff3");
    std::remove("selector_children_type.gff3");
    std::remove("selector_query_children_type.gff3");
    std::remove("selector_parents.gff3");
    std::remove("selector_parents_alias.gff3");
    std::remove("selector_query_parents.gff3");
    std::remove("selector_query_parents_alias.gff3");
    std::remove("selector_parents_gene.gff3");
    std::remove("selector_query_parents_gene.gff3");
    std::remove("selector_parents_children.gff3");
    std::remove("selector_parents_bad.out");
    std::remove("selector_parents_bad.err");
    std::remove("selector_query_parents_bad.out");
    std::remove("selector_query_parents_bad.err");
    std::remove("selector_model.gff3");
    std::remove("selector_gene_model_alias.gff3");
    std::remove("selector_query_model.gff3");
    std::remove("selector_query_gene_model_alias.gff3");
    std::remove("selector_model_gene.gff3");
    std::remove("selector_query_model_gene.gff3");
    std::remove("selector_model_orphan_children.gff3");
    std::remove("selector_model_bad.out");
    std::remove("selector_model_bad.err");
    std::remove("selector_query_model_bad.out");
    std::remove("selector_query_model_bad.err");
    std::remove("selector_nearest.gff3");
    std::remove("selector_nearest_alias.gff3");
    std::remove("selector_query_nearest.gff3");
    std::remove("selector_query_nearest_alias.gff3");
    std::remove("selector_nearest_overlap.gff3");
    std::remove("selector_nearest_children.gff3");
    std::remove("selector_nearest_seqid_keep.gff3");
    std::remove("selector_nearest_seqid_drop.gff3");
    std::remove("selector_nearest_summary.tsv");
    std::remove("selector_nearest_not_found.tsv");
    std::remove("selector_nearest_bad.out");
    std::remove("selector_nearest_bad.err");
    std::remove("selector_query_nearest_bad.out");
    std::remove("selector_query_nearest_bad.err");
    std::remove("selector_region_intersection.gff3");
    std::remove("selector_seqid_chr2.gff3");
    std::remove("selector_seqid_gene.gff3");
    std::remove("selector_seqid_id_intersection.gff3");
    std::remove("selector_source_other.gff3");
    std::remove("selector_source_gene.gff3");
    std::remove("selector_seqid_source_intersection.gff3");
    std::remove("selector_score_value.gff3");
    std::remove("selector_score_missing.gff3");
    std::remove("selector_score_feature.gff3");
    std::remove("selector_score_bad.out");
    std::remove("selector_score_bad.err");
    std::remove("selector_score_nan.out");
    std::remove("selector_score_nan.err");
    std::remove("selector_strand_minus.gff3");
    std::remove("selector_strand_dot.gff3");
    std::remove("selector_strand_gene.gff3");
    std::remove("selector_strand_bad.out");
    std::remove("selector_strand_bad.err");
    std::remove("selector_phase_zero.gff3");
    std::remove("selector_phase_dot.gff3");
    std::remove("selector_phase_cds.gff3");
    std::remove("selector_phase_bad.out");
    std::remove("selector_phase_bad.err");
    std::remove("selector_help.txt");
    std::remove("selector_bed_short.bed");
    std::remove("selector_bed_format.bed");
    std::remove("selector_bed_output_format.bed");
    std::remove("selector_window_top.gff3");
    std::remove("selector_window_top_short.gff3");
    std::remove("selector_window_command.gff3");
    std::remove("selector_window_command_short.gff3");
    std::remove("selector_qc_top.tsv");
    std::remove("selector_qc_command.tsv");
    std::remove("selector_qc_digit_coordinate.gff3");
    std::remove("selector_qc_digit_coordinate.tsv");
    std::remove("selector_query_help.txt");
    std::remove("selector_window_help.txt");
    std::remove("selector_qc_help.txt");
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: cli_selector_smoke <gffsub-executable>\n";
        return 2;
    }

    const std::string exe = std::string{"\""} + argv[1] + "\"";
    const std::string exe_raw = argv[1];
    const std::string gff{"cli_selector_smoke.gff3"};
    if (!write_test_annotation(gff)) {
        std::cerr << "cannot write test annotation\n";
        return 1;
    }

    if (run_command(exe + " query --help > selector_query_help.txt 2>&1") != 0 ||
        require_contains("selector_query_help.txt", "Most workflows can use the top-level form") != 0 ||
        require_contains("selector_query_help.txt", "--ids FILE") != 0 ||
        require_contains("selector_query_help.txt", "Verbose alias for --ids") != 0 ||
        require_contains("selector_query_help.txt", "--where KEY=VALUE") != 0 ||
        require_contains("selector_query_help.txt", "--children") != 0 ||
        require_contains("selector_query_help.txt", "--parents") != 0 ||
        require_contains("selector_query_help.txt", "--model") != 0 ||
        require_contains("selector_query_help.txt", "--nearest REGION") != 0 ||
        require_contains("selector_query_help.txt", "--summary FMT") != 0 ||
        require_contains("selector_query_help.txt", "Verbose alias for --summary") != 0) {
        return 1;
    }
    if (run_command(exe + " window --help > selector_window_help.txt 2>&1") != 0 ||
        require_contains("selector_window_help.txt", "Top-level equivalent") != 0 ||
        require_contains("selector_window_help.txt", "--up N") != 0 ||
        require_contains("selector_window_help.txt", "--down N") != 0) {
        return 1;
    }
    if (run_command(exe + " qc --help > selector_qc_help.txt 2>&1") != 0 ||
        require_contains("selector_qc_help.txt", "Top-level equivalent") != 0) {
        return 1;
    }
    if (run_command(exe + " --help > selector_help.txt 2>&1") != 0 ||
        require_contains("selector_help.txt", "--format FMT") != 0 ||
        require_contains("selector_help.txt", "--where KEY=VALUE") != 0 ||
        require_contains("selector_help.txt", "--seqid SEQID") != 0 ||
        require_contains("selector_help.txt", "--source SOURCE") != 0 ||
        require_contains("selector_help.txt", "--score SCORE") != 0 ||
        require_contains("selector_help.txt", "--strand STRAND") != 0 ||
        require_contains("selector_help.txt", "--phase PHASE") != 0 ||
        require_contains("selector_help.txt", "--grep FIELD:PATTERN") != 0 ||
        require_contains("selector_help.txt", "--grep-regex FIELD:REGEX") != 0 ||
        require_contains("selector_help.txt", "--grep-file FILE --grep-field FIELD") != 0 ||
        require_contains("selector_help.txt", "--include-expr EXPR") != 0 ||
        require_contains("selector_help.txt", "--exclude-expr EXPR") != 0 ||
        require_contains("selector_help.txt", "Invert --grep/--grep-regex/--grep-file matches") != 0 ||
        require_contains("selector_help.txt", "--model") != 0 ||
        require_contains("selector_help.txt", "--nearest REGION") != 0 ||
        require_contains("selector_help.txt", "--output-format remains as a verbose alias") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 > selector_id.gff3") != 0 ||
        require_contains("selector_id.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_id.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " query " + gff + " --id gene0001 > selector_query_id.gff3") != 0 ||
        compare_files("selector_id.gff3", "selector_query_id.gff3") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --attr ID=gene0001 > selector_attr_id.gff3") != 0 ||
        compare_files("selector_id.gff3", "selector_attr_id.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --where ID=gene0001 > selector_where_id.gff3") != 0 ||
        compare_files("selector_id.gff3", "selector_where_id.gff3") != 0) {
        return 1;
    }

    if (run_command(exe_raw + " " + gff + " --grep ID:gene000 > selector_grep_id.gff3") != 0 ||
        require_contains("selector_grep_id.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_id.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_grep_id.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_id.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep-regex ID:gene000[12] > selector_grep_regex_id.gff3") != 0 ||
        require_contains("selector_grep_regex_id.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_regex_id.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_grep_regex_id.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_regex_id.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep-regex seqid:chr[0-9]+ -f gene > selector_grep_regex_seqid.gff3") != 0 ||
        require_contains("selector_grep_regex_seqid.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_regex_seqid.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_grep_regex_seqid.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_regex_seqid.gff3", "ID=tx1") != 0) {
        return 1;
    }
    {
        std::ofstream patterns{"cli_selector_patterns.txt"};
        if (!patterns.is_open()) {
            std::cerr << "cannot write grep pattern list\n";
            return 1;
        }
        patterns << "gene0001\n"
                 << "gene0003\n";
    }
    if (run_command(exe_raw + " " + gff + " --grep-file cli_selector_patterns.txt --grep-field ID > selector_grep_file.gff3") != 0 ||
        require_contains("selector_grep_file.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_file.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_file.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_grep_file.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep ID:gene000 -v > selector_grep_invert.gff3") != 0 ||
        require_contains("selector_grep_invert.gff3", "ID=tx1") != 0 ||
        require_contains("selector_grep_invert.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_grep_invert.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep name:abc1 --ignore-case > selector_grep_ignore_case.gff3") != 0 ||
        require_contains("selector_grep_ignore_case.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_ignore_case.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_grep_ignore_case.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"type==gene && attr.biotype==protein_coding\" > selector_include_expr_biotype.gff3") != 0 ||
        require_contains("selector_include_expr_biotype.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_include_expr_biotype.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_include_expr_biotype.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"(type==gene && attr.biotype==protein_coding) || !seqid==chr1\" > selector_include_expr_logic.gff3") != 0 ||
        require_contains("selector_include_expr_logic.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_include_expr_logic.gff3", "ID=gene0003") != 0 ||
        require_contains("selector_include_expr_logic.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_include_expr_logic.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_include_expr_logic.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"attr.ID~\\\"gene000[13]\\\"\" > selector_include_expr_quoted_regex.gff3") != 0 ||
        require_contains("selector_include_expr_quoted_regex.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_include_expr_quoted_regex.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_include_expr_quoted_regex.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_include_expr_quoted_regex.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"type==gene && length>=101\" > selector_include_expr_numeric.gff3") != 0 ||
        require_contains("selector_include_expr_numeric.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_include_expr_numeric.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_include_expr_numeric.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_include_expr_numeric.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"score==42.5\" > selector_include_expr_score.gff3") != 0 ||
        require_contains("selector_include_expr_score.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_include_expr_score.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -E \"attr.Note~transposon|retroelement\" -f gene > selector_exclude_expr_note.gff3") != 0 ||
        require_contains("selector_exclude_expr_note.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_exclude_expr_note.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_exclude_expr_note.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I type > selector_expr_bad.out 2> selector_expr_bad.err") == 0 ||
        require_contains("selector_expr_bad.err", "Error: invalid include expression") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep-regex ID:[ > selector_grep_bad.out 2> selector_grep_bad.err") == 0 ||
        require_contains("selector_grep_bad.err", "Error: invalid regex in --grep-regex") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -v > selector_invert_bad.out 2> selector_invert_bad.err") == 0 ||
        require_contains("selector_invert_bad.err", "Error: --invert-match requires --grep, --grep-regex, or --grep-file") != 0) {
        return 1;
    }

    const std::string id_list{"cli_selector_ids.txt"};
    {
        std::ofstream out{id_list};
        if (!out.is_open()) {
            std::cerr << "cannot write ID list\n";
            return 1;
        }
        out << "gene0001\n"
            << "gene0002\n";
    }
    if (run_command(exe + " " + gff + " --ids " + id_list + " > selector_id_list.gff3") != 0 ||
        run_command(exe + " " + gff + " --id-list " + id_list + " > selector_id_list_verbose.gff3") != 0 ||
        run_command(exe + " query " + gff + " --ids " + id_list + " > selector_query_id_list.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id-list " + id_list + " > selector_query_id_list_verbose.gff3") != 0 ||
        compare_files("selector_id_list.gff3", "selector_id_list_verbose.gff3") != 0 ||
        compare_files("selector_id_list.gff3", "selector_query_id_list.gff3") != 0 ||
        compare_files("selector_query_id_list.gff3", "selector_query_id_list_verbose.gff3") != 0 ||
        require_contains("selector_id_list.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_id_list.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --ids " + id_list + " --children > selector_id_list_children.gff3") != 0 ||
        run_command(exe + " " + gff + " --ids " + id_list + " --include-children > selector_id_list_children_alias.gff3") != 0 ||
        run_command(exe + " " + gff + " --id-list " + id_list + " --include-children > selector_id_list_children_verbose.gff3") != 0 ||
        run_command(exe + " query " + gff + " --ids " + id_list + " --children > selector_query_id_list_children.gff3") != 0 ||
        run_command(exe + " query " + gff + " --ids " + id_list + " --include-children > selector_query_id_list_children_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id-list " + id_list + " --include-children > selector_query_id_list_children_verbose.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_id_list_children_alias.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_id_list_children_verbose.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_query_id_list_children.gff3") != 0 ||
        compare_files("selector_query_id_list_children.gff3", "selector_query_id_list_children_alias.gff3") != 0 ||
        compare_files("selector_query_id_list_children.gff3", "selector_query_id_list_children_verbose.gff3") != 0 ||
        require_contains("selector_id_list_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_id_list_children.gff3", "ID=exon1") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name ABC1 > selector_name.gff3") != 0 ||
        require_contains("selector_name.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_name.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --name ABC1 > selector_query_name.gff3") != 0 ||
        compare_files("selector_name.gff3", "selector_query_name.gff3") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name LegacyABC > selector_alias.gff3") != 0 ||
        require_contains("selector_alias.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_alias.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name GeneID:123 > selector_dbxref.gff3") != 0 ||
        require_contains("selector_dbxref.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_dbxref.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name ABC1 --summary tsv > selector_name_summary.tsv") != 0 ||
        run_command(exe + " " + gff + " --name ABC1 --summary-format tsv > selector_name_summary_verbose.tsv") != 0 ||
        compare_files("selector_name_summary.tsv", "selector_name_summary_verbose.tsv") != 0 ||
        require_contains("selector_name_summary.tsv", "ABC1\tgene0001\tName") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name G1 --summary tsv > selector_gene_id_summary.tsv") != 0 ||
        require_contains("selector_gene_id_summary.tsv", "G1\tgene0001\tgene_id") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name Locus1 --summary tsv > selector_locus_tag_summary.tsv") != 0 ||
        require_contains("selector_locus_tag_summary.tsv", "Locus1\tgene0001\tlocus_tag") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name LegacyABC --summary tsv > selector_alias_summary.tsv") != 0 ||
        require_contains("selector_alias_summary.tsv", "LegacyABC\tgene0001\tAlias") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name GeneID:123 --summary tsv > selector_dbxref_summary.tsv") != 0 ||
        require_contains("selector_dbxref_summary.tsv", "GeneID:123\tgene0001\tDbxref") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --where Parent=tx1 > selector_parent.gff3") != 0 ||
        require_contains("selector_parent.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_parent.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --attr Parent=tx1 > selector_parent_attr.gff3") != 0 ||
        run_command(exe + " query " + gff + " --where Parent=tx1 > selector_query_parent.gff3") != 0 ||
        run_command(exe + " query " + gff + " --attr Parent=tx1 > selector_query_parent_attr.gff3") != 0 ||
        compare_files("selector_parent.gff3", "selector_parent_attr.gff3") != 0 ||
        compare_files("selector_parent.gff3", "selector_query_parent.gff3") != 0 ||
        compare_files("selector_query_parent.gff3", "selector_query_parent_attr.gff3") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --children > selector_children.gff3") != 0 ||
        require_contains("selector_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_children.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --include-children > selector_children_alias.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_children_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --children > selector_query_children.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_query_children.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --include-children > selector_query_children_alias.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_query_children_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 -C > selector_query_children_short.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_query_children_short.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 -C > selector_children_short.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_children_short.gff3") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --include-children -f mRNA > selector_children_type.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id gene0001 --include-children --type mRNA > selector_query_children_type.gff3") != 0 ||
        compare_files("selector_children_type.gff3", "selector_query_children_type.gff3") != 0 ||
        require_contains("selector_children_type.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_children_type.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --parents > selector_parents.gff3") != 0 ||
        require_contains("selector_parents.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_parents.gff3", "ID=tx1") != 0 ||
        require_contains("selector_parents.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_parents.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --include-parents > selector_parents_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --parents > selector_query_parents.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --include-parents > selector_query_parents_alias.gff3") != 0 ||
        compare_files("selector_parents.gff3", "selector_parents_alias.gff3") != 0 ||
        compare_files("selector_parents.gff3", "selector_query_parents.gff3") != 0 ||
        compare_files("selector_query_parents.gff3", "selector_query_parents_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --parents -f gene > selector_parents_gene.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --parents --type gene > selector_query_parents_gene.gff3") != 0 ||
        compare_files("selector_parents_gene.gff3", "selector_query_parents_gene.gff3") != 0 ||
        require_contains("selector_parents_gene.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_parents_gene.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_parents_gene.gff3", "ID=exon1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id tx1 --parents -C > selector_parents_children.gff3") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=exon1") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_parents_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --parents > selector_parents_bad.out 2> selector_parents_bad.err") == 0 ||
        require_contains("selector_parents_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --parents > selector_query_parents_bad.out 2> selector_query_parents_bad.err") == 0 ||
        require_contains("selector_query_parents_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --model > selector_model.gff3") != 0 ||
        require_contains("selector_model.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_model.gff3", "ID=tx1") != 0 ||
        require_contains("selector_model.gff3", "ID=exon1") != 0 ||
        require_contains("selector_model.gff3", "ID=cds1") != 0 ||
        require_contains("selector_model.gff3", "ID=cds2") != 0 ||
        require_contains("selector_model.gff3", "ID=cds3") != 0 ||
        require_not_contains("selector_model.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --gene-model > selector_gene_model_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --model > selector_query_model.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --gene-model > selector_query_gene_model_alias.gff3") != 0 ||
        compare_files("selector_model.gff3", "selector_gene_model_alias.gff3") != 0 ||
        compare_files("selector_model.gff3", "selector_query_model.gff3") != 0 ||
        compare_files("selector_query_model.gff3", "selector_query_gene_model_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --model -f CDS > selector_model_gene.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --model --type CDS > selector_query_model_gene.gff3") != 0 ||
        compare_files("selector_model_gene.gff3", "selector_query_model_gene.gff3") != 0 ||
        require_contains("selector_model_gene.gff3", "ID=cds1") != 0 ||
        require_contains("selector_model_gene.gff3", "ID=cds2") != 0 ||
        require_contains("selector_model_gene.gff3", "ID=cds3") != 0 ||
        require_not_contains("selector_model_gene.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_model_gene.gff3", "ID=exon1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id orphan_tx --model -C --seqid chr2 > selector_model_orphan_children.gff3") != 0 ||
        require_contains("selector_model_orphan_children.gff3", "ID=orphan_tx") != 0 ||
        require_contains("selector_model_orphan_children.gff3", "ID=orphan_exon") != 0 ||
        require_not_contains("selector_model_orphan_children.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --model > selector_model_bad.out 2> selector_model_bad.err") == 0 ||
        require_contains("selector_model_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --model > selector_query_model_bad.out 2> selector_query_model_bad.err") == 0 ||
        require_contains("selector_query_model_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 > selector_nearest.gff3") != 0 ||
        require_contains("selector_nearest.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_nearest.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest-gene chr1:450-500 > selector_nearest_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --nearest chr1:450-500 > selector_query_nearest.gff3") != 0 ||
        run_command(exe + " query " + gff + " --nearest-gene chr1:450-500 > selector_query_nearest_alias.gff3") != 0 ||
        compare_files("selector_nearest.gff3", "selector_nearest_alias.gff3") != 0 ||
        compare_files("selector_nearest.gff3", "selector_query_nearest.gff3") != 0 ||
        compare_files("selector_query_nearest.gff3", "selector_query_nearest_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:610-620 > selector_nearest_overlap.gff3") != 0 ||
        require_contains("selector_nearest_overlap.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_nearest_overlap.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 -C > selector_nearest_children.gff3") != 0 ||
        require_contains("selector_nearest_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_nearest_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_nearest_children.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_nearest_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 --seqid chr1 > selector_nearest_seqid_keep.gff3") != 0 ||
        require_contains("selector_nearest_seqid_keep.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_nearest_seqid_keep.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 --seqid chr2 > selector_nearest_seqid_drop.gff3") != 0 ||
        require_not_contains("selector_nearest_seqid_drop.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_nearest_seqid_drop.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 --summary tsv > selector_nearest_summary.tsv") != 0 ||
        require_contains("selector_nearest_summary.tsv", "chr1:450-500\tgene0001\tnearest") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr9:1-100 --summary tsv > selector_nearest_not_found.tsv") != 0 ||
        require_contains("selector_nearest_not_found.tsv", "chr9:1-100\t\tnearest") != 0 ||
        require_contains("selector_nearest_not_found.tsv", "not_found") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1-450-500 > selector_nearest_bad.out 2> selector_nearest_bad.err") == 0 ||
        require_contains("selector_nearest_bad.err", "Error: invalid nearest region format chr1-450-500") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --nearest chr1-450-500 > selector_query_nearest_bad.out 2> selector_query_nearest_bad.err") == 0 ||
        require_contains("selector_query_nearest_bad.err", "Error: invalid nearest region format chr1-450-500") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --region chr1:130-140 > selector_region_intersection.gff3") != 0 ||
        require_contains("selector_region_intersection.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --seqid chr2 > selector_seqid_chr2.gff3") != 0 ||
        require_contains("selector_seqid_chr2.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_seqid_chr2.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --seqid chr1 -f gene > selector_seqid_gene.gff3") != 0 ||
        require_contains("selector_seqid_gene.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_seqid_gene.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_seqid_gene.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_seqid_gene.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --seqid chr2 > selector_seqid_id_intersection.gff3") != 0 ||
        require_not_contains("selector_seqid_id_intersection.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_seqid_id_intersection.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --source other > selector_source_other.gff3") != 0 ||
        require_contains("selector_source_other.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_source_other.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --source src -f gene > selector_source_gene.gff3") != 0 ||
        require_contains("selector_source_gene.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_source_gene.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_source_gene.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_source_gene.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --seqid chr2 --source src > selector_seqid_source_intersection.gff3") != 0 ||
        require_not_contains("selector_seqid_source_intersection.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_seqid_source_intersection.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score 42.5 > selector_score_value.gff3") != 0 ||
        require_contains("selector_score_value.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_score_value.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score . > selector_score_missing.gff3") != 0 ||
        require_contains("selector_score_missing.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_score_missing.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score . -f gene > selector_score_feature.gff3") != 0 ||
        require_contains("selector_score_feature.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_score_feature.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_score_feature.gff3", "ID=exon1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score abc > selector_score_bad.out 2> selector_score_bad.err") == 0 ||
        require_contains("selector_score_bad.err", "Error: --score expects a finite floating point number or .") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score nan > selector_score_nan.out 2> selector_score_nan.err") == 0 ||
        require_contains("selector_score_nan.err", "Error: --score expects a finite floating point number or .") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand - > selector_strand_minus.gff3") != 0 ||
        require_contains("selector_strand_minus.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_strand_minus.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand . > selector_strand_dot.gff3") != 0 ||
        require_contains("selector_strand_dot.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_strand_dot.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand - -f gene > selector_strand_gene.gff3") != 0 ||
        require_contains("selector_strand_gene.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_strand_gene.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_strand_gene.gff3", "ID=exon2") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand forward > selector_strand_bad.out 2> selector_strand_bad.err") == 0 ||
        require_contains("selector_strand_bad.err", "Error: --strand expects one of +, -, ., ?") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase 0 > selector_phase_zero.gff3") != 0 ||
        require_contains("selector_phase_zero.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_phase_zero.gff3", "ID=cds2") != 0 ||
        require_not_contains("selector_phase_zero.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase . > selector_phase_dot.gff3") != 0 ||
        require_contains("selector_phase_dot.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_phase_dot.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_phase_dot.gff3", "ID=cds1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase 1 -f CDS > selector_phase_cds.gff3") != 0 ||
        require_contains("selector_phase_cds.gff3", "ID=cds2") != 0 ||
        require_not_contains("selector_phase_cds.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_phase_cds.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase 3 > selector_phase_bad.out 2> selector_phase_bad.err") == 0 ||
        require_contains("selector_phase_bad.err", "Error: --phase expects one of 0, 1, 2, .") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " -r chr1:100-400 -t bed > selector_bed_short.bed") != 0 ||
        run_command(exe + " " + gff + " -r chr1:100-400 --format bed > selector_bed_format.bed") != 0 ||
        run_command(exe + " " + gff + " -r chr1:100-400 --output-format bed > selector_bed_output_format.bed") != 0 ||
        compare_files("selector_bed_short.bed", "selector_bed_format.bed") != 0 ||
        compare_files("selector_bed_short.bed", "selector_bed_output_format.bed") != 0 ||
        require_contains("selector_bed_format.bed", "chr1\t99\t400\tgene0001") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --upstream 50 --downstream 10 > selector_window_top.gff3") != 0 ||
        run_command(exe + " " + gff + " --id gene0001 --up 50 --down 10 > selector_window_top_short.gff3") != 0 ||
        run_command(exe + " window " + gff + " --id gene0001 --upstream 50 --downstream 10 > selector_window_command.gff3") != 0 ||
        run_command(exe + " window " + gff + " --id gene0001 --up 50 --down 10 > selector_window_command_short.gff3") != 0 ||
        compare_files("selector_window_top.gff3", "selector_window_top_short.gff3") != 0 ||
        compare_files("selector_window_top.gff3", "selector_window_command.gff3") != 0 ||
        compare_files("selector_window_command.gff3", "selector_window_command_short.gff3") != 0) {
        return 1;
    }

    const std::string qc_gff{"cli_selector_qc.gff3"};
    if (!write_qc_annotation(qc_gff)) {
        std::cerr << "cannot write QC test annotation\n";
        return 1;
    }
    if (run_command(exe + " " + qc_gff + " --qc > selector_qc_top.tsv") != 0 ||
        run_command(exe + " qc " + qc_gff + " > selector_qc_command.tsv") != 0 ||
        compare_files("selector_qc_top.tsv", "selector_qc_command.tsv") != 0 ||
        require_contains("selector_qc_top.tsv", "duplicate_id") != 0 ||
        require_not_contains("selector_qc_top.tsv", "disc_cds\tID appears more than once") != 0 ||
        require_not_contains("selector_qc_top.tsv", "disc_poly\tID appears more than once") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_gff_version") != 0 ||
        require_contains("selector_qc_top.tsv", "##gff-version must declare exactly one version beginning with 3") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_column_count") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_sequence_region") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid ##sequence-region coordinates for chr_signed") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid ##sequence-region seqid chr#bad") != 0 ||
        require_contains("selector_qc_top.tsv", "malformed ##sequence-region directive") != 0 ||
        require_contains("selector_qc_top.tsv", "duplicate_sequence_region") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_score") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_source") != 0 ||
        require_contains("selector_qc_top.tsv", "source field must not be empty; use . when source is unknown") != 0 ||
        require_contains("selector_qc_top.tsv", "attribute tag Name must not contain comma-separated values") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_attribute_syntax") != 0 ||
        require_contains("selector_qc_top.tsv", "attributes must be semicolon-separated tag=value fields") != 0 ||
        require_contains("selector_qc_top.tsv", "attributes field must be . or semicolon-separated tag=value fields") != 0 ||
        require_contains("selector_qc_top.tsv", "duplicate_attribute_tag") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_attribute_multivalue") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_attribute_value") != 0 ||
        require_contains("selector_qc_top.tsv", "attribute tag Alias must have a non-empty value") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_seqid") != 0 ||
        require_contains("selector_qc_top.tsv", "seqid contains unescaped character #") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_feature_type") != 0 ||
        require_contains("selector_qc_top.tsv", "Sequence Ontology accession must be SO: followed by seven digits") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_coordinate") != 0 ||
        require_contains("selector_qc_top.tsv", "start and end must be integer 1-based coordinates") != 0 ||
        require_contains("selector_qc_top.tsv", "outside_sequence_region") != 0 ||
        require_not_contains("selector_qc_top.tsv", "circular_region\tfeature is outside ##sequence-region") != 0 ||
        require_not_contains("selector_qc_top.tsv", "circular_gene\tfeature is outside ##sequence-region") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_is_circular") != 0 ||
        require_contains("selector_qc_top.tsv", "duplicate_parent") != 0 ||
        require_contains("selector_qc_top.tsv", "parent_cycle") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_target") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_gap") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_dbxref") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_ontology_term") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_percent_encoding") != 0 ||
        require_contains("selector_qc_top.tsv", "source percent escapes must be % followed by two hex digits") != 0 ||
        require_contains("selector_qc_top.tsv", "feature type percent escapes must be % followed by two hex digits") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_attribute_escape") != 0 ||
        require_contains("selector_qc_top.tsv", "double quote in attributes must be percent-escaped as %22") != 0 ||
        require_contains("selector_qc_top.tsv", "missing_derives_from") != 0 ||
        require_contains("selector_qc_top.tsv", "missing_parent") != 0 ||
        require_contains("selector_qc_top.tsv", "Parent missing_extra_parent was not found") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_strand") != 0 ||
        require_contains("selector_qc_top.tsv", "strand field must be +, -, ., or ?") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_phase") != 0 ||
        require_contains("selector_qc_top.tsv", "phase field must be ., 0, 1, or 2") != 0 ||
        require_contains("selector_qc_top.tsv", "non-CDS phase 1 must be .") != 0 ||
        require_contains("selector_qc_top.tsv", "invalid_cds_phase") != 0) {
        return 1;
    }
    if (require_not_contains("selector_qc_top.tsv", "fasta_parent") != 0) {
        return 1;
    }
    const std::string qc_digit_coordinate_gff{"selector_qc_digit_coordinate.gff3"};
    if (!write_qc_digit_coordinate_annotation(qc_digit_coordinate_gff)) {
        std::cerr << "cannot write digit coordinate QC test annotation\n";
        return 1;
    }
    if (run_command(exe + " " + qc_digit_coordinate_gff + " --qc > selector_qc_digit_coordinate.tsv") != 0 ||
        require_contains("selector_qc_digit_coordinate.tsv", "invalid_coordinate") != 0 ||
        require_contains("selector_qc_digit_coordinate.tsv", "start and end must be integer 1-based coordinates") != 0) {
        return 1;
    }

    cleanup_outputs();
    std::cout << "cli_selector_smoke OK\n";
    return 0;
}
