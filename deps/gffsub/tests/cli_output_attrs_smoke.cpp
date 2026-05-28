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
        << "chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=gene0001;Name=ABC1;Alias=ABC-1;Dbxref=GeneID:123\n"
        << "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene0001;Name=ABC1.1\n";
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

static int expect_command_failure(const std::string& command) {
    const int status = std::system(command.c_str());
    if (status == 0) {
        std::cerr << "command unexpectedly succeeded: " << command << '\n';
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

static int require_contains(const std::string& path, const std::string& needle) {
    const auto text = read_file(path);
    if (text.find(needle) == std::string::npos) {
        std::cerr << "missing '" << needle << "' in " << path << '\n';
        return 1;
    }
    return 0;
}

static void cleanup_outputs() {
    std::remove("cli_output_attrs_smoke.gff3");
    std::remove("output_attrs_main.tsv");
    std::remove("output_attrs_verbose.tsv");
    std::remove("output_attrs_old.tsv");
    std::remove("output_attrs_json.json");
    std::remove("query_output_attrs_main.tsv");
    std::remove("query_output_attrs_verbose.tsv");
    std::remove("query_output_attrs_old.tsv");
    std::remove("output_attrs_help.txt");
    std::remove("query_output_attrs_help.txt");
    std::remove("output_attrs_regions.bed");
    std::remove("output_attrs_bad.out");
    std::remove("output_attrs_bad.err");
    std::remove("output_attrs_bad.gff3");
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: cli_output_attrs_smoke <gffsub-executable>\n";
        return 2;
    }

    const std::string exe = std::string{"\""} + argv[1] + "\"";
    const std::string gff{"cli_output_attrs_smoke.gff3"};
    if (!write_test_annotation(gff)) {
        std::cerr << "cannot write test annotation\n";
        return 1;
    }
    {
        std::ofstream bed{"output_attrs_regions.bed"};
        if (!bed.is_open()) {
            std::cerr << "cannot write BED test regions\n";
            return 1;
        }
        bed << "chr1\t99\t400\n";
    }

    const std::string keys{"ID,Name,Alias,Dbxref"};
    if (run_command(exe + " " + gff + " --help > output_attrs_help.txt 2>&1") != 0 ||
        require_contains("output_attrs_help.txt", "--out-attrs KEYS") != 0 ||
        require_contains("output_attrs_help.txt", "--output-attrs KEYS") != 0) {
        return 1;
    }
    if (run_command(exe + " query --help > query_output_attrs_help.txt 2>&1") != 0 ||
        require_contains("query_output_attrs_help.txt", "--out-attrs KEYS") != 0 ||
        require_contains("query_output_attrs_help.txt", "Output selected attributes") != 0 ||
        require_contains("query_output_attrs_help.txt", "--output-attrs KEYS") != 0 ||
        require_contains("query_output_attrs_help.txt", "Verbose alias for --out-attrs") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --out-attrs " + keys + " > output_attrs_main.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --output-attrs " + keys + " > output_attrs_verbose.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --attrs " + keys + " > output_attrs_old.tsv") != 0) {
        return 1;
    }
    if (compare_files("output_attrs_main.tsv", "output_attrs_verbose.tsv") != 0 ||
        compare_files("output_attrs_main.tsv", "output_attrs_old.tsv") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --summary json --out-attrs " + keys +
                    " > output_attrs_json.json") != 0 ||
        require_contains("output_attrs_json.json", "\"attrs\":{\"ID\":\"gene0001\",\"Name\":\"ABC1\"") != 0 ||
        require_contains("output_attrs_json.json", "\"Alias\":\"ABC-1\"") != 0 ||
        require_contains("output_attrs_json.json", "\"Dbxref\":\"GeneID:123\"") != 0) {
        return 1;
    }

    if (run_command(exe + " query " + gff + " --id gene0001 --out-attrs " + keys + " > query_output_attrs_main.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --output-attrs " + keys + " > query_output_attrs_verbose.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --attrs " + keys + " > query_output_attrs_old.tsv") != 0) {
        return 1;
    }
    if (compare_files("query_output_attrs_main.tsv", "query_output_attrs_verbose.tsv") != 0 ||
        compare_files("query_output_attrs_main.tsv", "query_output_attrs_old.tsv") != 0) {
        return 1;
    }

    const std::string bad_redirect{" > output_attrs_bad.out 2> output_attrs_bad.err"};
    if (expect_command_failure(exe + " " + gff + " --id gene0001 --out-attrs , " + bad_redirect) != 0 ||
        require_contains("output_attrs_bad.err", "Error: --out-attrs expects a comma-separated list of keys") != 0) {
        return 1;
    }

    if (expect_command_failure(exe + " " + gff + " --id gene0001 --out-attrs " + keys +
                               " --bed output_attrs_regions.bed" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --out-attrs " + keys +
                               " --longest" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --out-attrs " + keys +
                               " --threads 2" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --out-attrs " + keys +
                               " --output-format gtf3" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --out-attrs " + keys +
                               " --output output_attrs_bad.gff3" + bad_redirect) != 0) {
        return 1;
    }
    if (require_contains("output_attrs_bad.err", "Error: --summary/--summary-format/--out-attrs only supports query-style selectors") != 0) {
        return 1;
    }

    if (expect_command_failure(exe + " " + gff + " --id gene0001 --summary xml" + bad_redirect) != 0 ||
        require_contains("output_attrs_bad.err", "Error: --summary expects tsv or json") != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --summary-format xml" + bad_redirect) != 0 ||
        require_contains("output_attrs_bad.err", "Error: --summary-format expects tsv or json") != 0) {
        return 1;
    }

    if (expect_command_failure(exe + " " + gff + " --id gene0001 --summary tsv" +
                               " --bed output_attrs_regions.bed" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --summary tsv" +
                               " --longest" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --summary tsv" +
                               " --threads 2" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --summary tsv" +
                               " --output-format gtf3" + bad_redirect) != 0 ||
        expect_command_failure(exe + " " + gff + " --id gene0001 --summary tsv" +
                               " --output output_attrs_bad.gff3" + bad_redirect) != 0) {
        return 1;
    }

    cleanup_outputs();
    std::cout << "cli_output_attrs_smoke OK\n";
    return 0;
}
