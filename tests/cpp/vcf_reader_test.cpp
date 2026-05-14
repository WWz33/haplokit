#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "selector.h"
#include "vcf_reader.h"

namespace {

std::filesystem::path repo_root() {
    return std::filesystem::current_path();
}

int run() {
    const auto vcf = repo_root() / "data" / "var.sorted.vcf.gz";
    if (!std::filesystem::exists(vcf)) {
        throw std::runtime_error("missing test fixture: " + vcf.string());
    }

    haptools::VcfReader reader(vcf.string());

    const auto region_data = reader.fetch(haptools::parse_region("scaffold_1:4300-5000"));
    if (region_data.samples.size() != 37U) {
        std::cerr << "expected 37 samples, got " << region_data.samples.size() << "\n";
        return 1;
    }
    if (region_data.variants.size() != 5U) {
        std::cerr << "expected 5 region variants, got " << region_data.variants.size() << "\n";
        return 1;
    }
    if (region_data.variants.front().pos != 4300 || region_data.variants.back().pos != 4950) {
        std::cerr << "unexpected region bounds " << region_data.variants.front().pos << " -> "
                  << region_data.variants.back().pos
                  << "\n";
        return 1;
    }
    if (region_data.variants.front().genotypes.size() != region_data.samples.size()) {
        std::cerr << "expected one genotype call per sample on first variant\n";
        return 1;
    }

    const auto site_data = reader.fetch(haptools::parse_region("scaffold_1:4300"));
    if (site_data.variants.size() != 1U) {
        std::cerr << "expected 1 site variant, got " << site_data.variants.size() << "\n";
        return 1;
    }
    if (site_data.variants.front().pos != 4300) {
        std::cerr << "unexpected site position " << site_data.variants.front().pos << "\n";
        return 1;
    }

    const auto subset_data = reader.fetch(haptools::parse_region("scaffold_1:4300-5000"), {"C1", "C5", "C16"});
    if (subset_data.samples.size() != 3U) {
        std::cerr << "expected 3 subset samples, got " << subset_data.samples.size() << "\n";
        return 1;
    }
    if (subset_data.variants.front().genotypes.size() != subset_data.samples.size()) {
        std::cerr << "expected subset genotype width to match requested samples\n";
        return 1;
    }

    return 0;
}

}  // namespace

int main() {
    try {
        return run();
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return 2;
    }
}
