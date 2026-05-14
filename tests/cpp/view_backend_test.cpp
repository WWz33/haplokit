#include <iostream>
#include <stdexcept>

#include "selector.h"
#include "vcf_reader.h"
#include "view_backend.h"

namespace {

int run() {
    haplokit::VcfReader reader("data/var.sorted.vcf.gz");

    const auto region_data = reader.fetch(haplokit::parse_region("scaffold_1:4300-5000"));
    haplokit::ViewOptions strict_region;
    strict_region.by = haplokit::GroupBy::Region;
    strict_region.output_mode = haplokit::OutputMode::Summary;
    const auto summary = haplokit::build_view_result(region_data, strict_region);
    if (summary.grouping_mode != "strict-region") {
        std::cerr << "unexpected grouping mode: " << summary.grouping_mode << "\n";
        return 1;
    }
    if (summary.variant_count != 5 || summary.sample_count != 37) {
        std::cerr << "unexpected counts: variants=" << summary.variant_count
                  << " samples=" << summary.sample_count << "\n";
        return 1;
    }
    if (summary.haplotype_count <= 0 || summary.haplotypes.empty()) {
        std::cerr << "expected non-empty summary haplotypes\n";
        return 1;
    }
    if (summary.sites.size() != 5U || summary.sites.front().pos != 4300) {
        std::cerr << "expected site metadata for summary plot\n";
        return 1;
    }

    haplokit::ViewOptions strict_site;
    strict_site.by = haplokit::GroupBy::Site;
    strict_site.output_mode = haplokit::OutputMode::Detail;
    const auto site_data = reader.fetch(haplokit::parse_region("scaffold_1:4300"));
    const auto detail = haplokit::build_view_result(site_data, strict_site);
    if (detail.grouping_mode != "strict-site") {
        std::cerr << "unexpected site grouping mode: " << detail.grouping_mode << "\n";
        return 1;
    }
    if (detail.variant_count != 1 || detail.accessions.empty()) {
        std::cerr << "expected one site variant and non-empty detail rows\n";
        return 1;
    }

    haplokit::ViewOptions approx_region;
    approx_region.by = haplokit::GroupBy::Region;
    approx_region.output_mode = haplokit::OutputMode::Summary;
    approx_region.max_diff = 0.2;
    const auto approx = haplokit::build_view_result(region_data, approx_region);
    if (approx.grouping_mode != "approx-region" || approx.grouping_method != "max-diff") {
        std::cerr << "unexpected approx metadata\n";
        return 1;
    }

    const auto payload = haplokit::serialize_view_result_json(summary);
    if (payload.find("\"grouping_mode\":\"strict-region\"") == std::string::npos) {
        std::cerr << "serialized payload missing grouping_mode\n";
        return 1;
    }
    if (payload.find("\"haplotypes\"") == std::string::npos) {
        std::cerr << "serialized payload missing haplotypes\n";
        return 1;
    }
    if (payload.find("\"sites\"") == std::string::npos) {
        std::cerr << "serialized payload missing sites\n";
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
