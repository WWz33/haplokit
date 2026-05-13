#pragma once

#include <string>
#include <vector>

#include "selector.h"

namespace haptools {

struct GenotypeCall {
    int allele1 = -1;
    int allele2 = -1;
    bool phased = false;
};

struct VariantRecord {
    std::string chrom;
    int pos = 0;
    std::string ref;
    std::vector<std::string> alt;
    std::vector<GenotypeCall> genotypes;
};

struct RegionData {
    std::vector<std::string> samples;
    std::vector<VariantRecord> variants;
};

class VcfReader {
public:
    explicit VcfReader(std::string path);

    RegionData fetch(const Region& region, const std::vector<std::string>& samples = {}) const;

private:
    std::string path_;
};

}  // namespace haptools
