#include "gff3.hpp"
#include <fstream>

namespace gffsub {

static std::vector<std::string> split_line(const std::string& line, char delimiter) {
    std::vector<std::string> cols;
    size_t start = 0;
    while (true) {
        auto pos = line.find(delimiter, start);
        if (pos == std::string::npos) {
            cols.emplace_back(line.substr(start));
            break;
        }
        cols.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return cols;
}

void filter_by_region(GffData& data, const Region& region) {
    for (auto& rec : data) {
        if (rec.seqid != region.seqid || rec.end < region.start || rec.start > region.end) {
            rec.kept = false;
        }
    }
}

static std::vector<Region> load_regions(const std::string& filename, bool is_bed) {
    std::vector<Region> regions;
    if (is_bed) {
        std::ifstream file(filename);
        if (!file.is_open()) return regions;

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) continue;

            Region r;
            try {
                r.seqid = cols[0];
                r.start = std::stoll(cols[1]) + 1;
                r.end = std::stoll(cols[2]);
            } catch (const std::exception&) {
                continue;
            }
            regions.push_back(r);
        }
    }
    return regions;
}

void filter_by_regions_from_file(GffData& data, const std::string& bed_file) {
    auto regions = load_regions(bed_file, true);
    for (auto& rec : data) {
        bool in_region = false;
        for (const auto& r : regions) {
            if (rec.seqid == r.seqid && rec.end >= r.start && rec.start <= r.end) {
                in_region = true;
                break;
            }
        }
        if (!in_region) rec.kept = false;
    }
}

void filter_by_feature(GffData& data, std::string_view feature_type) {
    for (auto& rec : data) {
        if (rec.kept && rec.type != feature_type) {
            rec.kept = false;
        }
    }
}

void filter_by_seqid(GffData& data, std::string_view seqid) {
    for (auto& rec : data) {
        if (rec.kept && rec.seqid != seqid) {
            rec.kept = false;
        }
    }
}

void filter_by_source(GffData& data, std::string_view source) {
    for (auto& rec : data) {
        if (rec.kept && rec.source != source) {
            rec.kept = false;
        }
    }
}

void filter_by_score(GffData& data, std::optional<double> score) {
    for (auto& rec : data) {
        if (rec.kept && rec.score != score) {
            rec.kept = false;
        }
    }
}

void filter_by_strand(GffData& data, char strand) {
    for (auto& rec : data) {
        if (rec.kept && rec.strand != strand) {
            rec.kept = false;
        }
    }
}

void filter_by_phase(GffData& data, char phase) {
    for (auto& rec : data) {
        if (rec.kept && rec.phase != phase) {
            rec.kept = false;
        }
    }
}

}  // namespace gffsub
