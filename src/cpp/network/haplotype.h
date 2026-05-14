#ifndef HAPLOKIT_NETWORK_HAPLOTYPE_H
#define HAPLOKIT_NETWORK_HAPLOTYPE_H

#include <string>
#include <vector>

namespace haplokit {
namespace network {

// Common haplotype structure used by all network algorithms
struct Haplotype {
    std::string sequence;
    std::vector<std::string> samples;
};

} // namespace network
} // namespace haplokit

#endif // HAPLOKIT_NETWORK_HAPLOTYPE_H
