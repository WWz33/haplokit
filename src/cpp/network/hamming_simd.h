#ifndef HAPLOKIT_HAMMING_SIMD_H
#define HAPLOKIT_HAMMING_SIMD_H

#include <cstdint>
#include <string>
#include <vector>

namespace haplokit {
namespace network {

// SIMD-accelerated Hamming distance calculation
// Packs 16 nucleotides into 32-bit integers (2-bit encoding: A=00, C=01, G=10, T=11)
// Processes 8 integers in parallel using AVX2 (256-bit registers)

class HammingSIMD {
public:
    // Encode nucleotide to 2-bit representation
    static inline uint8_t encode_base(char base) {
        switch (base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0; // treat N/gap as A for distance calculation
        }
    }

    // Pack sequence into 32-bit integers (16 bases per int)
    static std::vector<uint32_t> pack_sequence(const std::string& seq);

    // Compute Hamming distance between two sequences
    // Returns number of differing positions
    static int hamming_distance(const std::string& seq1, const std::string& seq2);

    // SIMD-optimized version (requires AVX2)
    static int hamming_distance_simd(const uint32_t* packed1, const uint32_t* packed2, size_t packed_len);

    // Fallback scalar version
    static int hamming_distance_scalar(const std::string& seq1, const std::string& seq2);

    // Check if CPU supports AVX2
    static bool has_avx2_support();

private:
    static bool avx2_checked_;
    static bool avx2_available_;
};

} // namespace network
} // namespace haplokit

#endif // HAPLOKIT_HAMMING_SIMD_H
