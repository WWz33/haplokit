#include "hamming_simd.h"
#include <algorithm>
#include <cstring>

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace haplokit {
namespace network {

bool HammingSIMD::avx2_checked_ = false;
bool HammingSIMD::avx2_available_ = false;

bool HammingSIMD::has_avx2_support() {
    if (avx2_checked_) return avx2_available_;

#ifdef __AVX2__
    // Runtime CPU feature detection
    unsigned int eax, ebx, ecx, edx;
    __asm__ __volatile__(
        "cpuid"
        : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
        : "a"(7), "c"(0)
    );
    avx2_available_ = (ebx & (1 << 5)) != 0; // AVX2 bit
#else
    avx2_available_ = false;
#endif

    avx2_checked_ = true;
    return avx2_available_;
}

std::vector<uint32_t> HammingSIMD::pack_sequence(const std::string& seq) {
    size_t len = seq.size();
    size_t packed_len = (len + 15) / 16; // 16 bases per uint32_t
    std::vector<uint32_t> packed(packed_len, 0);

    for (size_t i = 0; i < len; ++i) {
        uint8_t encoded = encode_base(seq[i]);
        size_t word_idx = i / 16;
        size_t bit_pos = (i % 16) * 2;
        packed[word_idx] |= (static_cast<uint32_t>(encoded) << bit_pos);
    }

    return packed;
}

int HammingSIMD::hamming_distance_scalar(const std::string& seq1, const std::string& seq2) {
    size_t len = std::min(seq1.size(), seq2.size());
    int dist = 0;

    for (size_t i = 0; i < len; ++i) {
        if (seq1[i] != seq2[i]) {
            ++dist;
        }
    }

    // Add length difference as mismatches
    dist += std::abs(static_cast<int>(seq1.size()) - static_cast<int>(seq2.size()));

    return dist;
}

#ifdef __AVX2__
int HammingSIMD::hamming_distance_simd(const uint32_t* packed1, const uint32_t* packed2, size_t packed_len) {
    int dist = 0;
    size_t i = 0;

    // Process 8 uint32_t at a time (256 bits = 128 bases)
    for (; i + 8 <= packed_len; i += 8) {
        __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(packed1 + i));
        __m256i v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(packed2 + i));

        // XOR to find differing bits
        __m256i xor_result = _mm256_xor_si256(v1, v2);

        // Count set bits (each 2-bit difference = 1 mismatch)
        // Extract to scalar and count
        uint32_t temp[8];
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(temp), xor_result);

        for (int j = 0; j < 8; ++j) {
            uint32_t val = temp[j];
            // Count 2-bit groups that are non-zero
            for (int k = 0; k < 16; ++k) {
                if ((val & 3) != 0) ++dist;
                val >>= 2;
            }
        }
    }

    // Handle remaining elements
    for (; i < packed_len; ++i) {
        uint32_t xor_val = packed1[i] ^ packed2[i];
        for (int k = 0; k < 16; ++k) {
            if ((xor_val & 3) != 0) ++dist;
            xor_val >>= 2;
        }
    }

    return dist;
}
#else
int HammingSIMD::hamming_distance_simd(const uint32_t* packed1, const uint32_t* packed2, size_t packed_len) {
    // Fallback to scalar when AVX2 not available
    int dist = 0;
    for (size_t i = 0; i < packed_len; ++i) {
        uint32_t xor_val = packed1[i] ^ packed2[i];
        for (int k = 0; k < 16; ++k) {
            if ((xor_val & 3) != 0) ++dist;
            xor_val >>= 2;
        }
    }
    return dist;
}
#endif

int HammingSIMD::hamming_distance(const std::string& seq1, const std::string& seq2) {
    // Use SIMD if available and sequences are long enough
    if (has_avx2_support() && seq1.size() >= 128 && seq1.size() == seq2.size()) {
        auto packed1 = pack_sequence(seq1);
        auto packed2 = pack_sequence(seq2);
        return hamming_distance_simd(packed1.data(), packed2.data(), packed1.size());
    }

    return hamming_distance_scalar(seq1, seq2);
}

} // namespace network
} // namespace haplokit
