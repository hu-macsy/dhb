#pragma once

#include <cstdint>

namespace dhb {

// Fast modulo reduction
// https://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
inline uint32_t reduce(uint32_t x, uint32_t operand) {
    return uint32_t(uint64_t(x) * uint64_t(operand) >> 32);
}

inline uint32_t hash32(uint32_t x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

// Thanks to David Stafford for his further research on Austin Appleby's
// MurmurHash3 specifically for input values with low entropy which are more
// close to counting numbers than random ones.
//
// The one we are using here is Mix13.
//
// Please go to
// http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html for
// further reading.
//
// Also thanks to Thomas Mueller for his answer on good integer hashing:
// https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
//
// The final function definition comes from Sebastiano Vigna:
// https://xorshift.di.unimi.it/splitmix64.c
inline uint64_t hash64(uint64_t x) {
    constexpr uint64_t op_a = 0xbf58476d1ce4e5b9;
    constexpr uint64_t op_b = 0x94d049bb133111eb;

    x = (x ^ (x >> 30)) * op_a;
    x = (x ^ (x >> 27)) * op_b;
    x = x ^ (x >> 31);
    return x;
}

} // namespace dhb