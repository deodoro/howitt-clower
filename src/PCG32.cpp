#include "PCG32.h"
#include <iostream>

#define DEBUG 0

void PCG32::seed(uint64_t initstate, uint64_t initseq) {
    state = 0U;
    inc = (initseq << 1u) | 1u;
    next();
    state += initstate;
    next();
}

uint32_t PCG32::next() {
    uint64_t oldstate = state;
    state = oldstate * 6364136223846793005ULL + (inc | 1ULL);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

int PCG32::uniform_int(int n) {
    uint32_t res = next();
    if (DEBUG) {
        std::cout << "RND: " << res << std::endl;
    }
    return (int)(res % (uint32_t)n);
}

double PCG32::uniform01_inclusive() {
    // Match original granularity of 0.001 if desired; here use full precision
    // Return (0,1] by excluding exact 0.0
    int x = uniform_int(1000);
    return double(1 + x)/1000.0;
}
