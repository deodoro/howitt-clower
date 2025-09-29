#ifndef PCG32_H
#define PCG32_H

#include <cstdint>

class PCG32 {
public:
    uint64_t state{0};
    uint64_t inc{0}; // must be odd

    void seed(uint64_t initstate, uint64_t initseq);
    uint32_t next();
    int uniform_int(int n);
    double uniform01_inclusive();
};

#endif // PCG32_H
