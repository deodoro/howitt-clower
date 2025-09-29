#include <stdint.h>
#include <stdio.h>

typedef struct {
    uint64_t state;
    uint64_t inc;   // must be odd
} pcg32_t;

static uint32_t pcg32_next(pcg32_t *rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1ULL);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

static void pcg32_seed(pcg32_t *rng, uint64_t initstate, uint64_t initseq) {
    rng->state = 0U;
    rng->inc   = (initseq << 1u) | 1u;
    pcg32_next(rng);
    rng->state += initstate;
    pcg32_next(rng);
}

int main(void) {
    pcg32_t rng;
    pcg32_seed(&rng, 42ULL, 54ULL);
    for (int i = 0; i < 20; ++i) {
        printf("%u%s", pcg32_next(&rng), (i == 19 ? "\n" : " "));
    }
    return 0;
}
