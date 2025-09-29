class PCG32:
    def __init__(self, initstate: int, initseq: int):
        self.state = 0
        self.inc   = ((initseq & ((1<<64)-1)) << 1) | 1
        self._next()
        self.state = (self.state + (initstate & ((1<<64)-1))) & ((1<<64)-1)
        self._next()

    def _next(self) -> int:
        oldstate = self.state
        self.state = (oldstate * 6364136223846793005 + self.inc) & ((1<<64)-1)
        xorshifted = (((oldstate >> 18) ^ oldstate) >> 27) & 0xFFFFFFFF
        rot = (oldstate >> 59) & 31
        return ((xorshifted >> rot) | ((xorshifted << ((-rot) & 31)) & 0xFFFFFFFF)) & 0xFFFFFFFF

    def next_u32(self) -> int:
        return self._next()

if __name__ == "__main__":
    rng = PCG32(42, 54)
    nums = [rng.next_u32() for _ in range(20)]
    print(" ".join(str(x) for x in nums))
