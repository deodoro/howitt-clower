/*
 * PCG32.h
 *
 * Author: Jose Deodoro <deodoro.filho@gmail.com> <jdeoliv@gmu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
