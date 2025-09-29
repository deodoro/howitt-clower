// HorwittClower_OOP.cpp
// Object-oriented C++17 refactor of the Howitt–Clower simulation described in the attached source,
// preserving algorithmic behavior while encapsulating state and utilities into classes.
//
// Notes:
// - Indices remain 1-based where the original program relied on that convention.
// - Constants mirror the attached file; pricing F(tr0,tr1,f_other) and overhead f(i) follow the same rules.
// - RNG uses a PCG32 class mirroring the original functions.
// - "Shop" and "Trader" are explicit types; simulation state is isolated in Simulation.
// - File output to evol.fil and periodic reporting preserved (guarded by prtoscr).
//
// Build: g++ -std=gnu++17 -O2 -pipe -o howclow_oop main.cpp PCG32.cpp Shop.cpp Trader.cpp Simulation.cpp
// Run:   ./howclow_oop
//
// Based on: horwitt-clower.c
// --------------------------------------------------------------------

#include "Simulation.h"

// -------------------------
// Entry point
// -------------------------
int main() {
    Simulation sim;
    sim.run_all();
    return 0;
}
