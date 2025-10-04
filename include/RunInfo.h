/*
 * RunInfo.h
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

#ifndef RUNINFO_H
#define RUNINFO_H

#include <ctime>
#include <cstdio>

/**
 * @brief Class to store run-scoped variables for the simulation.
 *
 * This encapsulates the state variables that are specific to each simulation run,
 * including statistics, timing, and IO handles.
 */
class RunInfo {
public:
    // Active shops and related counts
    int NumberOfShops{0};        // active shops
    int BS{0};                   // active non-money shops
    int Nshop{0};                // additional shop count
    int Slope{0};

    // Development and monetary emergence tracking
    int devyear{-1};
    int monyear{-1};
    int endcount{0}, devcount{0};
    int monetary{0}, fulldev{0}, moneygood{0};

    // Economic aggregates
    double Fmon{0.0};
    double W{0.0}, SurpSME{0.0};
    double part{0.0}, moneytraders{0.0}, usingmax{0.0};
    double Csurp{0.0}, Psurp{0.0};
    double R[2]{-1.0, -1.0};

    // Run and time variables
    int run{1};
    int t{0};
};

#endif // RUNINFO_H
