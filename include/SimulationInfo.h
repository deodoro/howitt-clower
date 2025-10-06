/*
 * SimulationInfo.h
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

#ifndef SIMULATIONINFO_H
#define SIMULATIONINFO_H

#include <string>
#include <sstream>

/**
 * @brief Class to store simulation level parameters.
 */
class SimulationInfo {
public:
    int n;        // number of goods
    int bsize;    // each (i!=j) type count
    int K;        // shop locations
    double f1;    // fixed overhead cost
    int xMax;     // maximum target value
    double lambda; // matching probability
    double alpha;  // target update rate
    double theta;  // exit probability
    double C;      // variable cost
    int persist;   // persistence for equilibrium
    int RANDSEED; // random seed
    int m;        // number of traders
    int numruns;  // number of runs
    int FirstSlope; // first slope value
    int LastSlope;  // last slope value
    int T;        // weeks

    /**
     * @brief Default constructor with default values.
     */
    SimulationInfo()
        : n(10), bsize(24), K(200), f1(0.0), xMax(200), lambda(0.05),
        // : n(3), bsize(2), K(6), f1(0.0), xMax(200), lambda(0.05),
          alpha(0.25), theta(0.01), C(5.0), persist(10), RANDSEED(1),
          m(bsize * n * (n - 1)), numruns(1), FirstSlope(16), LastSlope(18),
          T(20000) {}

    /**
     * @brief Convert the simulation info to a JSON string.
     * @return A JSON string representation of the simulation parameters.
     */
    std::string to_json() const;
};

#endif // SIMULATIONINFO_H
