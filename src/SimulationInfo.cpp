/*
 * SimulationInfo.cpp
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

#include "SimulationInfo.h"

std::string SimulationInfo::to_json() const {
    std::ostringstream oss;
    oss << "{";
    oss << "\"n\":" << n << ",";
    oss << "\"bsize\":" << bsize << ",";
    oss << "\"K\":" << K << ",";
    oss << "\"f1\":" << f1 << ",";
    oss << "\"xMax\":" << xMax << ",";
    oss << "\"lambda\":" << lambda << ",";
    oss << "\"alpha\":" << alpha << ",";
    oss << "\"theta\":" << theta << ",";
    oss << "\"C\":" << C << ",";
    oss << "\"persist\":" << persist << ",";
    oss << "\"RANDSEED\":" << RANDSEED << ",";
    oss << "\"m\":" << m << ",";
    oss << "\"numruns\":" << numruns << ",";
    oss << "\"FirstSlope\":" << FirstSlope << ",";
    oss << "\"LastSlope\":" << LastSlope << ",";
    oss << "\"T\":" << T;
    oss << "}";
    return oss.str();
}
