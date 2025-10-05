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
#include <string>
#include <sstream>
#include <vector>

/**
 * @brief Class to store run-scoped variables for the simulation.
 *
 * This encapsulates the state variables that are specific to each simulation run,
 * including statistics, timing, and IO handles.
 */
class RunInfo {
public:
    RunInfo(int n) : info_n(n) {}

    /**
     * @brief Create a copy of this RunInfo instance.
     * @return A pointer to a new RunInfo object with the same state.
     */
    RunInfo* clone() const {
        RunInfo* copy = new RunInfo(info_n);
        copy->NumberOfShops = NumberOfShops;
        copy->BS = BS;
        copy->Nshop = Nshop;
        copy->Slope = Slope;
        copy->info_n = info_n;
        copy->devyear = devyear;
        copy->monyear = monyear;
        copy->endcount = endcount;
        copy->devcount = devcount;
        copy->monetary = monetary;
        copy->fulldev = fulldev;
        copy->moneygood = moneygood;
        copy->Fmon = Fmon;
        copy->W = W;
        copy->SurpSME = SurpSME;
        copy->part = part;
        copy->moneytraders = moneytraders;
        copy->usingmax = usingmax;
        copy->Csurp = Csurp;
        copy->Psurp = Psurp;
        copy->R[0] = R[0];
        copy->R[1] = R[1];
        copy->time_spent = time_spent;
        copy->usingmoney = usingmoney;
        copy->run = run;
        copy->t = t;
        return copy;
    }

    // Active shops and related counts
    int NumberOfShops{0};        // active shops
    int BS{0};                   // active non-money shops
    int Nshop{0};                // additional shop count
    int Slope{0};
    int info_n{0};

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
    double time_spent{0.0};

    // Usage of money goods
    std::vector<double> usingmoney;

    // Run and time variables
    int run{1};
    int t{0};

    void report() {
        // Reports simulation progress and statistics at specified intervals.

        if (t == -1)
            std::printf("***");
        std::printf("%6.0f %6.0f ", part, moneytraders);
        for (int b = 1; b <= 5 && b <= info_n; ++b) {
            std::printf("%6.0f ", usingmoney[b]);
        }
        std::printf("%6d %4d\n", (t == -1 ? t : t) / 50, NumberOfShops);
    }

    /**
     * @brief Convert the run info to a JSON string.
     * @return A JSON string representation of the run parameters.
     */
    std::string to_json() const {
        std::ostringstream oss;
        oss << "{";
        oss << "\"NumberOfShops\":" << NumberOfShops << ",";
        oss << "\"BS\":" << BS << ",";
        oss << "\"Nshop\":" << Nshop << ",";
        oss << "\"Slope\":" << Slope << ",";
        oss << "\"devyear\":" << devyear << ",";
        oss << "\"monyear\":" << monyear << ",";
        oss << "\"endcount\":" << endcount << ",";
        oss << "\"devcount\":" << devcount << ",";
        oss << "\"monetary\":" << monetary << ",";
        oss << "\"fulldev\":" << fulldev << ",";
        oss << "\"moneygood\":" << moneygood << ",";
        oss << "\"Fmon\":" << Fmon << ",";
        oss << "\"W\":" << W << ",";
        oss << "\"SurpSME\":" << SurpSME << ",";
        oss << "\"part\":" << part << ",";
        oss << "\"moneytraders\":" << moneytraders << ",";
        oss << "\"usingmax\":" << usingmax << ",";
        oss << "\"Csurp\":" << Csurp << ",";
        oss << "\"Psurp\":" << Psurp << ",";
        oss << "\"R\":[" << R[0] << "," << R[1] << "],";
        oss << "\"run\":" << run << ",";
        oss << "\"t\":" << t << ",";
        oss << "\"time_spent\":" << time_spent << ",";
        oss << "\"usingmoney\":[";
        for (size_t i = 0; i < usingmoney.size(); ++i) {
            oss << usingmoney[i];
            if (i < usingmoney.size() - 1) oss << ",";
        }
        oss << "]";
        oss << "}";
        return oss.str();
    }
};

#endif // RUNINFO_H
