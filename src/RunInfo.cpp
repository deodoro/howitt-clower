/*
 * RunInfo.cpp
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

#include "RunInfo.h"

RunInfo::RunInfo(int n) : info_n(n) {}

RunInfo* RunInfo::clone() const {
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

void RunInfo::report() {
    // Reports simulation progress and statistics at specified intervals.

    if (t == -1)
        std::printf("***");
    std::printf("%6.0f %6.0f ", part, moneytraders);
    for (int b = 1; b <= 5 && b <= info_n; ++b) {
        std::printf("%6.0f ", usingmoney[b]);
    }
    std::printf("%6d %4d\n", (t == -1 ? t : t) / 50, NumberOfShops);
}

std::string RunInfo::to_json() const {
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
    for (size_t i = 1; i < usingmoney.size(); ++i) {
        oss << usingmoney[i];
        if (i < usingmoney.size() - 1) oss << ",";
    }
    oss << "]";
    oss << "}";
    return oss.str();
}
