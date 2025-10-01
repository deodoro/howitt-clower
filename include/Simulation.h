/*
 * Simulation.h
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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "PCG32.h"
#include "Shop.h"
#include "Trader.h"
#include "ResearchResults.h"
#include <vector>
#include <array>
#include <ctime>
#include <cstdio>
#include <functional>

struct MatchEvaluation {
    double Ucomp = 0.0;
    double Ubarter = 0.0;

    Shop *barter = nullptr;
    Shop *candidate_0 = nullptr;
    Shop *candidate_1 = nullptr;
};

class Simulation {
public:
    // Configuration constants (mirror attached code)
    static constexpr int numruns   = 1;
    static constexpr int FirstSlope = 16;
    static constexpr int LastSlope  = 18;
    static constexpr double f1     = 0.0;

    static constexpr int n      = 10;   // goods
    static constexpr int bsize  = 24;   // each (i!=j) type count
    static constexpr int K      = 200;  // shop locations

    static constexpr int xMax   = 200;
    static constexpr double lambda = 0.05;
    static constexpr double alpha  = 0.25;
    static constexpr double theta  = 0.01;
    static constexpr double C      = 5.0;
    static constexpr int persist   = 10;

    static constexpr int T       = 20000; // weeks
    static constexpr int RANDSEED = 1;
    static constexpr int RptPer   = 1;
    static constexpr int prtoscr  = 1;

    static constexpr int m = bsize * n * (n - 1); // number of traders

    // Constructor builds the static economy (agents and lists)
    Simulation();

    void run_all();

private:
    // Derived/utility
    static constexpr int PRINT_LOOP_N = 6;

    double slope{16.0};

    // State
    std::vector<Trader> traders;
    std::vector<Shop> shops;
    std::vector<Trader*> trader_line; // Weekly lineup - permutation of traders

    std::vector<std::vector<int>> produces; // by good i: trader ids producing i
    std::vector<std::vector<int>> consumes; // by good j: trader ids desiring j

    std::vector<double> usingmoney; // per good
    std::vector<double> Pinv; // SME inverse retail per good
    std::vector<double> vol[2], avp[2];

    // run-scoped variables
    int NS{0};        // active shops
    int BS{0};        // active non-money shops
    int devyear{-1};
    int monyear{-1};
    int endcount{0}, devcount{0};
    int monetary{0}, fulldev{0}, moneygood{0};
    double Fmon{0.0};
    double W{0.0}, SurpSME{0.0};
    double part{0.0}, moneytraders{0.0}, usingmax{0.0};
    double Csurp{0.0}, Psurp{0.0};
    double R[2]{-1.0, -1.0};

    // time & IO
    std::time_t firstbegin{}, finish{};
    std::clock_t clock_begin{}, clock_finish{};
    std::FILE* stream{nullptr};

    // run and time variables
    int run{1};
    int t{0};
    int Nshop{0};

    // Helper methods
    std::function<double(int)> make_overhead(double f1, double slope);
    void init_static();
    void init_run();
    void lineup();
    void weekly_entry();
    void weekly_matching();
    void report_trader(Trader const& trader);
    void weekly_trade_and_exit();
    void weekly_update_prices();
    ResearchResults research(Trader& trader);
    Trader &random_consumer(int good);
    Trader &random_producer(int good);
    Shop &random_shop();
    void addshop(const Trader *trader, Shop *shop, std::vector<int> &cand);
    void try_barter(Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval);
    void try_one(const Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval);
    void try_two(const Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval);
    int calc1();
    void calc2();
    void rmse();
    void report(int tt);
    void print_debug(std::string title) const;
    void print_traders() const;
    void print_shops() const;
};

#endif // SIMULATION_H
