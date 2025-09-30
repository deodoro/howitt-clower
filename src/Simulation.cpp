/*
 * Simulation.cpp
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

#include "Simulation.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <stack>
#include <array>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <functional>
#include <assert.h>

#define DEBUG 0

// Global RNG instance
PCG32 rng;

Simulation::Simulation() {
    traders.resize(m + 1); // 1-based
    shops.resize(K + 1);   // 1-based
    for (size_t i = 1; i < traders.size(); ++i) {
        traders[i].set_shops(&shops);
    }
    produces.assign(n + 1, {});
    consumes.assign(n + 1, {});
    usingmoney.assign(n + 1, 0.0);
    Pinv.assign(n + 1, 0.0);
    for (int h = 0; h < 2; ++h) {
        vol[h].assign(n + 1, 0.0);
        avp[h].assign(n + 1, 0.0);
    }
    init_static();
}

void Simulation::run_all() {
    // Record the start time of the entire simulation
    time(&firstbegin);
    
    // Initialize the output file "evol.fil" with headers for evolutionary data
    stream = std::fopen("evol.fil", "w+");
    if (stream) {
        std::fprintf(stream, " run slope dev part mon mtrd usmx Mgd Csur Psur Esur Nsh NS BS RMS_0 RMS_1 dyr myr\n");
        std::fclose(stream);
    }

    // Loop over different slope values for parameter sweeps
    for (int _s = FirstSlope; _s <= LastSlope; _s += 2) {
        slope = _s;
        
        // For each slope, run multiple simulation runs
        for (run = 1; run <= numruns; ++run) {
            // Initialize the run-specific state
            init_run();
            clock_begin = std::clock();

            // Main simulation loop: run for T time steps
            for (t = 1; t <= T; ++t) {
                // Ensure at least one shop exists before proceeding
                do {
                    weekly_entry();
                } while (NS == 0);
                
                // Perform weekly activities: matching, trading, exit, and price updates
                weekly_matching();
                weekly_trade_and_exit();
                weekly_update_prices();

                // Periodically report progress and check for monetary equilibrium
                if (t % (PRINT_LOOP_N * RptPer) == 0) {
                    report(t);
                    monetary = calc1();
                    if (monetary == 1) {
                        break;  // Exit early if monetary equilibrium is reached
                    }
                }
            }

            // Calculate final statistics for the run
            calc2();

            // Log timing information
            clock_finish = std::clock();
            std::printf("Run number %d. Time elapsed: %.2f seconds.\n",
                        run,
                        (clock_finish - clock_begin) / (double)CLOCKS_PER_SEC);
            std::printf("Slope equals %-.0f, xMax equals %d\n\n", slope, xMax);

            // Append run results to the output file
            stream = std::fopen("evol.fil", "a");
            if (stream) {
                std::fprintf(stream,
                    "%5d %5.0f %3d %4.0f %3d %4.0f %4.0f %3d %4.0f %6.0f %5.0f %3d %3d %3d %6.3f %6.3f %4d %4d\n",
                    run, slope, fulldev, part, monetary, moneytraders, usingmax,
                    moneygood, Csurp, Psurp, SurpSME, Nshop, NS, BS, R[0], R[1], devyear, monyear);
                std::fclose(stream);
            }

            // Add a blank line between runs in the output file
            stream = std::fopen("evol.fil", "a");
            if (stream) {
                std::fprintf(stream, "\n");
                std::fclose(stream);
            }
        }

        // Add blank lines between slope sweeps in the output file
        stream = std::fopen("evol.fil", "a");
        if (stream) {
            std::fprintf(stream, "\n\n n bsiz K xMax lambda alpha theta C Nrun Pst(yr) T RND\n");
            std::fprintf(stream, "%3d %4d %3d %5d %6.3f %6.3f %6.3f %2.0f %5d %3d %6d %3d\n",
                         n, bsize, K, xMax, lambda, alpha, theta, C, numruns, persist, T, RANDSEED);
            time(&finish);
            std::fprintf(stream, "\n Started : %s", std::ctime(&firstbegin));
            std::fprintf(stream, " Finished: %s\n", std::ctime(&finish));
            std::fclose(stream);
        }
    }
}

std::function<double(int)> Simulation::make_overhead(double f1, double slope) {
    return [f1, slope](int i) -> double {
        return f1 + (i - 1) * slope;
    };
}

void Simulation::init_static() {
    // Build all (i,j) types with multiplicity bsize, i != j
    // Traders are 1..m enumerating all pairs with bsize copies
    int r = 0;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (i == j) continue;
            for (int k = 1; k <= bsize; ++k) {
                ++r;
                traders[r].idx = r;  // PROVISIONAL: to make the code compatible while it's  refactored
                traders[r].set_supplies(i);
                traders[r].set_demands(j);
                produces[i].push_back(r);
                consumes[j].push_back(r);
            }
        }
    }
    // Theoretical max number of "money traders": bsize*(n-2)*(n-1)
    Fmon = bsize * (n - 2.0) * (n - 1.0);
}

void Simulation::init_run() {
    std::printf("Number Using Using Using Using Using Using \n");
    std::printf("Active Money good1 good2 good3 good4 good5 Year NS \n");

    endcount = 0; fulldev = 0; devyear = -1; monyear = -1; devcount = 0;

    rng.seed(RANDSEED + run - 1, RANDSEED + run - 1);

    clock_begin = std::clock();
    t = 0;
    NS = 0;

    for (Trader& trader : traders) {
        trader.set_seller_idx(0);
        trader.set_buyer_idx(0);
        trader.set_familyshop(0);
    }
    int i = 0;
    for (Shop& shop : shops) {
        shop.clear();
        shop.idx = i++;
    }

    lineup();
}

// Entry process: potential entrepreneur tries to open a shop
void Simulation::weekly_entry() {
    auto overhead_f = make_overhead(f1, slope);
    int r = rng.uniform_int(m) + 1; // prospective owner
    Trader&  trader = traders[r];
    if (NS < K && trader.get_familyshop() == 0) {
        ResearchResults ok = research(r);
        if (ok.enter > 0) {
            for (Shop& shop : shops) {
                if (!shop.active && (&shop != &shops.front())) { // skip index 0
                    NS++;
                    trader.open_shop(shop);
                    shop.set_targets(ok.targ0, ok.targ1);
                    shop.update_prices(C, overhead_f);
                    break;
                }
            }
        }
    }
    print_debug("Weekly entry completed");
}

// Matching: agents sample a small set of shops and adopt best links
void Simulation::weekly_matching() {
    for (Trader* trader_p : trader_line) {
        Trader& trader = *trader_p;
        double U = trader.utility(shops);
        double psearch = (U > 0.0 ? lambda : 1.0);
        // Skip condition: random or already owns a shop
        if (rng.uniform01_inclusive() < psearch && trader.get_family_shop() == nullptr) {
            struct MatchEvaluation eval;
            eval.Ucomp = U;
            eval.candidate_0 = trader.get_seller_shop();
            eval.candidate_1 = trader.get_buyer_shop();
            // candidate initialization with current links
            std::stack<Shop*> cand_stack;
            cand_stack.push(trader.get_seller_shop());
            cand_stack.push(trader.get_buyer_shop());
            
            std::vector<int> cand;
            cand.reserve(8);
            cand.push_back(trader.get_seller_idx()); // c[0]
            cand.push_back(trader.get_buyer_idx()); // c[1]

            // add friend outlets/sources and one random shop
            // Trader& comrade_ = traders[comrade(trader)];
            Trader& comrade_ = traders[trader.trade_comrade(produces)];
            addshop(&comrade_, comrade_.get_seller_shop(), cand);

            Trader& soulmate_ = traders[trader.soulmate(consumes)];
            addshop(&soulmate_, soulmate_.get_buyer_shop(), cand);

            addshop(&trader, &random_shop(), cand);

            if (cand.size() > 2) {
                try_barter(trader, cand, eval);
                if ((shops[cand[0]].get_good(trader.get_supplies()) != trader.get_demands() || shops[cand[0]].get_price_supply(trader.get_supplies()) == 0.0)) {
                    try_one(trader, cand, eval);
                }
                try_two(trader, cand, eval);

                if (eval.Ucomp < eval.Ubarter && eval.barter != nullptr) {
                    trader.set_buyer_shop(nullptr);
                    trader.set_seller_shop(eval.barter);
                } else {
                    // adopt c[0], c[1] as improved chain if any
                    trader.set_seller_shop(eval.candidate_0);
                    trader.set_buyer_shop(eval.candidate_1);
                }
            }
        }
        report_trader(trader);
    }
}

void Simulation::report_trader(Trader const& trader) {
    if (DEBUG) {
        std::cout << trader.to_string() << std::endl;
    }
}

// Trade/accounting and stochastic exit of unprofitable shops
void Simulation::weekly_trade_and_exit() {
    print_debug("Weekly trade begins");

    // reset shop weekly incomes
    for (Shop& shop : shops) {
        shop.reset_weekly_incomes();
    }
    // tally incomes from adopted relationships
    for (Trader& trader : traders) {
        Shop *shop_a = trader.get_seller_shop();
        if (shop_a) {
            if (shop_a->get_good(trader.get_supplies()) == trader.get_demands()) {
                // direct barter
                shop_a->add_income(trader.get_supplies(), 1.0, true);
            } else {
                Shop* shop_b = trader.get_buyer_shop();
                if (shop_b && shop_a->get_good(trader.get_supplies()) == shop_b->get_good(trader.get_demands())) {
                        // indirect via common intermediary
                        shop_a->add_income(trader.get_supplies(), 1.0, true);
                        shop_b->add_income(trader.get_demands(), shop_a->get_price_supply(trader.get_supplies()));
                }
            }
        }
    }
    // exit if unprofitable
    for (Shop& shop: shops) {
        if (shop.active) {
            if ((!shop.is_profitable(make_overhead(f1, slope))) && rng.uniform01_inclusive() <= theta) {
                // sever links
                for (Trader& trader: traders) { trader.sever_links(shop); }
                shop.clear();
                NS--;
            }
        }
    }

    print_debug("Weekly trade and exit completed");
}

// Update targets adaptively and recompute posted prices
void Simulation::weekly_update_prices() {
    auto overhead_f = make_overhead(f1, slope);
    for (Shop& shop: shops) {
        if (shop.active) {
            shop.update_targets(alpha);
            shop.update_prices(C, overhead_f);
        }
    }
}

// Research process for a prospective owner r
// Should return an object ResearchResults, with targ0 and targ1 set as the local variables, and enter set as the return value of the function
ResearchResults Simulation::research(int idx) {
    Trader& trader = traders[idx];
    
    auto overhead_f = make_overhead(f1, slope);
    auto priceF = [this](double tr0, double tr1, double f_other) {
        return (tr1 - f_other - C > 0.0) ? ((tr1 - C - f_other) / tr0) : 0.0;
    };
    
    // Prepare targets and compute prices
    ResearchResults res{rng.uniform_int(xMax) + 1.0, rng.uniform_int(xMax) + 1.0, 0};
    double P0 = priceF(res.targ0, res.targ1, overhead_f(trader.get_demands()));
    double P1 = priceF(res.targ1, res.targ0, overhead_f(trader.get_supplies()));

    /*
        NOTE: Should be checking if comrade is same buyer or seller shop. 
        Since these are conditions for trading, simulation is considering
        non-matching roles as potential trade partners.
    */
    // Test with a comrade
    Trader& partner = traders[trader.any_comrade(produces)];
    double U = 0.0;
    double Ucomp = partner.utility(shops);
    // direct or indirect reachability check through fr's links
    if (trader.wants_to_trade_in(partner.get_demands())) {
        U = P0;
    } else {
        if (partner.get_buyer_shop() && trader.wants_to_trade_in(partner.get_buyer_shop()->get_good(partner.get_demands()))) {
            U = partner.get_buyer_shop()->get_price(partner.get_demands()) * P0;
        }
    }
    // Test with a soulmate
    if (U < Ucomp) {
        Trader& partner = traders[trader.soulmate(consumes)];
        // partner = traders[fr];
        Ucomp = partner.utility(shops);
        U = 0.0;
        if (trader.wants_to_trade_out(partner.get_supplies())) {
            U = P0;
        } else {
            if (partner.get_seller_shop() && trader.wants_to_trade_out(partner.get_seller_shop()->get_good(partner.get_supplies()))) {
                U = partner.get_seller_shop()->get_price_supply(partner.get_supplies()) * P0;
            }
        }
        if (U < Ucomp) {
            res.enter = 0;
            return res;
        }
    }

    {
        U = 0.0;
        Trader& partner = random_consumer(trader.get_supplies());
        Ucomp = partner.utility(shops);
        if (trader.wants_to_trade_in(partner.get_supplies())) {
            U = P1;
        } else {
            if (partner.get_seller_shop() && trader.wants_to_trade_in(partner.get_seller_shop()->get_good(partner.get_supplies()))) {
                U = partner.get_seller_shop()->get_price(trader.get_demands()) * P1;
            }
        }
        // Stranger who produces d[r]
        if (U < Ucomp) {
            Trader& partner = random_producer(trader.get_demands());
            Ucomp = partner.utility(shops);
            U = 0.0;
            if (trader.wants_to_trade_out(partner.get_demands())) {
                U = P1;
            } else {
                if (partner.get_buyer_shop() && trader.wants_to_trade_out(partner.get_buyer_shop()->get_good(partner.get_demands()))) {
                    U = P1 * partner.get_buyer_shop()->get_price(partner.get_demands());
                }
            }
            if (U < Ucomp) {
                res.enter = 0;
                return res;
            }
        }
    }
    res.enter = 1;
    return res;
}

Trader& Simulation::random_consumer(int good) {
    int k = 1 + rng.uniform_int(std::max(1, (int)consumes[good].size()) - 1);
    return traders[consumes[good][k - 1]];
}

Trader& Simulation::random_producer(int good) {
    int k = 1 + rng.uniform_int(std::max(1, (int)produces[good].size()) - 1);
    return traders[produces[good][k - 1]];
}

Shop& Simulation::random_shop() {
    return shops[rng.uniform_int(K) + 1];
}

void Simulation::addshop(const Trader* trader, Shop* shop, std::vector<int>& cand) {
    if (shop && shop->active && trader->is_compatible_with(shop)) {
        if (std::find(cand.begin(), cand.end(), shop->idx) == cand.end()) {
            // if (!(shop->idx != trader->get_buyer_idx() && shop->idx != trader->get_seller_idx())) {
            //     printf("hit!\n");
            // }
            cand.push_back(shop->idx);
        }
    }
}

/* NOTE: Line is generated only at the beginning of the run. Should it be randomized? */
void Simulation::lineup() {
    // random permutation of 1..m
    int line[m + 1];
    int start[m + 1];
    for (int i = 1; i <= m; ++i)
        start[i] = i;
    for (int j = 1; j <= m; ++j) {
        int k = rng.uniform_int(m - j + 1) + 1;
        line[j] = start[k];
        for (int i = k; i <= m - j; ++i) {
            start[i] = start[i + 1];
        }
    }
    trader_line.clear();
    for (int i = 0; i < m; ++i) {
        trader_line.push_back(&traders[line[i + 1]]);
    }
}

void Simulation::try_barter(Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    // iterate candidates from index 2 onward
    for (size_t idx = 2; idx < c.size(); ++idx) {
        Shop& shop = shops[c[idx]];
        if (trader.allows_barter_with(shop)) {
            double val = shop.get_price_supply(trader.get_supplies());
            if (val > std::max(eval.Ucomp, eval.Ubarter)) {
                eval.barter = &shop;
                eval.Ubarter = val;
                c.erase(c.begin() + idx);
                --idx;
            }
        }
    }
}

void Simulation::try_one(const Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    int s = trader.get_supplies();
    int d = trader.get_demands();
    // Track which side matches for current c[0] and c[1]
    for (size_t idx = 2; idx < c.size(); ++idx) {
        Shop& shop = shops[c[idx]];
        Shop& candidate_0 = shops[c[0]];
        Shop& candidate_1 = shops[c[1]];
        // improve outlet (sell s)
        if (shop.provides(s)) {
            if ((shop.get_good(s) == candidate_1.get_good(d)) || (candidate_0.get_price_supply(s) == 0.0)) {
                if (candidate_0.get_price_supply(s) < shop.get_price_supply(s)) {
                    if (shop.get_good(s) == candidate_1.get_good(d)) {
                        eval.Ucomp = shop.get_price_supply(s) * candidate_1.get_price_supply(s);
                    }
                    c[0] = c[idx];
                    eval.candidate_0 = &shop;
                    c.erase(c.begin() + idx);
                    --idx;
                }
            }
        }
        else {
            // improve source (buy d)
            if (shop.provides(d)) {
                if ((shop.get_good(d) == candidate_0.get_good(s)) || (candidate_1.get_price(d) == 0.0)) {
                    if (candidate_1.get_price(d) < shop.get_price(d)) {
                        if (shop.get_good(d) == candidate_0.get_good(s)) {
                            eval.Ucomp = (candidate_0.get_price_supply(s) * shop.get_price(d));
                        }
                        c[1] = c[idx];
                        eval.candidate_1 = &shop;
                        c.erase(c.begin() + idx);
                        --idx;
                    }
                }
            }
        }
    }
}

void Simulation::try_two(const Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    for (size_t ia = 2; ia < c.size(); ++ia) {
        int a = c[ia];
        if (shops[a].provides(trader.get_supplies())) {
            for (size_t ib = 2; ib < c.size(); ++ib) {
                if (ia == ib) continue;
                int b = c[ib];
                if (shops[b].provides(trader.get_demands())) {
                    // common intermediary condition
                    if (shops[a].get_good(trader.get_supplies()) == shops[b].get_good(trader.get_demands())) {
                        double val = shops[a].get_price_supply(trader.get_supplies()) * shops[b].get_price(trader.get_demands());
                        if (eval.Ucomp < val) {
                            eval.Ucomp = val;
                            eval.candidate_0 = &shops[a];
                            eval.candidate_1 = &shops[b];
                            c[0] = a;
                            c[1] = b;
                        }
                    }
                }
            }
        }
    }
}

int Simulation::calc1() {
    part = 0.0;
    moneytraders = 0.0;
    std::fill(usingmoney.begin(), usingmoney.end(), 0.0);

    for (Trader& trader : traders) {
        if (trader.get_seller_idx() > 0) {
            int a = trader.get_seller_idx();
            int b = trader.get_buyer_idx();
            int ma = (shops[a].g[0] == trader.get_supplies());
            int mb = (b > 0 && shops[b].g[0] == trader.get_demands());

            part += (shops[a].g[ma] == trader.get_demands()) ||
                    (shops[a].g[ma] == shops[b].g[mb]);

            if (b > 0 && shops[a].g[ma] == shops[b].g[mb]) {
                moneytraders += 1.0;
                usingmoney[shops[b].g[mb]] += 1.0;
            }
        }
    }

    if (fulldev == 0) {
        if (part >= 0.99 * m) {
            if (devcount == 0) devyear = t / 50;
            devcount++;
        } else {
            devcount = 0;
        }
        if (devcount >= persist) fulldev = 1;
    }

    usingmax = 0.0;
    moneygood = 0;
    for (int i = 1; i <= n; ++i) {
        if (usingmoney[i] > usingmax) {
            usingmax = usingmoney[i];
            moneygood = i;
        }
    }

    if (usingmoney[moneygood] >= 0.99 * Fmon) {
        if (endcount == 0) {
            monyear = t / 50;
        }
        endcount++;
    } else {
        endcount = 0;
    }

    return (endcount >= persist);
}

void Simulation::calc2() {
    auto overhead_f = make_overhead(f1, slope);
    if (monetary == 0) monyear = -1;
    if (fulldev == 0) devyear = -1;

    // count non-money active shops
    BS = 0;
    for (Shop& shop : shops) {
        if (shop.active) {
            if (shop.g[0] != moneygood && shop.g[1] != moneygood) BS++;
        }
    }

    W = 1.0 - (overhead_f(moneygood) + C) / bsize;
    if (W > 0.0) {
        for (int i = 1; i <= n; ++i) {
            Pinv[i] = ((m / n) - (overhead_f(i) + C)) /
                      ((m / n) - (n - 2) * (overhead_f(moneygood) + C));
        }
    }
    SurpSME = m - n * f1 - (slope / 2.0) * n * (n - 1) - (n - 2) * overhead_f(moneygood);

    rmse();

    Csurp = 0.0;
    for (Trader& trader: traders) Csurp += trader.utility(shops);

    Psurp = 0.0;
    Nshop = 0;
    for (Shop& shop : shops) {
        if (shop.active) {
            for (int h = 0; h < 2; ++h) {
                Psurp += (shop.y[h] - overhead_f(shop.g[h]) - shop.P[1 - h] * shop.y[1 - h]);
            }
            int own = shop.owner;
            int qown = (own > 0) ? traders[own].q : 0;
            if ((own > 0) && (shop.y[qown] > 1.0 || shop.y[1 - qown] > 0.0)) Nshop++;
        }
    }
}

void Simulation::rmse() {
    if (W > 0.0) {
        for (int h = 0; h < 2; ++h) {
            std::fill(vol[h].begin(), vol[h].end(), 0.0);
            std::fill(avp[h].begin(), avp[h].end(), 0.0);
        }

        // aggregate by pairs with moneygood
        for (Shop& shop: shops) {
            if (shop.active) {
                if (shop.g[0] == moneygood || shop.g[1] == moneygood) {
                    int ma = (shop.g[1] == moneygood);
                    int i = shop.g[1 - ma];
                    vol[0][i] += shop.y[1 - ma];
                    if (vol[0][i] > 0.0) {
                        avp[0][i] += (shop.y[1 - ma] / vol[0][i]) * (shop.P[1 - ma] - avp[0][i]);
                    }
                    vol[1][i] += shop.y[ma];
                    if (vol[1][i] > 0.0) {
                        avp[1][i] += (shop.y[ma] / vol[1][i]) * (shop.P[ma] - avp[1][i]);
                    }
                }
            }
        }

        for (int i = 1; i <= n; ++i) {
            if (i == moneygood) continue;
            avp[0][i] = (W != 0.0) ? (avp[0][i] / W) : 0.0;
            avp[1][i] = (Pinv[i] != 0.0) ? (avp[1][i] / Pinv[i]) : 0.0;
        }

        for (int h = 0; h < 2; ++h) {
            double sse = 0.0;
            int count = 0;
            for (int i = 1; i <= n; ++i) {
                if (i == moneygood) continue;
                double e = 1.0 - avp[h][i];
                sse += e * e;
                ++count;
            }
            R[h] = (count > 0) ? std::sqrt(sse / count) : -1.0;
        }
    } else {
        R[0] = R[1] = -1.0;
    }
}

void Simulation::report(int tt) {
    if (prtoscr != 0) {
        if (tt == -1) std::printf("***");
        std::printf("%6.0f %6.0f ", part, moneytraders);
        for (int b = 1; b <= 5 && b <= n; ++b) {
            std::printf("%6.0f ", usingmoney[b]);
        }
        std::printf("%6d %4d\n", (tt == -1 ? t : tt) / 50, NS);
    }
}

void Simulation::print_debug(std::string title) const {
    if (DEBUG) {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << title << std::endl;
        print_shops();
        print_traders();
        std::cout << "--------------------------------------------------" << std::endl;
    }
}

void Simulation::print_traders() const {
    std::cout << "Traders:" << std::endl;
    for (size_t i = 1; i < traders.size(); ++i) {
        std::cout << traders[i].to_string() << std::endl;
    }
}

void Simulation::print_shops() const {
    std::cout << "Shops:" << std::endl;
    for (size_t i = 1; i < shops.size(); ++i) {
        std::cout << "Shop " << i << ": " << shops[i].to_string() << std::endl;
    }
}
