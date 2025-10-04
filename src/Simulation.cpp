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
#define CONTAINS(_COLLECTION, _ITEM)  (std::find(_COLLECTION.begin(), _COLLECTION.end(), _ITEM) != _COLLECTION.end())
 // Global RNG instance
PCG32 rng;

struct MatchEvaluation {
    double Ucomp = 0.0;
    double Ubarter = 0.0;
    Shop* barter = nullptr;
    Shop* candidate_seller = nullptr;
    Shop* candidate_buyer = nullptr;
};

Simulation::Simulation(SimulationInfo info) : info(info), RptPer(1), prtoscr(1), PRINT_LOOP_N(6) {
    traders.resize(info.m + 1); // 1-based
    shops.resize(info.K + 1);   // 1-based
    produces.assign(info.n + 1, {});
    consumes.assign(info.n + 1, {});
    usingmoney.assign(info.n + 1, 0.0);
    // Pinv.assign(info.n + 1, 0.0);
    init_static();
}

void Simulation::run_all() {
    // Main simulation loop: runs all parameter sweeps and time steps.
    // Simulation rule: Orchestrates initialization, weekly activities, reporting, and statistics collection.

    // Time and IO variables
    std::time_t firstbegin{}, finish{};
    std::clock_t clock_begin{}, clock_finish{};
    std::FILE* stream{nullptr};

    // Record the start time of the entire simulation
    time(&firstbegin);

    // Initialize the output file "evol.fil" with headers for evolutionary data
    stream = std::fopen("evol.fil", "w+");
    if (stream) {
        std::fprintf(stream, " run slope dev part mon mtrd usmx Mgd Csur Psur Esur Nsh NS BS RMS_0 RMS_1 dyr myr\n");
        std::fclose(stream);
    }

    // Loop over different slope values for parameter sweeps
    for (int _s = info.FirstSlope; _s <= info.LastSlope; _s += 2) {
        slope = _s;

        // For each slope, run multiple simulation runs
        for (int run = 1; run <= info.numruns; ++run) {
            runInfo.run = run;
            runInfo.Slope = slope;
            // Initialize the run-specific state
            init_run();
            clock_begin = std::clock();

            // Main simulation loop: run for T time steps
            for (int t = 1; t <= info.T; ++t) {
                runInfo.t = t;
                // Ensure at least one shop exists before proceeding
                do {
                    weekly_entry();
                } while (runInfo.NumberOfShops == 0);

                // Perform weekly activities: matching, trading, exit, and price updates
                weekly_matching();
                weekly_trade_and_exit();
                weekly_update_prices();

                // Periodically report progress and check for monetary equilibrium
                if (runInfo.t % (PRINT_LOOP_N * RptPer) == 0) {
                    runInfo.monetary = calc1();
                    report(runInfo.t);
                    if (runInfo.monetary == 1) {
                        printf("Monetary equilibrium reached\n");
                        break; // Exit early if monetary equilibrium is reached
                    }
                }
            }

            // Calculate final statistics for the run
            calc2();

            // Log timing information
            clock_finish = std::clock();
            std::printf("Run number %d. Time elapsed: %.2f seconds.\n",
                runInfo.run,
                (clock_finish - clock_begin) / (double)CLOCKS_PER_SEC);
            std::printf("Slope equals %-.0f, xMax equals %d\n\n", slope, info.xMax);

            // Append run results to the output file
            stream = std::fopen("evol.fil", "a");
            if (stream) {
                std::fprintf(stream,
                    "%5d %5.0f %3d %4.0f %3d %4.0f %4.0f %3d %4.0f %6.0f %5.0f %3d %3d %3d %6.3f %6.3f %4d %4d\n",
                    runInfo.run, slope, runInfo.fulldev, runInfo.part, runInfo.monetary, runInfo.moneytraders, runInfo.usingmax,
                    runInfo.moneygood, runInfo.Csurp, runInfo.Psurp, runInfo.SurpSME, runInfo.Nshop, runInfo.NumberOfShops, runInfo.BS, runInfo.R[0], runInfo.R[1], runInfo.devyear, runInfo.monyear);
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
                info.n, info.bsize, info.K, info.xMax, info.lambda, info.alpha, info.theta, info.C, info.numruns, info.persist, info.T, info.RANDSEED);
            time(&finish);
            std::fprintf(stream, "\n Started : %s", std::ctime(&firstbegin));
            std::fprintf(stream, " Finished: %s\n", std::ctime(&finish));
            std::fclose(stream);
        }
    }
}

std::function<double(int)> Simulation::make_overhead(double f1, double slope) {
    // Computes the overhead cost for a shop based on good index and slope.
    // Simulation rule: Overhead increases with good index and slope, affecting shop profitability.

    return [f1, slope](int i) -> double {
        return f1 + (i - 1) * slope;
        };
}

void Simulation::init_static() {
    // Initializes the static structure of the economy: traders, shops, and goods assignments.
    // Simulation rule: Sets up agent populations and their supply/demand relationships.

    // Build all (i,j) types with multiplicity bsize, i != j
    // Traders are 1..m enumerating all pairs with bsize copies
    int r = 0;
    for (int i = 1; i <= info.n; ++i) {
        for (int j = 1; j <= info.n; ++j) {
            if (i == j)
                continue;
            for (int k = 1; k <= info.bsize; ++k) {
                ++r;
                traders[r].idx = r; // PROVISIONAL: to make the code compatible while it's  refactored
                traders[r].set_supplies(i);
                traders[r].set_demands(j);
                produces[i].push_back(r);
                consumes[j].push_back(r);
            }
        }
    }
    // Theoretical max number of "money traders": bsize*(n-2)*(n-1)
    runInfo.Fmon = info.bsize * (info.n - 2.0) * (info.n - 1.0);
}

void Simulation::init_run() {
    // Initializes the state for a new simulation run, resetting agents and shops.
    // Simulation rule: Prepares the market for a fresh round of trading and adaptation.

    std::printf("Number Using Using Using Using Using Using \n");
    std::printf("Active Money good1 good2 good3 good4 good5 Year NS \n");

    runInfo.endcount = 0;
    runInfo.fulldev = 0;
    runInfo.devyear = -1;
    runInfo.monyear = -1;
    runInfo.devcount = 0;

    rng.seed(info.RANDSEED + runInfo.run - 1, info.RANDSEED + runInfo.run - 1);

    runInfo.t = 0;
    runInfo.NumberOfShops = 0;

    // TODO: Traders can be static, since it is a full set, but Shops could be a dynamic set.
    for (Trader& trader : traders) {
        trader.sever_links();
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
    // Handles entry of new shops by prospective entrepreneurs.
    // Simulation rule: Traders may open shops if conditions are favorable.

    auto overhead_f = make_overhead(info.f1, slope);
    int r = rng.uniform_int(info.m) + 1; // prospective owner
    Trader& trader = traders[r];
    if (runInfo.NumberOfShops < info.K && trader.get_familyshop() == nullptr) {
        ResearchResults ok = research(trader);
        if (ok.enter > 0) {
            for (Shop& shop : shops) {
                if (!shop.active && (&shop != &shops.front())) { // skip index 0
                    runInfo.NumberOfShops++;
                    trader.open_shop(shop);
                    shop.set_targets(ok.targ0, ok.targ1);
                    shop.update_prices(info.C, overhead_f);
                    break;
                }
            }
        }
    }
    print_debug("Weekly entry completed");
}

// Matching: agents sample a small set of shops and adopt best links
std::vector<MatchEvaluation>* Simulation::weekly_matching() {
    // Matches traders to shops based on utility and compatibility, updating links.
    // Simulation rule: Agents seek optimal trading partners each week.

    std::vector<MatchEvaluation>* response = new std::vector<MatchEvaluation>();

    for (Trader* trader : trader_line) {
        double U = trader->utility();
        double psearch = (U > 0.0 ? info.lambda : 1.0);
        // Skip condition: random or already owns a shop
        if (rng.uniform01_inclusive() < psearch && trader->get_familyshop() == nullptr) {
            struct MatchEvaluation eval;
            eval.Ucomp = U;
            eval.candidate_seller = trader->get_outlet();
            eval.candidate_buyer = trader->get_source();
            // candidate initialization with current links

            std::vector<int> cand;
            Shop* comrade_shop = traders[trader->trade_comrade(produces)].get_outlet(); // same production
            Shop* soulmate_shop = traders[trader->soulmate(consumes)].get_source(); // same consumption
            Shop* random_shop_ptr = &random_shop(); // random trader
            std::vector<Shop*> v = { eval.candidate_seller, eval.candidate_buyer };
            for (Shop* shop : { comrade_shop, soulmate_shop, random_shop_ptr }) {
                if (shop && shop->active && trader->is_compatible_with(shop) && !CONTAINS(v, shop) && !CONTAINS(cand, shop->idx))
                    cand.push_back(shop->idx);
            }

            if (cand.size() > 0) {
                // NOTE: Hack to emulate the fact that the original code addresses the item at index zero
                Shop zero;
                Shop* temp;
                if (eval.candidate_seller == nullptr)
                    temp = &zero;
                else
                    temp = eval.candidate_seller;
                try_barter(trader, cand, eval);
                if ((temp->get_the_other_good(trader->get_supplied_good()) != trader->get_demand_good()) ||
                    (temp->get_price_supply(trader->get_supplied_good()) == 0.0)) {
                    try_one(trader, cand, eval);
                }
                try_two(trader, cand, eval);

                if (eval.Ucomp < eval.Ubarter && eval.barter != nullptr) { // is barter better than trade?
                    trader->set_outlet(eval.barter);
                    trader->set_source(nullptr);
                }
                else { // trade wins
                    trader->set_outlet(eval.candidate_seller);
                    trader->set_source(eval.candidate_buyer);
                }
            }

            struct MatchEvaluation* temp = new MatchEvaluation(eval);
            temp->barter = eval.barter;
            temp->candidate_seller = eval.candidate_seller;
            temp->candidate_buyer = eval.candidate_buyer;
            temp->Ubarter = eval.Ubarter;
            temp->Ucomp = eval.Ucomp;
            response->push_back(*temp);
        }
        report_trader(trader);
    }
    return response;
}

void Simulation::report_trader(Trader const* trader) {
    if (DEBUG) {
        std::cout << trader->to_string() << std::endl;
    }
}

// Trade/accounting and stochastic exit of unprofitable shops
void Simulation::weekly_trade_and_exit() {
    // Executes weekly trading, income accounting, and stochastic exit of unprofitable shops.
    // Simulation rule: Shops may exit if not profitable, severing links with traders.

    print_debug("Weekly trade begins");

    // reset shop weekly incomes
    for (Shop& shop : shops) {
        shop.reset_weekly_incomes();
    }
    // tally incomes from adopted relationships
    for (Trader& trader : traders) {
        Shop* shop_a = trader.get_outlet();
        if (shop_a) {
            if (shop_a->get_the_other_good(trader.get_supplied_good()) == trader.get_demand_good()) {
                // direct barter
                shop_a->add_income(trader.get_supplied_good(), 1.0, true);
            }
            else {
                Shop* shop_b = trader.get_source();
                if (shop_b && shop_a->get_the_other_good(trader.get_supplied_good()) == shop_b->get_the_other_good(trader.get_demand_good())) {
                    // indirect via common intermediary
                    shop_a->add_income(trader.get_supplied_good(), 1.0, true);
                    shop_b->add_income(trader.get_demand_good(), shop_a->get_price_supply(trader.get_supplied_good()));
                }
            }
        }
    }
    // exit if unprofitable
    for (Shop& shop : shops) {
        if (shop.active) {
            if ((!shop.is_profitable(make_overhead(info.f1, slope))) && rng.uniform01_inclusive() <= info.theta) {
                // sever links
                for (Trader& trader : traders) {
                    trader.sever_links(&shop);
                }
                shop.clear();
                runInfo.NumberOfShops--;
            }
        }
    }

    print_debug("Weekly trade and exit completed");
}

// Update targets adaptively and recompute posted prices
void Simulation::weekly_update_prices() {
    // Updates shop targets and posted prices adaptively based on recent performance.
    // Simulation rule: Shops adjust strategies to maximize profit and survive.

    auto overhead_f = make_overhead(info.f1, slope);
    for (Shop& shop : shops) {
        if (shop.active) {
            shop.update_targets(info.alpha);
            shop.update_prices(info.C, overhead_f);
        }
    }
}

// NOTE: why is it that traders researched are not the same as traders matched later? If it's a random/statistical thing, should it be a larger sample?
// NOTE: is this process supposed to be paired two by two?
// Research process for a prospective owner r
// Should return an object ResearchResults, with targ0 and targ1 set as the local variables, and enter set as the return value of the function
ResearchResults Simulation::research(Trader& trader) {
    // Simulates the research process for a prospective shop owner, determining entry viability.
    // Simulation rule: Entry depends on utility comparisons with existing trading options.

    auto overhead_f = make_overhead(info.f1, slope);
    auto priceF = [this](double tr0, double tr1, double f_other) {
        return (tr1 - f_other - info.C > 0.0) ? ((tr1 - info.C - f_other) / tr0) : 0.0;
        };

    // Prepare targets and compute prices
    ResearchResults res{ rng.uniform_int(info.xMax) + 1.0, rng.uniform_int(info.xMax) + 1.0, 0 };
    double P0 = priceF(res.targ0, res.targ1, overhead_f(trader.get_demand_good()));
    double P1 = priceF(res.targ1, res.targ0, overhead_f(trader.get_supplied_good()));

    bool enter = true;

    // Test with a comrade
    Trader& partner = traders[trader.any_comrade(produces)];
    double U = 0.0;
    // NOTE: why partner utility not own?
    double Ucomp = partner.utility();
    // direct or indirect reachability check through fr's links
    if (trader.wants_to_trade_in(partner.get_demand_good())) {
        U = P0;
    }
    else {
        if (partner.get_source() &&
            trader.wants_to_trade_in(partner.get_source()->get_the_other_good(partner.get_demand_good()))) {
            U = partner.get_source()->get_price_demand(partner.get_demand_good()) * P0;
        }
    }

    if (U < Ucomp) {
        // Test with a soulmate
        Trader& partner2 = traders[trader.soulmate(consumes)];
        Ucomp = partner2.utility();
        U = 0.0;
        if (trader.wants_to_trade_out(partner2.get_supplied_good())) {
            U = P0;
        }
        else {
            if (partner2.get_outlet() &&
                trader.wants_to_trade_out(partner2.get_outlet()->get_the_other_good(partner2.get_supplied_good()))) {
                U = partner2.get_outlet()->get_price_supply(partner2.get_supplied_good()) * P0;
            }
        }
        if (U < Ucomp) {
            enter = false;
        }
    }

    if (enter) {
        Trader& partner3 = random_consumer(trader.get_supplied_good());
        Ucomp = partner3.utility();
        U = 0.0;
        if (trader.wants_to_trade_in(partner3.get_supplied_good())) {
            U = P1;
        }
        else {
            if (partner3.get_outlet() &&
                trader.wants_to_trade_in(partner3.get_outlet()->get_the_other_good(partner3.get_supplied_good()))) {
                U = partner3.get_outlet()->get_price_demand(trader.get_demand_good()) * P1;
            }
        }
    }

    if (enter && (U < Ucomp)) {
        // Stranger who produces d[r]
        Trader& partner4 = random_producer(trader.get_demand_good());
        Ucomp = partner4.utility();
        U = 0.0;
        if (trader.wants_to_trade_out(partner4.get_demand_good())) {
            U = P1;
        }
        else {
            if (partner4.get_source() &&
                trader.wants_to_trade_out(partner4.get_source()->get_the_other_good(partner4.get_demand_good()))) {
                U = P1 * partner4.get_source()->get_price_demand(partner4.get_demand_good());
            }
        }
        if (U < Ucomp) {
            enter = false;
        }
    }

    res.enter = enter ? 1 : 0;
    return res;
}

Trader& Simulation::random_consumer(int good) {
    // Selects a random consumer of a given good from the trader population.

    int k = 1 + rng.uniform_int(std::max(1, (int)consumes[good].size()) - 1);
    return traders[consumes[good][k - 1]];
}

Trader& Simulation::random_producer(int good) {
    // Selects a random producer of a given good from the trader population.

    int k = 1 + rng.uniform_int(std::max(1, (int)produces[good].size()) - 1);
    return traders[produces[good][k - 1]];
}

Shop& Simulation::random_shop() {
    // Selects a random shop from the shop population.

    return shops[rng.uniform_int(info.K) + 1];
}

/* NOTE: Line is generated only at the beginning of the run. Should it be randomized? */
void Simulation::lineup() {
    // random permutation of 1..m
    int line[info.m + 1];
    int start[info.m + 1];
    for (int i = 1; i <= info.m; ++i)
        start[i] = i;
    for (int j = 1; j <= info.m; ++j) {
        int k = rng.uniform_int(info.m - j + 1) + 1;
        line[j] = start[k];
        for (int i = k; i <= info.m - j; ++i) {
            start[i] = start[i + 1];
        }
    }
    trader_line.clear();
    for (int i = 0; i < info.m; ++i) {
        trader_line.push_back(&traders[line[i + 1]]);
    }
}

void Simulation::try_barter(const Trader* trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    // Attempts to match a trader with shops for barter opportunities.
    // Simulation rule: Barter is considered if compatible shops are available.

    // iterate candidates from index 2 onward
    for (size_t idx = 0; idx < c.size(); ++idx) {
        Shop& shop = shops[c[idx]];
        if (trader->allows_barter_with(shop)) {
            double val = shop.get_price_supply(trader->get_supplied_good());
            if (val > std::max(eval.Ucomp, eval.Ubarter)) {
                eval.barter = &shop;
                eval.Ubarter = val;
                c.erase(c.begin() + idx);
                --idx;
            }
        }
    }
}

// NOTE: Should NULL store be evaluated? PS: c[.] may be zero in the loop
// TODO: replace array for a proper queue
void Simulation::try_one(const Trader* trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    // Attempts to improve a trader's outlet or source by matching with candidate shops.
    // Simulation rule: Agents seek better trading terms through shop selection.

    Shop zero;
    int s = trader->get_supplied_good();
    int d = trader->get_demand_good();
    for (size_t idx = 0; idx < c.size(); ++idx) {
        Shop& shop = shops[c[idx]];
        Shop* candidate_0 = eval.candidate_seller ? eval.candidate_seller : &zero;
        Shop* candidate_1 = eval.candidate_buyer ? eval.candidate_buyer : &zero;

        // improve outlet (sell s)
        if (shop.provides(s)) {
            if (candidate_0->get_price_supply(s) < shop.get_price_supply(s)) {
                if ((shop.get_the_other_good(s) == candidate_1->get_the_other_good(d))) {
                    eval.Ucomp = shop.get_price_supply(s) * candidate_1->get_price_supply(s);
                    eval.candidate_seller = c[idx] > 0 ? &shop : nullptr;
                    c.erase(c.begin() + idx--); // Remove candidate from queue
                }
                else if (candidate_0->get_price_supply(s) == 0.0) {
                    eval.candidate_seller = c[idx] > 0 ? &shop : nullptr;
                    c.erase(c.begin() + idx--); // Remove candidate from queue
                }
            }
        }
        else if (shop.provides(d)) { // improve source (buy d)
            if (candidate_1->get_price_demand(d) < shop.get_price_demand(d)) {
                if ((shop.get_the_other_good(d) == candidate_0->get_the_other_good(s))) {
                    eval.Ucomp = candidate_0->get_price_supply(s) * shop.get_price_demand(d);
                    eval.candidate_buyer = c[idx] > 0 ? &shop : nullptr;
                    c.erase(c.begin() + idx--);
                }
                else if (candidate_1->get_price_demand(d) == 0.0) {
                    eval.candidate_buyer = c[idx] > 0 ? &shop : nullptr;
                    c.erase(c.begin() + idx--); // Remove candidate from queue
                }
            }
        }
    }
}

void Simulation::try_two(const Trader* trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    // Attempts to match a trader with two shops for indirect trade via a common intermediary.
    // Simulation rule: Indirect trade is considered for maximizing utility.

    for (size_t ia = 0; ia < c.size(); ++ia) {
        int a = c[ia];
        if (shops[a].provides(trader->get_supplied_good())) {
            for (size_t ib = 0; ib < c.size(); ++ib) {
                if (ia == ib)
                    continue;
                int b = c[ib];
                if (shops[b].provides(trader->get_demand_good())) {
                    // common intermediary condition
                    if (shops[a].get_the_other_good(trader->get_supplied_good()) == shops[b].get_the_other_good(trader->get_demand_good())) {
                        double val = shops[a].get_price_supply(trader->get_supplied_good()) * shops[b].get_price_demand(trader->get_demand_good());
                        if (eval.Ucomp < val) {
                            eval.Ucomp = val;
                            eval.candidate_seller = &shops[a];
                            eval.candidate_buyer = &shops[b];
                        }
                    }
                }
            }
        }
    }
}

int Simulation::calc1() {
    // Calculates monetary equilibrium and tracks the emergence of money in the simulation.
    // Simulation rule: Aggregates statistics to detect monetary phases and equilibrium.

    runInfo.part = 0.0;
    runInfo.moneytraders = 0.0;
    std::fill(usingmoney.begin(), usingmoney.end(), 0.0);

    for (Trader& trader : traders) {
        if (trader.get_outlet() != nullptr) {
            Shop zero;

            Shop* shop_a = trader.get_outlet();
            Shop* shop_b = trader.get_source();
            int ma = (shop_a->g[0] == trader.get_supplied_good());
            int mb = (shop_b != nullptr && shop_b->g[0] == trader.get_demand_good());

            runInfo.part += (shop_a->g[ma] == trader.get_demand_good()) ||
                    (shop_b && (shop_a->g[ma] == shop_b->g[mb])) ||
                    (shop_b == nullptr && shop_a->g[ma] == 0);

            if (shop_b != nullptr && shop_a->g[ma] == shop_b->g[mb]) {
                runInfo.moneytraders += 1.0;
                usingmoney[shop_b->g[mb]] += 1.0;
            }
        }
    }

    if (runInfo.fulldev == 0) {
        if (runInfo.part >= 0.99 * info.m) {
            if (runInfo.devcount == 0)
                runInfo.devyear = runInfo.t / 50;
            runInfo.devcount++;
        }
        else {
            runInfo.devcount = 0;
        }
        if (runInfo.devcount >= info.persist)
            runInfo.fulldev = 1;
    }

    runInfo.usingmax = 0.0;
    runInfo.moneygood = 0;
    for (int i = 1; i <= info.n; ++i) {
        if (usingmoney[i] > runInfo.usingmax) {
            runInfo.usingmax = usingmoney[i];
            runInfo.moneygood = i;
        }
    }

    if (usingmoney[runInfo.moneygood] >= 0.99 * runInfo.Fmon) {
        if (runInfo.endcount == 0) {
            runInfo.monyear = runInfo.t / 50;
        }
        runInfo.endcount++;
    }
    else {
        runInfo.endcount = 0;
    }

    return (runInfo.endcount >= info.persist);
}

void Simulation::calc2() {
    // Calculates final statistics for the simulation run, including surplus and shop counts.
    // Simulation rule: Summarizes market outcomes for analysis.
    std::vector<double> Pinv(info.n + 1, 0.0);

    auto overhead_f = make_overhead(info.f1, slope);
    if (runInfo.monetary == 0)
        runInfo.monyear = -1;
    if (runInfo.fulldev == 0)
        runInfo.devyear = -1;

    // count non-money active shops
    runInfo.BS = 0;
    for (Shop& shop : shops) {
        if (shop.active) {
            if (shop.g[0] != runInfo.moneygood && shop.g[1] != runInfo.moneygood)
                runInfo.BS++;
        }
    }

    runInfo.W = 1.0 - (overhead_f(runInfo.moneygood) + info.C) / info.bsize;
    if (runInfo.W > 0.0) {
        for (int i = 1; i <= info.n; ++i) {
            Pinv[i] = ((info.m / info.n) - (overhead_f(i) + info.C)) /
                ((info.m / info.n) - (info.n - 2) * (overhead_f(runInfo.moneygood) + info.C));
        }
    }
    runInfo.SurpSME = info.m - info.n * info.f1 - (slope / 2.0) * info.n * (info.n - 1) - (info.n - 2) * overhead_f(runInfo.moneygood);

    rmse(Pinv);

    runInfo.Csurp = 0.0;
    for (Trader& trader : traders)
        runInfo.Csurp += trader.utility();

    runInfo.Psurp = 0.0;
    runInfo.Nshop = 0;
    for (Shop& shop : shops) {
        if (shop.active) {
            for (int h = 0; h < 2; ++h) {
                runInfo.Psurp += (shop.y[h] - overhead_f(shop.g[h]) - shop.P[1 - h] * shop.y[1 - h]);
            }
            int own = shop.owner;
            int qown = (own > 0) ? traders[own].q : 0;
            if ((own > 0) && (shop.y[qown] > 1.0 || shop.y[1 - qown] > 0.0))
                runInfo.Nshop++;
        }
    }
}

void Simulation::rmse(std::vector<double> &Pinv) {
    std::vector<double> vol[2], avp[2];

    // Computes root mean square error for price and volume statistics.
    // Simulation rule: Measures market efficiency and price dispersion.
    for (int h = 0; h < 2; ++h) {
        vol[h].assign(info.n + 1, 0.0);
        avp[h].assign(info.n + 1, 0.0);
    }

    if (runInfo.W > 0.0) {
        for (int h = 0; h < 2; ++h) {
            std::fill(vol[h].begin(), vol[h].end(), 0.0);
            std::fill(avp[h].begin(), avp[h].end(), 0.0);
        }

        // aggregate by pairs with moneygood
        for (Shop& shop : shops) {
            if (shop.active) {
                if (shop.g[0] == runInfo.moneygood || shop.g[1] == runInfo.moneygood) {
                    int ma = (shop.g[1] == runInfo.moneygood);
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

        for (int i = 1; i <= info.n; ++i) {
            if (i == runInfo.moneygood)
                continue;
            avp[0][i] = (runInfo.W != 0.0) ? (avp[0][i] / runInfo.W) : 0.0;
            avp[1][i] = (Pinv[i] != 0.0) ? (avp[1][i] / Pinv[i]) : 0.0;
        }

        for (int h = 0; h < 2; ++h) {
            double sse = 0.0;
            int count = 0;
            for (int i = 1; i <= info.n; ++i) {
                if (i == runInfo.moneygood)
                    continue;
                double e = 1.0 - avp[h][i];
                sse += e * e;
                ++count;
            }
            runInfo.R[h] = (count > 0) ? std::sqrt(sse / count) : -1.0;
        }
    }
    else {
        runInfo.R[0] = runInfo.R[1] = -1.0;
    }
}

void Simulation::report(int tt) {
    // Reports simulation progress and statistics at specified intervals.

    if (prtoscr != 0) {
        if (tt == -1)
            std::printf("***");
        std::printf("%6.0f %6.0f ", runInfo.part, runInfo.moneytraders);
        for (int b = 1; b <= 5 && b <= info.n; ++b) {
            std::printf("%6.0f ", usingmoney[b]);
        }
        std::printf("%6d %4d\n", (tt == -1 ? runInfo.t : tt) / 50, runInfo.NumberOfShops);
    }
}

void Simulation::print_debug(std::string title) const {
    // Prints debug information about the simulation state.

    if (DEBUG) {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << title << std::endl;
        print_shops();
        print_traders();
        std::cout << "--------------------------------------------------" << std::endl;
    }
}

void Simulation::print_traders() const {
    // Prints the state of all traders for debugging.

    std::cout << "Traders:" << std::endl;
    for (size_t i = 1; i < traders.size(); ++i) {
        std::cout << traders[i].to_string() << std::endl;
    }
}

void Simulation::print_shops() const {
    // Prints the state of all shops for debugging.

    std::cout << "Shops:" << std::endl;
    for (size_t i = 1; i < shops.size(); ++i) {
        std::cout << "Shop " << i << ": " << shops[i].to_string() << std::endl;
    }
}
