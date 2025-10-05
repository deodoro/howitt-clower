/*
 * Simulation.cpp
 *
 * Author: Jose Deodoro <deodoro.filho@gmail.com> <jdeoliv@gmu.edu>
 *
 * This program is free software: you can    // NOTE[AUTOMATED]: Paper's overhead function is f(i) = s*(i-1) with no        // NOTE[AUTOMATED]: Paper Sectio            // NOTE[AUTOMATED]: Paper Section 7.2 sa    // NOTE[AUTOMATED]: Paper Section 7.2 i    // NOTE[AUTOMAT    // NOTE[AUTOMATED]: Paper Section 5.3 exit     // NOTE[AUTOMATED]: Paper Section 7.1 market research samples 4 specific t    // NOTE[AUTOMATED]: Paper Section 7.1 entry criterion: open shop if ≥1 prospective customer on EACH side
    // would choose this shop using Section 5 utility maximization rule.
    // CRITICAL: Entry requires customers on BOTH sides, not just one side.
    // Current logic tracks 'enter' flag - verify it implements this bilateral requirement correctly.ansactors for shop (i,j):
    // Prospective customers side 1: one with production=i, one with consumption=j
    // Prospective customers side 2: one with production=j, one with consumption=i
    // CURRENT DEVIATION: This code tests (comrade, soulmate, random consumer of i, random producer of j)
    // IMPACT: May affect entry decisions and thus emergence patterns. Consider refactoring for exact compliance.: shop exits with probability θ if operating surplus ≤ 0 in either commodity.
    // VERIFICATION: Current is_profitable() checks both π₀ₖ > 0 AND π₁ₖ > 0.
    // Paper suggests exit if either π₀ₖ ≤ 0 OR π₁ₖ ≤ 0. Check if logic matches intended behavior.D]: Verify income accounting matches paper's Section 5.3 operating surplus formulas:
    // π₀ₖ = y₀ₖ - p₁ₖ*y₁ₖ - f(g₀ₖ)  and  π₁ₖ = y₁ₖ - p₀ₖ*y₀ₖ - f(g₁ₖ)
    // Critical for correct exit decisions. Shop.is_profitable() should implement these exactly.cludes fallback for "widowed" transactors:
    // If no profitable relationship found, switch to any outlet offering positive price for production commodity.
    // VERIFICATION NEEDED: Current code has some fallback logic but requires review against paper specification.
    // This is critical for handling shop exits that leave traders stranded.pling order is:
            // 1. Random shop location k∈{1,...,K} (if occupied)
            // 2. Comrade's outlet (transactor with same production commodity)
            // 3. Soulmate's source (transactor with same consumption commodity)
            // CURRENT: Code includes all three but order in vector may not match paper sequence.
            // IMPACT: Minor - affects which shop is tested first, but all valid shops are considered. 7.1 - market research tests 4 specific transactor types:
        // Side 1: one with production=i, one with consumption=j
        // Side 2: one with production=j, one with consumption=i
        // CURRENT ISSUE: research() tests different types (comrade, soulmate, random consumer/producer).
        // RECOMMENDATION: Verify research() logic matches paper's exact 4-type specification.intercept term.
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

 // Loop interval between state dumps
#define PRINT_LOOP_N 6
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

Simulation::Simulation(SimulationInfo info) : info(info) {
    traders.resize(info.m + 1); // 1-based
    shops.resize(info.K + 1);   // 1-based
    produces.assign(info.n + 1, {});
    consumes.assign(info.n + 1, {});
    init_static();
}

void Simulation::run_all() {
    // Main simulation loop: runs all parameter sweeps and time steps.
    // Simulation rule: Orchestrates initialization, weekly activities, reporting, and statistics collection.

    // Time variables
    std::clock_t clock_begin{}, clock_finish{};

    runs_per_slope.clear(); // clear any previous data

    // Loop over different slope values for parameter sweeps
    for (int _s = info.FirstSlope; _s <= info.LastSlope; _s += 2) {
        slope = _s;

        std::vector<RunInfo*> runs_for_slope;

        // For each slope, run multiple simulation runs
        for (int run = 1; run <= info.numruns; ++run) {
            RunInfo runInfo(info.n);
            runInfo.run = run;
            runInfo.Slope = slope;
            // Initialize the run-specific state
            init_run(runInfo);
            clock_begin = std::clock();

            // Main simulation loop: run for T time steps
            for (int t = 1; t <= info.T; ++t) {
                runInfo.t = t;
                // Ensure at least one shop exists before proceeding
                do {
                    weekly_entry(runInfo);
                } while (runInfo.NumberOfShops == 0);

                // Perform weekly activities: matching, trading, exit, and price updates
                  // NOTE[AUTOMATED] Paper’s Stage 7.1 "Entry" is once per week; Stage 7.2 "Shopping" processes all traders in a fixed weekly order. Keep this order.
                weekly_matching();
                weekly_trade_and_exit(runInfo);
                weekly_update_prices();

                // Periodically report progress and check for monetary equilibrium
                if (runInfo.t % PRINT_LOOP_N == 0) {
                    runInfo.monetary = calc1(runInfo);
                    runs_for_slope.push_back(runInfo.clone());
                    runInfo.report();
                    if (runInfo.monetary) {
                        printf("Monetary equilibrium reached\n");
                        break; // Exit early if monetary equilibrium is reached
                    }
                }
            }

            // Calculate final statistics for the run
            calc2(runInfo);

            // Log timing information
            clock_finish = std::clock();
            runInfo.time_spent = (clock_finish - clock_begin) / (double)CLOCKS_PER_SEC;
            std::printf("Run number %d. Time elapsed: %.2f seconds.\n",
                runInfo.run,
                (clock_finish - clock_begin) / (double)CLOCKS_PER_SEC);
            std::printf("Slope equals %-.0f, xMax equals %d\n\n", slope, info.xMax);

            // Store the runInfo for this run
            runs_for_slope.push_back(runInfo.clone());
        }

        // Store all runs for this slope
        runs_per_slope.push_back(runs_for_slope);
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
                traders[r].idx = r; // Trader's unique key
                traders[r].set_supplies(i);
                traders[r].set_demands(j);
                produces[i].push_back(r);
                consumes[j].push_back(r);
            }
        }
    }
    int u =0;
    for (Shop& shop : shops) {
        shop.idx = u++; // Shop's unique key
    }
}

void Simulation::init_run(RunInfo& runInfo) {
    // Initializes the state for a new simulation run, resetting agents and shops.
    // Simulation rule: Prepares the market for a fresh round of trading and adaptation.

    rng.seed(info.RANDSEED + runInfo.run - 1, info.RANDSEED + runInfo.run - 1);

    std::printf("Number  Using  Using  Using  Using  Using  Using \n");
    std::printf("Active  Money  good1  good2  good3  good4  good5   Year   NS \n");

    // Theoretical max number of "money traders": bsize*(n-2)*(n-1)
    runInfo.Fmon = info.bsize * (info.n - 2.0) * (info.n - 1.0); // move to Info
    runInfo.endcount = 0; // Move to info
    runInfo.fulldev = 0;
    runInfo.devyear = -1; // Move to info
    runInfo.monyear = -1; // Move to info
    runInfo.devcount = 0;  // Move to info
    runInfo.t = 0;
    runInfo.NumberOfShops = 0;

    // TODO: Traders can be static, since it is a full set, but Shops could be a dynamic set.
    for (Trader& trader : traders) { trader.sever_links(); }
    for (Shop& shop : shops) { shop.clear(); }

    // Trader order is shuffled every run
    lineup();
}

// Entry process: potential entrepreneur tries to open a shop
void Simulation::weekly_entry(RunInfo& runInfo) {
    // Handles entry of new shops by prospective entrepreneurs.
    // Simulation rule: Traders may open shops if conditions are favorable.

    auto overhead_f = make_overhead(info.f1, slope);
    int r = rng.uniform_int(info.m) + 1; // prospective owner
    Trader& trader = traders[r];
    if (runInfo.NumberOfShops < info.K && trader.get_familyshop() == nullptr) {
        ResearchResults ok = research(trader);
        // NOTE[AUTOMATED] Paper’s 7.1 requires success with “at least one prospective customer on each side” using the 4 fixed types (prod-i, cons-j, prod-j, cons-i). Confirm research() replicates this exact rule before opening.
        if (ok.enter > 0) {
            for (Shop& shop : shops) {
                if (!shop.active && (&shop != &shops.front())) { // skip index 0
                    runInfo.NumberOfShops++;
                    trader.open_shop(shop);
                    shop.set_targets(ok.targ0, ok.targ1);
                    shop.update_prices(info.C, overhead_f);
                    // NOTE[AUTOMATED] Verify posted prices use p0 = ((tr1 - f(g1) - C)/tr0)^+ and p1 analogously (positive-part). Also check that implied p0*p1 < 1 where relevant.
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

            // NOTE: order of test does not match paper (random should be first)
            // NOTE[AUTOMATED] Paper’s 7.2 sampling order: (1) a random location k (if occupied), (2) comrade’s outlet (same production), (3) soulmate’s source (same consumption). Reorder push to cand accordingly for fidelity.
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

    // NOTE[AUTOMATED] Paper also forces a switch to any outlet with positive price if current relationships are unprofitable (“widowed” cases). Ensure your fallback covers this.
    return response;
}

void Simulation::report_trader(Trader const* trader) {
    if (DEBUG) {
        std::cout << trader->to_string() << std::endl;
    }
}

// Trade/accounting and stochastic exit of unprofitable shops
void Simulation::weekly_trade_and_exit(RunInfo& runInfo) {
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

    // NOTE[AUTOMATED] Confirm add_income and internal accounting ensure the shop’s two side-flows imply
    // π0k = y0k - p1k*y1k - f(g0k) and π1k = y1k - p0k*y0k - f(g1k) as in §5.3. Otherwise exit decisions will differ.

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

    // NOTE[AUTOMATED] Paper’s exit rule: if either operating surplus ≤ 0, exit with probability θ that week. Ensure is_profitable() implements “positive in both commodities” exactly.

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

// NOTE: ASK PROFESSOR AXTELL: Why is it that traders evaluate the market positively given their utility compared to a partner? I understand
// it is the trade in essence, but why does it require the utility of trader be greater than the partner's?

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


    // NOTE: Double check the evaluation order. In the paper, it says:
    // (prod-i), (cons-j) on one side; (prod-j), (cons-i) on the other.
    // Code uses (comrade), (soulmate), (random consumer of i), (random producer of j).
    // Currently assuming they are equivalent because the code requires AND condition using both sides


    /* Testing if anyone would be willing to trade trader's demanded good*/
    // Test with a comrade
    Trader& comrade = traders[trader.any_comrade(produces)];
    double U = 0.0;
    // NOTE: why comrade utility not own?
    double Ucomp = comrade.utility();
    // direct or indirect reachability check through fr's links
    if (trader.wants_to_trade_in(comrade.get_demand_good())) {  // Comrade produces and consumes the same as  trader
        U = P0;
    }
    else { // Comrade produces the same, but consumes another. Can they trade indirectly?
        if (comrade.get_source() &&
            trader.wants_to_trade_in(comrade.get_source()->get_the_other_good(comrade.get_demand_good()))) {
            U = comrade.get_source()->get_price_demand(comrade.get_demand_good()) * P0;
        }
    }

    if (U < Ucomp) { // Comrades  don't buy the demand, would soulmates supply something at least (same consumption)?
        Trader& soulmate = traders[trader.soulmate(consumes)];
        Ucomp = soulmate.utility();
        U = 0.0;
        if (trader.wants_to_trade_out(soulmate.get_supplied_good())) { // Soulmate produces and consumes the same
            U = P0;
        }
        else { // Soulmate produces the same, but consumes another. Can they trade indirectly?
            if (soulmate.get_outlet() &&
                trader.wants_to_trade_out(soulmate.get_outlet()->get_the_other_good(soulmate.get_supplied_good()))) {
                U = soulmate.get_outlet()->get_price_supply(soulmate.get_supplied_good()) * P0;
            }
        }
        if (U < Ucomp) {
            enter = false;
        }
    }

    /* Testing if anyone would be willing to trade trader's supplied good*/
    if (enter) { // Utility does not increase with the traders, would it work with a random consumer?
        Trader& any_consumer = random_consumer(trader.get_supplied_good());
        Ucomp = any_consumer.utility();
        U = 0.0;
        if (trader.wants_to_trade_in(any_consumer.get_supplied_good())) {
            U = P1;
        }
        else {
            if (any_consumer.get_outlet() &&
                trader.wants_to_trade_in(any_consumer.get_outlet()->get_the_other_good(any_consumer.get_supplied_good()))) {
                U = any_consumer.get_outlet()->get_price_demand(trader.get_demand_good()) * P1;
            }
        }

        if (U < Ucomp) {
            // Stranger who produces d[r]
            Trader& any_producer = random_producer(trader.get_demand_good());
            Ucomp = any_producer.utility();
            U = 0.0;
            if (trader.wants_to_trade_out(any_producer.get_demand_good())) {
                U = P1;
            }
            else {
                if (any_producer.get_source() &&
                    trader.wants_to_trade_out(any_producer.get_source()->get_the_other_good(any_producer.get_demand_good()))) {
                    U = P1 * any_producer.get_source()->get_price_demand(any_producer.get_demand_good());
                }
            }
            if (U < Ucomp) {
                enter = false;
            }
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
    // NOTE[AUTOMATED] Paper picks a random location k∈{1,…,K} and includes it iff occupied. random_shop() should emulate that (possibly returning an inactive slot and being filtered out above).
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

    // NOTE[AUTOMATED] Paper does not *prioritize* direct barter; it simply considers all feasible sets and chooses the best via the Section 5 rule. Ensure this heuristic does not bias away from monetary patterns.
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

    // NOTE[AUTOMATED] Ensure “candidate improvement” respects the Section 5 choice rule (maximize weekly consumption), not a greedy local swap that could violate global best among the sample.
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

int Simulation::calc1(RunInfo& runInfo) {
    // Calculates monetary equilibrium and tracks the emergence of money in the simulation.
    // Simulation rule: Aggregates statistics to detect monetary phases and equilibrium.
    runInfo.part = 0.0;
    runInfo.moneytraders = 0.0;
    runInfo.usingmoney.assign(info.n + 1, 0.0);

    for (Trader& trader : traders) {
        if (trader.get_outlet() != nullptr) {
            Shop zero;

            Shop* shop_a = trader.get_outlet();
            Shop* shop_b = trader.get_source();
            int ma = (shop_a->g[0] == trader.get_supplied_good());
            int mb = (shop_b != nullptr && shop_b->g[0] == trader.get_demand_good());

            runInfo.part += (shop_a->g[ma] == trader.get_demand_good()) ||
                    (shop_b != nullptr && (shop_a->g[ma] == shop_b->g[mb])) ||
                    (shop_b == nullptr && (shop_a->g[ma] == 0));

            if (shop_b != nullptr && shop_a->g[ma] == shop_b->g[mb]) {
                runInfo.moneytraders += 1.0;
                runInfo.usingmoney[shop_b->g[mb]] += 1.0;
            }
        }
    }

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

    runInfo.usingmax = 0.0;
    runInfo.moneygood = 0;
    for (int i = 1; i <= info.n; ++i) {
        if (runInfo.usingmoney[i] > runInfo.usingmax) {
            runInfo.usingmax = runInfo.usingmoney[i];
            runInfo.moneygood = i;
        }
    }

    if (runInfo.usingmoney[runInfo.moneygood] >= 0.99 * runInfo.Fmon) {
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

void Simulation::calc2(RunInfo& runInfo) {
    // Calculates final statistics for the simulation run, including surplus and shop counts.
    // Simulation rule: Summarizes market outcomes for analysis.
    std::vector<double> Pinv(info.n + 1, 0.0);

    auto overhead_f = make_overhead(info.f1, slope);

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

    rmse(runInfo, Pinv);

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

void Simulation::rmse(RunInfo& runInfo, std::vector<double> &Pinv) {
    // Computes root mean square error for price and volume statistics.
    // Simulation rule: Measures market efficiency and price dispersion.
    std::vector<double> vol[2], avp[2];
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
