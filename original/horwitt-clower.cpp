// HorwittClower_OOP.cpp
// Object-oriented C++17 refactor of the Howitt–Clower simulation described in the attached source,
// preserving algorithmic behavior while encapsulating state and utilities into classes.
//
// Build: g++ -std=gnu++17 -O2 -pipe -o howclow_oop HorwittClower_OOP.cpp
// Run:   ./howclow_oop
//
// Based on: horwitt-clower.c
// --------------------------------------------------------------------

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <functional>
#include <assert.h>

#define DEBUG 0

static int last_seen = -1;

struct PCG32 {
    uint64_t state{0};
    uint64_t inc{0}; // must be odd

    void seed(uint64_t initstate, uint64_t initseq) {
        state = 0U;
        inc = (initseq << 1u) | 1u;
        next();
        state += initstate;
        next();
    }

    uint32_t next() {
        uint64_t oldstate = state;
        state = oldstate * 6364136223846793005ULL + (inc | 1ULL);
        uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        uint32_t rot = (uint32_t)(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
    }

    // Uniform integer in [0, n-1]
    int uniform_int(int n) {
        uint32_t res = next();
        if (DEBUG) {
            std::cout << "RND: " << res << std::endl;
        }
        return (int)(res % (uint32_t)n);
    }

    // Uniform double in (0,1]
    double uniform01_inclusive() {
        // Match original granularity of 0.001 if desired; here use full precision
        // Return (0,1] by excluding exact 0.0
        int x = uniform_int(1000);
        return double(1 + x)/1000.0;
    }
};
// RNG
static PCG32 rng;

struct ResearchResults {
    double targ0;
    double targ1;
    int enter;
};

class Shop {
public:
    int idx{0};
    int active{0};         // 1 if active
    int g[2]{0,0};         // traded goods, ordered so g[0]<=g[1]
    int owner{0};          // trader id that owns this shop
    double P[2]{0.0,0.0};  // posted buying prices per side
    double y[2]{0.0,0.0};  // realized incomes this week per side
    double tr[2]{0.0,0.0}; // income targets per side

    std::string to_string() const {
        return "Shop{active=" + std::to_string(active) + ", g=[" + std::to_string(g[0]) + "," + std::to_string(g[1]) +
               "], owner=" + std::to_string(owner) + ", P=[" + std::to_string(P[0]) + "," + std::to_string(P[1]) +
               "], y=[" + std::to_string(y[0]) + "," + std::to_string(y[1]) + "], tr=[" + std::to_string(tr[0]) + "," + std::to_string(tr[1]) + "]}";
    }

    bool provides(int good) const {
        return g[0] == good || g[1] == good;
    }

    void clear() {
        active = 0;
        g[0] = g[1] = 0;
        owner = 0;
        P[0] = P[1] = 0.0;
        y[0] = y[1] = 0.0;
        tr[0] = tr[1] = 0.0;
    }

    int get_good(int side, bool opposite=false) const {
        if (opposite)
            return g[g[0] != side];
        else
            return g[g[0] == side];
    }

    double get_income(int side, bool opposite=false) const {
        if (opposite)
            return y[g[0] != side];
        else
            return y[g[0] == side];
    }

    void add_income(int side, double val, bool opposite=false) {
        int idx;
        if (opposite)
            idx = g[0] != side;
        else
            idx = g[0] == side;
        y[idx] += val;
    }

    double get_price(int side, bool opposite=false) const {
        if (opposite)
            return P[g[0] != side];
        else
            return P[g[0] == side];
    }

    void update_prices(double C, std::function<double(int)> overhead_f) {
        auto priceF = [C](double tr0, double tr1, double f_other) {
            return (tr1 - f_other - C > 0.0) ? ((tr1 - C - f_other) / tr0) : 0.0;
        };
        P[0]  = priceF(tr[0], tr[1], overhead_f(g[1]));
        P[1]  = priceF(tr[1], tr[0], overhead_f(g[0]));
    }

    void update_targets(double targ0, double targ1) {
        tr[0] = targ0;
        tr[1] = targ1;
    }

    // overhead is a pointer to a  function that returns a double given an integer
    bool is_profitable(std::function<double(int)> overhead = nullptr) const {
        return (y[0] - P[1] * y[1] - (overhead ? overhead(g[0]) : 0.0)) > 0 && (y[1] - P[0] * y[0] - (overhead ? overhead(g[1]) : 0.0)) > 0;
    }

    void reset_weekly_incomes() {
        if (active) {
            y[0] = y[1] = 0.0;
        }
    }
};

class Trader {
public:
    int idx{0};             // PROVISIONAL: to make the code compatible while it's  refactored
    int supplied_good{0};              // produced good
    int demanded_good{0};              // desired good
    int q{0};              // 0 if s<=d else 1 (orders shop's good pair)
    int outlet_idx{0};
    int source_idx{0};
    int familyshop{0};     // owned shop index (0 if none)

    std::string to_string() const {
        return "Trader{s=" + std::to_string(supplied_good) + ", d=" + std::to_string(demanded_good) + ", q=" + std::to_string(q) +
               ", sh=[" + std::to_string(outlet_idx) + "," + std::to_string(source_idx) + "], familyshop=" + std::to_string(familyshop) + "}";
    }

    void sever_links(Shop& shop) {
        if (outlet_idx == shop.idx) outlet_idx = 0;
        if (source_idx == shop.idx) source_idx = 0;
    }

    /*
     NOTE: ChatGPT says:
     Checking fr >= r is wrong: it skips any trader with id >= r (not just self) and can bias selection.
     The current index math (k+1) can go out of range when k is the last index.
     Simpler, correct approach: pick a random index in the full vector; if it equals r and size>1, pick the next (or wrap) — this guarantees a different trader without bias and avoids bounds errors.
    */
    int comrade(const std::vector<std::vector<int>>& produces) const {
        // someone else producing s[r] (the same production good)
        int k = rng.uniform_int(std::max(0, std::max(1, (int)produces[supplied_good].size()) - 1));
        int fr = produces[supplied_good][k];
        if (fr >= idx) fr = produces[supplied_good][k+1]; // skip self
        return fr;
    }

    int soulmate(const std::vector<std::vector<int>>& consumes) const {
        // someone else consuming d[r] (the  same consumption good)
        int k = rng.uniform_int(std::max(0, std::max(1, (int)consumes[demanded_good].size()) - 1));
        int fr = consumes[demanded_good][k];
        if (fr >= idx) fr = consumes[demanded_good][k+1]; // skip self
        return fr;
    }

    double utility(std::vector<Shop>& shops) const {
        // attainable consumption for r via current links
        double X = 0.0;
        if (outlet_idx > 0) {
            Shop& sell_shop = shops[outlet_idx];
            if (sell_shop.get_good(supplied_good) == demanded_good) {
                X = sell_shop.get_price(supplied_good, true);
            } else {
                if (source_idx > 0) {
                    Shop& buy_shop = shops[source_idx];
                    if (sell_shop.get_good(supplied_good) == buy_shop.get_good(demanded_good)) {
                        X = sell_shop.get_price(supplied_good, true) * buy_shop.get_price(demanded_good);
                    }
                }
            }
        }
        return X;
    }
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

    // static constexpr int n      = 10;   // goods
    // static constexpr int bsize  = 24;   // each (i!=j) type count
    // static constexpr int K      = 200;  // shop locations

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
    Simulation() {
        traders.resize(m + 1); // 1-based
        shops.resize(K + 1);   // 1-based
        produces.assign(n + 1, {});
        consumes.assign(n + 1, {});
        line.resize(m + 1, 0);
        usingmoney.assign(n + 1, 0.0);
        Pinv.assign(n + 1, 0.0);
        for (int h = 0; h < 2; ++h) {
            vol[h].assign(n + 1, 0.0);
            avp[h].assign(n + 1, 0.0);
        }
        init_static();
    }

    void run_all() {
        time(&firstbegin);
        stream = std::fopen("evol.fil", "w+");
        if (stream) {
            std::fprintf(stream, " run slope dev part mon mtrd usmx Mgd Csur Psur Esur Nsh NS BS RMS_0 RMS_1 dyr myr\n");
            std::fclose(stream);
        }

        for (int _s = FirstSlope; _s <= LastSlope; _s += 2) {
            slope = _s;
            for (run = 1; run <= numruns; ++run) {
                init_run();
                clock_begin = std::clock();

                for (t = 1; t <= T; ++t) {
                    do {
                        weekly_entry();
                    } while (NS == 0);
                    weekly_matching();
                    weekly_trade_and_exit();
                    weekly_update_prices();

                    if (t % (PRINT_LOOP_N * RptPer) == 0) {
                        report(t);
                        monetary = calc1();
                        if (monetary == 1) {
                            break;
                        }
                    }
                }

                calc2();

                clock_finish = std::clock();
                std::printf("Run number %d. Time elapsed: %.2f seconds.\n",
                            run,
                            (clock_finish - clock_begin) / (double)CLOCKS_PER_SEC);
                std::printf("Slope equals %-.0f, xMax equals %d\n\n", slope, xMax);

                stream = std::fopen("evol.fil", "a");
                if (stream) {
                    std::fprintf(stream,
                        "%5d %5.0f %3d %4.0f %3d %4.0f %4.0f %3d %4.0f %6.0f %5.0f %3d %3d %3d %6.3f %6.3f %4d %4d\n",
                        run, slope, fulldev, part, monetary, moneytraders, usingmax,
                        moneygood, Csurp, Psurp, SurpSME, Nshop, NS, BS, R[0], R[1], devyear, monyear);
                    std::fclose(stream);
                }

                // blank line between runs
                stream = std::fopen("evol.fil", "a");
                if (stream) {
                    std::fprintf(stream, "\n");
                    std::fclose(stream);
                }
            }

            // blank line between slope sweeps
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

private:
    // Derived/utility
    static constexpr int PRINT_LOOP_N = 6;

    double slope{16.0};

    // State
    std::vector<Trader> traders;
    std::vector<Shop> shops;

    std::vector<std::vector<int>> produces; // by good i: trader ids producing i
    std::vector<std::vector<int>> consumes; // by good j: trader ids desiring j

    std::vector<int> line; // weekly lineup (permutation of traders)
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

    // temps for search/comparison in a week
    // kept as locals in methods where possible

    // time & IO
    std::time_t firstbegin{}, finish{};
    std::clock_t clock_begin{}, clock_finish{};
    std::FILE* stream{nullptr};

    // Helpers that mirror macro logic
    // double overhead_f(int i) const {
    //     return f1 + (i - 1) * slope;
    // }

    std::function<double(int)> make_overhead(double f1, double slope) {
        return [f1, slope](int i) -> double {
            return f1 + (i - 1) * slope;
        };
    }

    double priceF(double tr0, double tr1, double f_other) const {
        if (DEBUG) {
            printf("priceF %.2f %.2f %.2f\n", tr0, tr1, f_other);
            printf("results -> %.2f\n", (tr1 - f_other - C > 0.0) ? ((tr1 - C - f_other) / tr0) : 0.0);
        }
        return (tr1 - f_other - C > 0.0) ? ((tr1 - C - f_other) / tr0) : 0.0;
    }

    void init_static() {
        // Build all (i,j) types with multiplicity bsize, i != j
        // Traders are 1..m enumerating all pairs with bsize copies
        int r = 0;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                if (i == j) continue;
                for (int k = 1; k <= bsize; ++k) {
                    ++r;
                    traders[r].idx = r;  // PROVISIONAL: to make the code compatible while it's  refactored
                    traders[r].supplied_good = i;
                    traders[r].demanded_good = j;
                    traders[r].q = (traders[r].supplied_good > traders[r].demanded_good);
                    produces[i].push_back(r);
                    consumes[j].push_back(r);
                }
            }
        }
        // Theoretical max number of "money traders": bsize*(n-2)*(n-1)
        Fmon = bsize * (n - 2.0) * (n - 1.0);
    }

    void init_run() {
        std::printf("Number Using Using Using Using Using Using \n");
        std::printf("Active Money good1 good2 good3 good4 good5 Year NS \n");

        endcount = 0; fulldev = 0; devyear = -1; monyear = -1; devcount = 0;

        rng.seed(RANDSEED + run - 1, RANDSEED + run - 1);

        clock_begin = std::clock();
        t = 0;
        NS = 0;

        for (Trader& trader : traders) {
            trader.outlet_idx = 0;
            trader.source_idx = 0;
            trader.familyshop = 0;
        }
        int i = 0;
        for (Shop& shop : shops) {
            shop.clear();
            shop.idx = i++;
        }

        lineup();
    }

    void lineup() {
        // random permutation of 1..m
        std::vector<int> start(m + 1);
        for (int i = 1; i <= m; ++i)
            start[i] = i;
        for (int j = 1; j <= m; ++j) {
            int k = rng.uniform_int(m - j + 1) + 1;
            line[j] = start[k];
            for (int i = k; i <= m - j; ++i) {
                start[i] = start[i + 1];
            }
        }
    }

    // Entry process: potential entrepreneur tries to open a shop
    void weekly_entry() {
        auto overhead_f = make_overhead(f1, slope);
        int r = rng.uniform_int(m) + 1; // prospective owner
        if (NS < K && traders[r].familyshop == 0) {
            ResearchResults ok = research(r);
            if (ok.enter > 0) {
                for (Shop& shop : shops) {
                    if (!shop.active && (&shop != &shops.front())) { // skip index 0
                        NS++;
                        shop.active = 1;
                        shop.g[1 - traders[r].q] = traders[r].demanded_good;
                        shop.g[traders[r].q]     = traders[r].supplied_good;
                        // initialize prices - targets are taken from previous research
                        shop.update_targets(ok.targ0, ok.targ1);
                        shop.update_prices(C, overhead_f);
                        // owner links
                        traders[r].outlet_idx = shop.idx;
                        traders[r].source_idx = 0;
                        traders[r].familyshop = shop.idx;
                        shop.owner = r;
                        break;
                    }
                }
            }
        }
        print_debug("Weekly entry completed");
    }

    // Matching: agents sample a small set of shops and adopt best links
    void weekly_matching() {
        if (DEBUG) {
            std::cout << "line: ";
            for (int i = 0; i <m; i++)
                std::cout << line[i + 1] << " ";
            std::cout << std::endl;
        }
        for (int i = 1; i <= m; i++) {
            int r = line[i];
            Trader& trader = traders[r];
            double U = trader.utility(shops);
            double psearch = (U > 0.0 ? lambda : 1.0);
            // Skip condition: random or already owns a shop
            if (rng.uniform01_inclusive() < psearch && trader.familyshop == 0) {
                // candidate initialization with current links
                std::vector<int> cand;
                cand.reserve(8);
                cand.push_back(trader.outlet_idx); // c[0]
                cand.push_back(trader.source_idx); // c[1]

                // add friend outlets/sources and one random shop
                // Trader& comrade_ = traders[comrade(trader)];
                Trader& comrade_ = traders[trader.comrade(produces)];
                addshop(comrade_, shops[comrade_.outlet_idx], cand);

                Trader& soulmate_ = traders[trader.soulmate(consumes)];
                addshop(soulmate_, shops[soulmate_.source_idx], cand);

                addshop(trader, shops[rng.uniform_int(K) + 1], cand);

                if (cand.size() > 2) {
                    double Ucomp = U;
                    int bestbarter = 0;
                    double Ubarter = 0.0;

                    try_barter(trader, cand, bestbarter, Ubarter, Ucomp);
                    if (cand.size() > 2 && (shops[cand[0]].g[1 - trader.q] != trader.demanded_good || shops[cand[0]].P[trader.q] == 0.0)) {
                        try_one(r, cand, Ucomp);
                    }
                    if (cand.size() > 2) {
                        try_two(r, cand, Ucomp);
                    }

                    if (Ucomp < Ubarter && bestbarter > 0) {
                        trader.outlet_idx = bestbarter;
                        trader.source_idx = 0;
                    } else {
                        // adopt c[0], c[1] as improved chain if any
                        trader.outlet_idx = cand[0];
                        trader.source_idx = cand[1];
                    }
                }
            }
            report_trader(trader);
        }
    }

    void report_trader(Trader const& trader) {
        if (DEBUG) {
            std::cout << trader.to_string() << std::endl;
        }
    }

    // Trade/accounting and stochastic exit of unprofitable shops
    void weekly_trade_and_exit() {
        print_debug("Weekly trade begins");

        // reset shop weekly incomes
        for (Shop& shop : shops) {
            shop.reset_weekly_incomes();
        }
        // tally incomes from adopted relationships
        for (Trader& trader : traders) {
            int a = trader.outlet_idx;
            if (a > 0) {
                int b = trader.source_idx;
                if (shops[a].get_good(trader.supplied_good) == trader.demanded_good) {
                    // direct barter
                    shops[a].add_income(trader.supplied_good, 1.0, true);
                } else if (b > 0 && shops[a].get_good(trader.supplied_good) == shops[b].get_good(trader.demanded_good)) {
                    // indirect via common intermediary
                    shops[a].add_income(trader.supplied_good, 1.0, true);
                    shops[b].add_income(trader.demanded_good, shops[a].get_price(trader.supplied_good, true));
                }
            }
        }
        // exit if unprofitable
        for (Shop& shop: shops) {
            if (shop.active) {
                if ((!shop.is_profitable(make_overhead(f1, slope))) && rng.uniform01_inclusive() <= theta) {
                    if (shop.owner) {
                        traders[shop.owner].familyshop = 0;
                    }
                    shop.clear();
                    NS--;

                    // sever links
                    for (Trader& trader: traders)
                        trader.sever_links(shop);
                }
            }
        }

        print_debug("Weekly trade and exit completed");
    }

    // Update targets adaptively and recompute posted prices
    void weekly_update_prices() {
        auto overhead_f = make_overhead(f1, slope);
        if (DEBUG) {
            for (Shop& shop: shops) {
                if (shop.active) {
                    printf("Shop %d: tr0=%.2f tr1=%.2f\n", shop.idx, shop.tr[0], shop.tr[1]);
                }
            }
        }
        for (Shop& shop: shops) {
            if (shop.active) {
                for (int h = 0; h < 2; ++h) {
                    shop.tr[h] += alpha * (shop.y[h] - shop.tr[h]);
                }
                shop.P[0] = priceF(shop.tr[0], shop.tr[1], overhead_f(shop.g[1]));
                shop.P[1] = priceF(shop.tr[1], shop.tr[0], overhead_f(shop.g[0]));
            }
        }
    }

    // Research process for a prospective owner r
    // Should return an object ResearchResults, with targ0 and targ1 set as the local variables, and enter set as the return value of the function
    ResearchResults research(int idx) {
        ResearchResults res;
        Trader& trader = traders[idx];
        auto overhead_f = make_overhead(f1, slope);
        // trial targets
        res.targ0 = rng.uniform_int(xMax) + 1.0;
        res.targ1 = rng.uniform_int(xMax) + 1.0;

        double P0 = priceF(res.targ0, res.targ1, overhead_f(trader.demanded_good));
        double P1 = priceF(res.targ1, res.targ0, overhead_f(trader.supplied_good));

        // Test with a comrade
        Trader& partner = traders[trader.comrade(produces)];
        double Ucomp = partner.utility(shops);
        double U = 0.0;
        // direct or indirect reachability check through fr’s links
        if (partner.demanded_good == trader.demanded_good) U = P0;
        else {
            int sh1 = partner.source_idx;
            if (sh1 > 0) {
                int m1 = (shops[sh1].g[0] == partner.demanded_good);
                if (shops[sh1].g[m1] == trader.demanded_good) U = P0 * shops[sh1].P[m1];
            }
        }
        // Test with a soulmate
        if (U < Ucomp) {
            U = 0;
            Trader& partner = traders[trader.soulmate(consumes)];
            // partner = traders[fr];
            Ucomp = partner.utility(shops);
            U = 0.0;
            if (partner.supplied_good == trader.supplied_good) {
                U = P0;
            } else {
                if (partner.outlet_idx > 0) {
                    int m0 = (shops[partner.outlet_idx].g[0] == partner.supplied_good);
                    if (shops[partner.outlet_idx].g[m0] == trader.supplied_good) U = shops[partner.outlet_idx].P[1 - m0] * P0;
                }
            }
            if (U < Ucomp) {
                res.enter = 0;
                return res;
            }
        }

        /*
        NOTE: random choices here look funky. I tried to remove the +1 -1,
        but random number sequence changes for reasons I cannot explain.
        I assume this is wrong. It's kept for now for comparison with former
        results, but it should be removed */
        // Stranger who likes s[r]
        {
            U = 0.0;
            int k = 1 + rng.uniform_int(std::max(1, (int)consumes[trader.supplied_good].size()) - 1);
            Trader &partner = traders[consumes[trader.supplied_good][k - 1]];
            Ucomp = partner.utility(shops);
            U = 0.0;
            if (partner.supplied_good == trader.demanded_good) {
                U = P1;
            } else {
                if (partner.outlet_idx > 0) {
                    int m0 = (shops[partner.outlet_idx].g[0] == partner.supplied_good);
                    if (shops[partner.outlet_idx].g[m0] == trader.demanded_good) U = shops[partner.outlet_idx].P[1 - m0] * P1;
                }
            }
            // Stranger who produces d[r]
            if (U < Ucomp) {
                int k = 1 + rng.uniform_int(std::max(1, (int)produces[trader.demanded_good].size()) - 1);
                Trader& partner = traders[produces[trader.demanded_good][k - 1]];
                Ucomp = partner.utility(shops);
                U = 0.0;
                if (partner.demanded_good == trader.supplied_good) {
                    U = P1;
                } else {
                    if (partner.source_idx > 0) {
                        int m1 = (shops[partner.source_idx].g[0] == partner.demanded_good);
                        if (shops[partner.source_idx].g[m1] == trader.supplied_good) U = P1 * shops[partner.source_idx].P[m1];
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

    void addshop(const Trader& trader, Shop& shop, std::vector<int>& cand) {
        if (shop.active && (shop.provides(trader.supplied_good) || shop.provides(trader.demanded_good))) {
            if (std::find(cand.begin(), cand.end(), shop.idx) == cand.end()) {
                cand.push_back(shop.idx);
            }
        }
    }

    void try_barter(Trader& trader, std::vector<int>& c, int& bestbarter, double& Ubarter, double& Ucomp) {
        bestbarter = 0;
        Ubarter = 0.0;
        // iterate candidates from index 2 onward
        for (size_t idx = 2; idx < c.size(); ++idx) {
            Shop& shop = shops[c[idx]];
            if (shop.provides(trader.supplied_good) && shop.provides(trader.demanded_good)) {
                double val = shop.P[trader.q];
                if (val > std::max(Ucomp, Ubarter)) {
                    bestbarter = shop.idx;
                    Ubarter = val;
                    c.erase(c.begin() + idx);
                    --idx;
                }
            }
        }
    }

    void try_one(int r, std::vector<int>& c, double& Ucomp) {
        int s = traders[r].supplied_good;
        int d = traders[r].demanded_good;
        int m0 = (shops[traders[r].outlet_idx].g[0] == s);
        int m1 = (shops[traders[r].source_idx].g[0] == d);
        // Track which side matches for current c[0] and c[1]
        for (size_t idx = 2; idx < c.size(); ++idx) {
            Shop& shop = shops[c[idx]];
            Shop& candidate_0 = shops[c[0]];
            Shop& candidate_1 = shops[c[1]];
            // improve outlet (sell s)
            if (shop.provides(s)) {
                if ((shop.get_good(s) == candidate_1.g[m1]) || (candidate_0.P[1 - m0] == 0.0)) {
                    if (candidate_0.P[1 - m0] < shop.get_price(s, true)) {
                        double candidate = 0;
                        if (c[1] > 0) {
                            candidate = shop.get_price(s, true) * candidate_1.P[1 - m1];
                        }
                        c[0] = c[idx];
                        m0 = (shop.g[0] == s);
                        Ucomp = (shop.get_good(s) == candidate_1.g[m1]) ? candidate : 0.0;
                        c.erase(c.begin() + idx);
                        --idx;
                    }
                }
            }
            else {
                // improve source (buy d)
                if (shop.provides(d)) {
                    if ((shop.get_good(d) == candidate_0.g[m0]) || (candidate_1.P[m1] == 0.0)) {
                        if (candidate_1.P[m1] < shop.get_price(d)) {
                            double candidate = (candidate_0.P[1 - m0]) * shop.get_price(d);
                            c[1] = c[idx];
                            m1 = (shop.g[0] == d);
                            last_seen = r;
                            if (DEBUG) {
                                printf("[3] SET m1 to %d\n", m1);
                            }
                            Ucomp = (shop.get_good(d) == candidate_0.g[m0]) ? candidate : 0.0;
                            c.erase(c.begin() + idx);
                            --idx;
                        }
                    }
                }
            }
        }
    }

    void try_two(int r, std::vector<int>& c, double& Ucomp) {
        int s = traders[r].supplied_good, d = traders[r].demanded_good;

        for (size_t ia = 2; ia < c.size(); ++ia) {
            int a = c[ia];
            if (shops[a].provides(s)) {
                for (size_t ib = 2; ib < c.size(); ++ib) {
                    if (ia == ib) continue;
                    int b = c[ib];
                    if (shops[b].provides(d)) {
                        // common intermediary condition
                        if (shops[a].get_good(s) == shops[b].get_good(d)) {
                            double val = shops[a].get_price(s, true) * shops[b].get_price(d);
                            if (Ucomp < val) {
                                Ucomp = val;
                                c[0] = a;
                                c[1] = b;
                            }
                        }
                    }
                }
            }
        }
    }

    int calc1() {
        part = 0.0;
        moneytraders = 0.0;
        std::fill(usingmoney.begin(), usingmoney.end(), 0.0);

        for (Trader& trader : traders) {
            if (trader.outlet_idx > 0) {
                int a = trader.outlet_idx;
                int b = trader.source_idx;
                int ma = (shops[a].g[0] == trader.supplied_good);
                int mb = (b > 0 && shops[b].g[0] == trader.demanded_good);

                part += (shops[a].g[ma] == trader.demanded_good) ||
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

    void calc2() {
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

    void rmse() {
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

    void report(int tt) {
        if (prtoscr != 0) {
            if (tt == -1) std::printf("***");
            std::printf("%6.0f %6.0f ", part, moneytraders);
            for (int b = 1; b <= 5 && b <= n; ++b) {
                std::printf("%6.0f ", usingmoney[b]);
            }
            std::printf("%6d %4d\n", (tt == -1 ? t : tt) / 50, NS);
        }
    }

    void print_debug(std::string title) const {
        if (DEBUG) {
            std::cout << "--------------------------------------------------" << std::endl;
            std::cout << title << std::endl;
            print_shops();
            print_traders();
            std::cout << "--------------------------------------------------" << std::endl;
        }
    }

    void print_traders() const {
        std::cout << "Traders:" << std::endl;
        for (size_t i = 1; i < traders.size(); ++i) {
            std::cout << traders[i].to_string() << std::endl;
        }
    }

    void print_shops() const {
        std::cout << "Shops:" << std::endl;
        for (size_t i = 1; i < shops.size(); ++i) {
            std::cout << "Shop " << i << ": " << shops[i].to_string() << std::endl;
        }
    }

private:
    int run{1};
    int t{0};
    int Nshop{0};
};

// -------------------------
// Entry point
// -------------------------
int main() {
    Simulation sim;
    sim.run_all();
    return 0;
}
