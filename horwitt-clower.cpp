// HorwittClower_OOP.cpp
// Object-oriented C++17 refactor of the Howitt–Clower simulation described in the attached source,
// preserving algorithmic behavior while encapsulating state and utilities into classes.
//
// Notes:
// - Indices remain 1-based where the original program relied on that convention.
// - Constants mirror the attached file; pricing F(tr0,tr1,f_other) and overhead f(i) follow the same rules.
// - RNG uses a PCG32 class mirroring the original functions.
// - "Shop" and "Trader" are explicit types; simulation state is isolated in Simulation.
// - File output to evol.fil and periodic reporting preserved (guarded by prtoscr).
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
#include <assert.h>

#define DEBUG 0
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

struct ResearchResults {
    double targ0;
    double targ1;
    int enter;
};

class Trader {
public:
    int idx{0};             // PROVISIONAL: to make the code compatible while it's  refactored
    int s{0};              // produced good
    int d{0};              // desired good
    int q{0};              // 0 if s<=d else 1 (orders shop's good pair)
    int sell{0};
    int buy{0};
    int familyshop{0};     // owned shop index (0 if none)

    std::string to_string() const {
        return "Trader{s=" + std::to_string(s) + ", d=" + std::to_string(d) + ", q=" + std::to_string(q) +
               ", sh=[" + std::to_string(sell) + "," + std::to_string(buy) + "], familyshop=" + std::to_string(familyshop) + "}";
    }
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
};

class Simulation {
public:
    // Configuration constants (mirror attached code)
    static constexpr int numruns   = 1;
    static constexpr int FirstSlope = 16;
    static constexpr int LastSlope  = 18;

    static constexpr int n      = 4;   // goods
    static constexpr int bsize  = 2;   // each (i!=j) type count
    static constexpr int K      = 4;  // shop locations

    // static constexpr int n      = 10;   // goods
    // static constexpr int bsize  = 24;   // each (i!=j) type count
    // static constexpr int K      = 200;  // shop locations

    static constexpr int xMax   = 200;
    static constexpr double lambda = 0.05;
    static constexpr double alpha  = 0.25;
    static constexpr double theta  = 0.01;
    static constexpr double C      = 5.0;
    static constexpr double f1     = 0.0;
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

    int global_m0; // HACK NECESSARY TO MAKE IT COMPATIBLE WITH FORMER VERSION
    int global_m1; // HACK NECESSARY TO MAKE IT COMPATIBLE WITH FORMER VERSION

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

    // RNG
    PCG32 rng;

    // Helpers that mirror macro logic
    double overhead_f(int i) const {
        return f1 + (i - 1) * slope;
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
                    traders[r].s = i;
                    traders[r].d = j;
                    traders[r].q = (traders[r].s > traders[r].d);
                    produces[i].push_back(r);
                    consumes[j].push_back(r);
                }
            }
        }
        // Theoretical max number of "money traders": bsize*(n-2)*(n-1)
        Fmon = bsize * (n - 2.0) * (n - 1.0);
        global_m0 = 0;
        global_m1 = 0;
    }

    void init_run() {
        std::printf("Number Using Using Using Using Using Using \n");
        std::printf("Active Money good1 good2 good3 good4 good5 Year NS \n");

        endcount = 0; fulldev = 0; devyear = -1; monyear = -1; devcount = 0;

        rng.seed(RANDSEED + run - 1, RANDSEED + run - 1);

        clock_begin = std::clock();
        t = 0;
        NS = 0;

        for (int r = 1; r <= m; ++r) {
            traders[r].sell = 0;
            traders[r].buy = 0;
            traders[r].familyshop = 0;
        }
        for (int k = 1; k <= K; ++k) {
            shops[k] = Shop{};
            shops[k].idx = k; // PROVISIONAL: to make the code compatible while it's  refactored
        }

        lineup();
    }

    void lineup() {
        // random permutation of 1..m
        std::vector<int> start(m + 1);
        for (int i = 1; i <= m; ++i) start[i] = i;
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
        int r = rng.uniform_int(m) + 1; // prospective owner
        if (NS < K && traders[r].familyshop == 0) {
            ResearchResults ok = research(r);
            if (ok.enter > 0) {
                // open new shop
                NS++;
                int k = 1;
                while (k <= K && shops[k].active == 1) k++;
                if (k > K) return; // no slot (shouldn't occur if NS<K)

                shops[k].active = 1;
                shops[k].g[1 - traders[r].q] = traders[r].d;
                shops[k].g[traders[r].q]     = traders[r].s;

                // initialize prices - targets are taken from previous research
                shops[k].tr[0] = ok.targ0;
                shops[k].tr[1] = ok.targ1;
                shops[k].P[0]  = priceF(shops[k].tr[0], shops[k].tr[1], overhead_f(shops[k].g[1]));
                shops[k].P[1]  = priceF(shops[k].tr[1], shops[k].tr[0], overhead_f(shops[k].g[0]));

                // owner links
                traders[r].sell = k;
                traders[r].buy = 0;
                traders[r].familyshop = k;
                shops[k].owner = r;
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
            double U = utility(traders[r]);
            double psearch = (U > 0.0 ? lambda : 1.0);
            // Skip condition: random or already owns a shop
            double rr = rng.uniform01_inclusive();
            if (rr < psearch && traders[r].familyshop == 0) {
                // candidate initialization with current links
                std::vector<int> cand;
                cand.reserve(8);
                cand.push_back(traders[r].sell); // c[0]
                cand.push_back(traders[r].buy); // c[1]

                // add friend outlets/sources and one random shop
                int fr1 = comrade(traders[r]);
                addshop(traders[fr1], shops[traders[fr1].sell], cand);

                int fr2 = soulmate(traders[r]);
                addshop(traders[fr2], shops[traders[fr2].buy], cand);

                addshop(traders[r], shops[rng.uniform_int(K) + 1], cand);

                if (cand.size() > 2) {
                    double Ucomp = U;
                    int bestbarter = 0;
                    double Ubarter = 0.0;

                    try_barter(r, cand, bestbarter, Ubarter, Ucomp);
                    if (cand.size() > 2 && (shops[cand[0]].g[1 - traders[r].q] != traders[r].d || shops[cand[0]].P[traders[r].q] == 0.0)) {
                        try_one(r, cand, Ucomp);
                    }
                    if (cand.size() > 2) {
                        try_two(r, cand, Ucomp);
                    }

                    if (Ucomp < Ubarter && bestbarter > 0) {
                        traders[r].sell = bestbarter;
                        traders[r].buy = 0;
                    } else {
                        // adopt c[0], c[1] as improved chain if any
                        traders[r].sell = cand[0];
                        traders[r].buy = cand[1];
                    }
                }
            }
            report_trader(traders[r]);
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
        for (int k = 1; k <= K; ++k) if (shops[k].active) {
            shops[k].y[0] = shops[k].y[1] = 0.0;
        }
        // tally incomes from adopted relationships
        for (int r = 1; r <= m; ++r) {
            int a = traders[r].sell;
            if (a > 0) {
                int b = traders[r].buy;
                int ma = (shops[a].g[0] == traders[r].s);
                int mb = (shops[b].g[0] == traders[r].d);

                if (shops[a].g[ma] == traders[r].d) {
                    // direct barter
                    shops[a].y[1 - ma] += 1.0;
                } else if (b > 0 && shops[a].g[ma] == shops[b].g[mb]) {
                    // indirect via common intermediary
                    shops[a].y[1 - ma] += 1.0;
                    shops[b].y[mb]     += shops[a].P[1 - ma];
                }
            }
        }
        // exit if unprofitable
        for (int k = 1; k <= K; ++k) if (shops[k].active) {
            bool unprof0 = (shops[k].y[0] - shops[k].P[1] * shops[k].y[1] - overhead_f(shops[k].g[0]) <= 0.0);
            bool unprof1 = (shops[k].y[1] - shops[k].P[0] * shops[k].y[0] - overhead_f(shops[k].g[1]) <= 0.0);
            if ((unprof0 || unprof1) && rng.uniform01_inclusive() <= theta) {
                int own = shops[k].owner;
                shops[k] = Shop{};
                shops[k].idx = k;
                NS--;
                if (own > 0) traders[own].familyshop = 0;

                // sever links
                for (int r = 1; r <= m; ++r) {
                    if (traders[r].sell == k) traders[r].sell = 0;
                    if (traders[r].buy == k) traders[r].buy = 0;
                }
            }
        }


        print_debug("Weekly trade and exit completed");
    }

    // Update targets adaptively and recompute posted prices
    void weekly_update_prices() {
        if (DEBUG) {
            for (int k = 1; k < K; ++k) if (shops[k].active) {
                printf("Shop %d: tr0=%.2f tr1=%.2f\n", k, shops[k].tr[0], shops[k].tr[1]);
            }
        }
        for (int k = 1; k <= K; ++k) if (shops[k].active) {
            for (int h = 0; h < 2; ++h) {
                shops[k].tr[h] += alpha * (shops[k].y[h] - shops[k].tr[h]);
            }
            shops[k].P[0] = priceF(shops[k].tr[0], shops[k].tr[1], overhead_f(shops[k].g[1]));
            shops[k].P[1] = priceF(shops[k].tr[1], shops[k].tr[0], overhead_f(shops[k].g[0]));
        }
    }

    // Research process for a prospective owner r
    // Should return an object ResearchResults, with targ0 and targ1 set as the local variables, and enter set as the return value of the function
    ResearchResults research(int r) {
        ResearchResults res;
        // trial targets
        res.targ0 = rng.uniform_int(xMax) + 1.0;
        res.targ1 = rng.uniform_int(xMax) + 1.0;

        double P0 = priceF(res.targ0, res.targ1, overhead_f(traders[r].d));
        double P1 = priceF(res.targ1, res.targ0, overhead_f(traders[r].s));

        // Test with a comrade
        int fr = comrade(traders[r]);
        Trader& partner = traders[fr];
        double Ucomp = u_sample(partner);
        double U = 0.0;
        // direct or indirect reachability check through fr’s links
        if (partner.d == traders[r].d) U = P0;
        else {
            int sh1 = partner.buy;
            if (sh1 > 0) {
                int m1 = (shops[sh1].g[0] == partner.d);
                if (shops[sh1].g[m1] == traders[r].d) U = P0 * shops[sh1].P[m1];
            }
        }
        // Test with a soulmate
        if (U < Ucomp) {
            U = 0;
            fr = soulmate(traders[r]);
            Trader& partner = traders[fr];
            // partner = traders[fr];
            Ucomp = u_sample(partner);
            U = 0.0;
            if (partner.s == traders[r].s) {
                U = P0;
            } else {
                int sh0 = partner.sell;
                if (sh0 > 0) {
                    int m0 = (shops[sh0].g[0] == partner.s);
                    if (shops[sh0].g[m0] == traders[r].s) U = shops[sh0].P[1 - m0] * P0;
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
            int k = 1 + rng.uniform_int(std::max(1, (int)consumes[traders[r].s].size()) - 1);
            int fr = consumes[traders[r].s][k - 1];
            Trader &partner = traders[fr];
            Ucomp = u_sample(partner);
            U = 0.0;
            if (partner.s == traders[r].d) {
                U = P1;
            } else {
                int sh0 = partner.sell;
                if (sh0 > 0) {
                    int m0 = (shops[sh0].g[0] == partner.s);
                    if (shops[sh0].g[m0] == traders[r].d) U = shops[sh0].P[1 - m0] * P1;
                }
            }
            // Stranger who produces d[r]
            if (U < Ucomp) {
                int k = 1 + rng.uniform_int(std::max(1, (int)produces[traders[r].d].size()) - 1);
                int  fr = produces[traders[r].d][k - 1];
                Trader& partner = traders[fr];
                Ucomp = u_sample(partner);
                U = 0.0;
                if (partner.d == traders[r].s) {
                    U = P1;
                } else {
                    int sh1 = partner.buy;
                    if (sh1 > 0) {
                        int m1 = (shops[sh1].g[0] == partner.d);
                        if (shops[sh1].g[m1] == traders[r].s) U = P1 * shops[sh1].P[m1];
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

    /*
     NOTE: ChatGPT says:
     Checking fr >= r is wrong: it skips any trader with id >= r (not just self) and can bias selection.
     The current index math (k+1) can go out of range when k is the last index.
     Simpler, correct approach: pick a random index in the full vector; if it equals r and size>1, pick the next (or wrap) — this guarantees a different trader without bias and avoids bounds errors.
    */
    int comrade(const Trader& trader) {
        // someone else producing s[r] (the same production good)
        int k = rng.uniform_int(std::max(0, std::max(1, (int)produces[trader.s].size()) - 1));
        int fr = produces[trader.s][k];
        if (fr >= trader.idx) fr = produces[trader.s][k+1]; // skip self
        return fr;
    }

    int soulmate(const Trader& trader) {
        // someone else consuming d[r] (the  same consumption good)
        int k = rng.uniform_int(std::max(0, std::max(1, (int)consumes[trader.d].size()) - 1));
        int fr = consumes[trader.d][k];
        if (fr >= trader.idx) fr = consumes[trader.d][k+1]; // skip self
        return fr;
    }

    double u_sample(const Trader& trader) {
        // attainable consumption for fr via current links
        double X = 0.0;
        Shop& sell_shop = shops[trader.sell];
        Shop& buy_shop = shops[trader.buy];
        int m0 = (sell_shop.g[0] == trader.s);
        int m1 = (buy_shop.g[0] == trader.d);
        global_m0 = m0;
        global_m1 = m1;
        if (DEBUG) {
            printf("[1] SET m1 to %d\n", global_m1);
        }
        if (sell_shop.g[m0] == trader.d) {
            X = sell_shop.P[1 - m0];
        } else {
            if (sell_shop.g[m0] == buy_shop.g[m1]) {
                X = sell_shop.P[1 - m0] * buy_shop.P[m1];
            }
        }

        return X;
    }

    double utility(const Trader& trader) {
        // attainable consumption for r via current links
        double X = 0.0;
        int a = trader.sell;
        int b = trader.buy;
        int m0 = (shops[a].g[0] == trader.s);
        int m1 = (shops[b].g[0] == trader.d);
        global_m0 = m0;
        global_m1 = m1;
        if (DEBUG) {
            printf("[2] SET m1 to %d\n", global_m1);
        }
        if (a > 0) {
            if (shops[a].g[m0] == trader.d) {
                X = shops[a].P[1 - m0];
            } else if (b > 0) {
                if (shops[a].g[m0] == shops[b].g[m1]) {
                    X = shops[a].P[1 - m0] * shops[b].P[m1];
                }
            }
        }
        return X;
    }

    void addshop(const Trader& trader, Shop& shop, std::vector<int>& cand) {
        if (shop.active && (shop.provides(trader.s) || shop.provides(trader.d))) {
            if (std::find(cand.begin(), cand.end(), shop.idx) == cand.end()) {
                cand.push_back(shop.idx);
            }
        }
    }

    void try_barter(int r, std::vector<int>& c, int& bestbarter, double& Ubarter, double& Ucomp) {
        bestbarter = 0;
        Ubarter = 0.0;
        // iterate candidates from index 2 onward
        for (size_t idx = 2; idx < c.size(); ++idx) {
            int k = c[idx];
            int qr = traders[r].q;
            if (shops[k].g[qr] == traders[r].s && shops[k].g[1 - qr] == traders[r].d) {
                double val = shops[k].P[qr];
                if (val > std::max(Ucomp, Ubarter)) {
                    bestbarter = k;
                    Ubarter = val;
                    c.erase(c.begin() + idx);
                    --idx;
                }
            }
        }
    }

    void try_one(int r, std::vector<int>& c, double& Ucomp) {
        int s = traders[r].s, d = traders[r].d;
        // Track which side matches for current c[0] and c[1]
        auto side_of = [&](int k, int good)->int { return shops[k].g[0] == good; };
        for (size_t idx = 2; idx < c.size(); ++idx) {
            Shop& shop = shops[c[idx]];
            // improve outlet (sell s)
            if (shop.provides(s)) {
                int ma = side_of(c[idx], s);
                if ((shop.g[ma] == shops[c[1]].g[global_m1]) || (shops[c[0]].P[1 - global_m0] == 0.0)) {
                    if (shops[c[0]].P[1 - global_m0] < shop.P[1 - ma]) {
                        double candidate = 0;
                        if (c[1] > 0) {
                            candidate = shop.P[1 - ma] * shops[c[1]].P[1 - global_m0];
                        }
                        c[0] = c[idx];
                        global_m0 = ma;
                        Ucomp = (shop.g[ma] == shops[c[1]].g[global_m1]) ? candidate : 0.0;
                        c.erase(c.begin() + idx);
                        --idx;
                    }
                }
            }
            else {
                // improve source (buy d)
                if (shop.provides(d)) {
                    int ma = (shop.g[0] == d);
                    if ((shop.g[ma] == shops[c[0]].g[global_m0]) || (shops[c[1]].P[global_m1] == 0.0)) {
                        if (shops[c[1]].P[global_m1] < shop.P[ma]) {
                            double candidate = (shops[c[0]].P[1 - global_m0]) * shop.P[ma];
                            c[1] = c[idx];
                            global_m1 = ma;
                            if (DEBUG) {
                                printf("[3] SET m1 to %d\n", global_m1);
                            }
                            Ucomp = (shop.g[ma] == shops[c[0]].g[global_m0]) ? candidate : 0.0;
                            c.erase(c.begin() + idx);
                            --idx;
                        }
                    }
                }
            }
        }
    }

    void try_two(int r, std::vector<int>& c, double& Ucomp) {
        int s = traders[r].s, d = traders[r].d;
        auto side_of = [&](int k, int good)->int { return shops[k].g[0] == good; };

        for (size_t ia = 2; ia < c.size(); ++ia) {
            int a = c[ia];
            if (shops[a].provides(s)) {
                for (size_t ib = 2; ib < c.size(); ++ib) {
                    if (ia == ib) continue;
                    int b = c[ib];
                    if (shops[b].provides(d)) {
                        int ma = side_of(a, s);
                        int mb = side_of(b, d);
                        // common intermediary condition
                        if (shops[a].g[ma] == shops[b].g[mb]) {
                            double val = shops[a].P[1 - ma] * shops[b].P[mb];
                            if (Ucomp < val) {
                                Ucomp = val;
                                global_m0 = ma;
                                global_m1 = mb;
                                if (DEBUG) {
                                    printf("[4] SET m1 to %d\n", global_m1);
                                }
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

        for (int r = 1; r <= m; ++r) {
            if (traders[r].sell <= 0) continue;
            int a = traders[r].sell;
            int b = traders[r].buy;
            int ma = (shops[a].g[0] == traders[r].s);
            int mb = (b > 0 && shops[b].g[0] == traders[r].d);

            part += (shops[a].g[ma] == traders[r].d) ||
                    (shops[a].g[ma] == shops[b].g[mb]);

            if (b > 0 && shops[a].g[ma] == shops[b].g[mb]) {
                moneytraders += 1.0;
                usingmoney[shops[b].g[mb]] += 1.0;
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
            if (endcount == 0) monyear = t / 50;
            endcount++;
        } else {
            endcount = 0;
        }

        return (endcount >= persist);
    }

    void calc2() {
        if (monetary == 0) monyear = -1;
        if (fulldev == 0) devyear = -1;

        // count non-money active shops
        BS = 0;
        for (int k = 1; k <= K; ++k) if (shops[k].active) {
            if (shops[k].g[0] != moneygood && shops[k].g[1] != moneygood) BS++;
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
        for (int r = 1; r <= m; ++r) Csurp += utility(traders[r]);

        Psurp = 0.0;
        Nshop = 0;
        for (int k = 1; k <= K; ++k) if (shops[k].active) {
            for (int h = 0; h < 2; ++h) {
                Psurp += (shops[k].y[h] - overhead_f(shops[k].g[h]) - shops[k].P[1 - h] * shops[k].y[1 - h]);
            }
            int own = shops[k].owner;
            int qown = (own > 0) ? traders[own].q : 0;
            if ((own > 0) && (shops[k].y[qown] > 1.0 || shops[k].y[1 - qown] > 0.0)) Nshop++;
        }
    }

    void rmse() {
        if (W > 0.0) {
            for (int h = 0; h < 2; ++h) {
                std::fill(vol[h].begin(), vol[h].end(), 0.0);
                std::fill(avp[h].begin(), avp[h].end(), 0.0);
            }

            // aggregate by pairs with moneygood
            for (int k = 1; k <= K; ++k) {
                if (!shops[k].active) continue;
                if (shops[k].g[0] == moneygood || shops[k].g[1] == moneygood) {
                    int ma = (shops[k].g[1] == moneygood);
                    int i = shops[k].g[1 - ma];
                    vol[0][i] += shops[k].y[1 - ma];
                    if (vol[0][i] > 0.0) {
                        avp[0][i] += (shops[k].y[1 - ma] / vol[0][i]) * (shops[k].P[1 - ma] - avp[0][i]);
                    }
                    vol[1][i] += shops[k].y[ma];
                    if (vol[1][i] > 0.0) {
                        avp[1][i] += (shops[k].y[ma] / vol[1][i]) * (shops[k].P[ma] - avp[1][i]);
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
