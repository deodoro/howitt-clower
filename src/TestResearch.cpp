#include "TestResearch.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstddef>
#include "PCG32.h"

PCG32 rng;

// Constructor
TestResearch::TestResearch(std::vector<Trader> &traders, std::vector<Shop> &shops, std::vector<std::vector<int>> &produces, std::vector<std::vector<int>> &consumes, int m, int K, double f1, double slope, double xMax, double C)
    : traders(traders), shops(shops), produces(produces), consumes(consumes), m(m), K(K), f1(f1), slope(slope), xMax(xMax), C(C)
{
    ::rng.seed(12345, 67890);
    this->rng.seed(12345, 67890);
}

std::function<double(int)> make_overhead(double f1, double slope)
{
    return [f1, slope](int i) -> double
    {
        return f1 + (i - 1) * slope;
    };
}
// Research process for a prospective owner r
// Should return an object ResearchResults, with targ0 and targ1 set as the local variables, and enter set as the return value of the function
ResearchResults TestResearch::research(Trader &trader)
{

    auto overhead_f = make_overhead(f1, slope);
    auto priceF = [this](double tr0, double tr1, double f_other)
    {
        return (tr1 - f_other - C > 0.0) ? ((tr1 - C - f_other) / tr0) : 0.0;
    };

    // Prepare targets and compute prices
    ResearchResults res{rng.uniform_int(xMax) + 1.0, rng.uniform_int(xMax) + 1.0, 0};
    double P0 = priceF(res.targ0, res.targ1, overhead_f(trader.get_supplied_good()));
    double P1 = priceF(res.targ1, res.targ0, overhead_f(trader.get_supplied_good()));

    /*
        NOTE: Should be checking if comrade is same buyer or seller shop.
        Since these are conditions for trading, TestResearch is considering
        non-matching roles as potential trade partners.
    */
    // Test with a comrade
    Trader &partner = this->traders[trader.any_comrade(this->produces)];
    double U = 0.0;
    double Ucomp = partner.utility(this->shops);
    // direct or indirect reachability check through fr's links
    if (trader.wants_to_trade_in(partner.get_supplied_good()))
    {
        U = P0;
    }
    else
    {
        if (partner.get_source() && trader.wants_to_trade_in(partner.get_source()->get_the_other_good(partner.get_supplied_good())))
        {
            U = partner.get_source()->get_price_demand(partner.get_supplied_good()) * P0;
        }
    }
    // Test with a soulmate
    if (U < Ucomp)
    {
        Trader &partner = this->traders[trader.soulmate(this->consumes)];
        // partner = traders[fr];
        Ucomp = partner.utility(this->shops);
        U = 0.0;
        if (trader.wants_to_trade_out(partner.get_supplied_good()))
        {
            U = P0;
        }
        else
        {
            if (partner.get_outlet() && trader.wants_to_trade_out(partner.get_outlet()->get_the_other_good(partner.get_supplied_good())))
            {
                U = partner.get_outlet()->get_price_supply(partner.get_supplied_good()) * P0;
            }
        }
        if (U < Ucomp)
        {
            res.enter = 0;
            return res;
        }
    }

    {
        U = 0.0;
        Trader &partner = random_consumer(trader.get_supplied_good());
        Ucomp = partner.utility(this->shops);
        if (trader.wants_to_trade_in(partner.get_supplied_good()))
        {
            U = P1;
        }
        else
        {
            if (partner.get_outlet() && trader.wants_to_trade_in(partner.get_outlet()->get_the_other_good(partner.get_supplied_good())))
            {
                U = partner.get_outlet()->get_price_demand(trader.get_supplied_good()) * P1;
            }
        }
        // Stranger who produces d[r]
        if (U < Ucomp)
        {
            Trader &partner = random_producer(trader.get_supplied_good());
            Ucomp = partner.utility(this->shops);
            U = 0.0;
            if (trader.wants_to_trade_out(partner.get_supplied_good()))
            {
                U = P1;
            }
            else
            {
                if (partner.get_source() && trader.wants_to_trade_out(partner.get_source()->get_the_other_good(partner.get_supplied_good())))
                {
                    U = P1 * partner.get_source()->get_price_demand(partner.get_supplied_good());
                }
            }
            if (U < Ucomp)
            {
                res.enter = 0;
                return res;
            }
        }
    }
    res.enter = 1;
    return res;
}

Trader &TestResearch::random_consumer(int good)
{
    auto &vec = this->consumes[good];
    int idx = rng.uniform_int(vec.size());
    return this->traders[vec[idx]];
}

Trader &TestResearch::random_producer(int good)
{
    auto &vec = this->produces[good];
    int idx = rng.uniform_int(vec.size());
    return this->traders[vec[idx]];
}

Shop &TestResearch::random_shop()
{
    return this->shops[rng.uniform_int(this->K) + 1];
}

void TestResearch::addshop(const Trader *trader, Shop *shop, std::vector<int> &cand)
{
    if (shop && shop->active && trader->is_compatible_with(shop))
    {
        if (std::find(cand.begin(), cand.end(), shop->idx) == cand.end())
        {
            if (!(shop->idx != trader->get_source_idx() && shop->idx != trader->get_outlet_idx()))
            {
                printf("hit!\n");
            }
            cand.push_back(shop->idx);
        }
    }
}

void TestResearch::lineup()
{
    // random permutation of 1..m
    int line[this->m + 1];
    int start[this->m + 1];
    for (int i = 1; i <= this->m; ++i)
        start[i] = i;
    for (int j = 1; j <= this->m; ++j)
    {
        int k = rng.uniform_int(this->m - j + 1) + 1;
        line[j] = start[k];
        for (int i = k; i <= this->m - j; ++i)
        {
            start[i] = start[i + 1];
        }
    }
    this->trader_line.clear();
    for (int i = 0; i < this->m; ++i)
    {
        this->trader_line.push_back(&this->traders[line[i + 1]]);
    }
}

void TestResearch::try_barter(Trader &trader, std::vector<int> &c, struct MatchEvaluation &eval)
{
    // iterate candidates from index 2 onward
    for (size_t idx = 2; idx < c.size(); ++idx)
    {
        Shop &shop = this->shops[c[idx]];
        if (trader.allows_barter_with(shop))
        {
            double val = shop.get_price_supply(trader.get_supplied_good());
            if (val > std::max(eval.Ucomp, eval.Ubarter))
            {
                eval.barter = &shop;
                eval.Ubarter = val;
                c.erase(c.begin() + idx);
                --idx;
            }
        }
    }
}

void TestResearch::try_one(const Trader &trader, std::vector<int> &c, struct MatchEvaluation &eval)
{
    int s = trader.get_supplied_good();
    int d = trader.get_supplied_good();
    // Track which side matches for current c[0] and c[1]
    for (size_t idx = 2; idx < c.size(); ++idx)
    {
        Shop &shop = this->shops[c[idx]];
        Shop &candidate_0 = this->shops[c[0]];
        Shop &candidate_1 = this->shops[c[1]];
        // improve outlet (sell s)
        if (shop.provides(s))
        {
            if ((shop.get_the_other_good(s) == candidate_1.get_the_other_good(d)) || (candidate_0.get_price_supply(s) == 0.0))
            {
                if (candidate_0.get_price_supply(s) < shop.get_price_supply(s))
                {
                    if (shop.get_the_other_good(s) == candidate_1.get_the_other_good(d))
                    {
                        eval.Ucomp = shop.get_price_supply(s) * candidate_1.get_price_supply(s);
                    }
                    c[0] = c[idx];
                    eval.candidate_seller = &shop;
                    c.erase(c.begin() + idx);
                    --idx;
                }
            }
        }
        else
        {
            // improve source (buy d)
            if (shop.provides(d))
            {
                if ((shop.get_the_other_good(d) == candidate_0.get_the_other_good(s)) || (candidate_1.get_price_demand(d) == 0.0))
                {
                    if (candidate_1.get_price_demand(d) < shop.get_price_demand(d))
                    {
                        if (shop.get_the_other_good(d) == candidate_0.get_the_other_good(s))
                        {
                            eval.Ucomp = (candidate_0.get_price_supply(s) * shop.get_price_demand(d));
                        }
                        c[1] = c[idx];
                        eval.candidate_buyer = &shop;
                        c.erase(c.begin() + idx);
                        --idx;
                    }
                }
            }
        }
    }
}

void TestResearch::try_two(const Trader &trader, std::vector<int> &c, struct MatchEvaluation &eval)
{
    for (size_t ia = 2; ia < c.size(); ++ia)
    {
        int a = c[ia];
        if (this->shops[a].provides(trader.get_supplied_good()))
        {
            for (size_t ib = 2; ib < c.size(); ++ib)
            {
                if (ia == ib)
                    continue;
                int b = c[ib];
                if (this->shops[b].provides(trader.get_supplied_good()))
                {
                    // common intermediary condition
                    if (this->shops[a].get_the_other_good(trader.get_supplied_good()) == this->shops[b].get_the_other_good(trader.get_supplied_good()))
                    {
                        double val = this->shops[a].get_price_supply(trader.get_supplied_good()) * this->shops[b].get_price_demand(trader.get_supplied_good());
                        if (eval.Ucomp < val)
                        {
                            eval.Ucomp = val;
                            eval.candidate_seller = &shops[a];
                            eval.candidate_buyer = &shops[b];
                            c[0] = a;
                            c[1] = b;
                        }
                    }
                }
            }
        }
    }
}

int main()
{
    std::vector<Trader> traders(7);
    std::vector<Shop> shops(5);

    // Map of producers/consumers by good id (size = max_good_id + 1)
    std::vector<std::vector<int>> produces(4); // 0..3
    std::vector<std::vector<int>> consumes(4); // 0..3

    // --- Trader endowments & desires (covering all cross-good pairs) ---
    // t1: 1->2, t2: 2->1 (comrade pair)
    // t3: 3->1, t4: 1->3 (soulmate pair)
    // t5: 2->3, t6: 3->2 (random consumer/producer branches)
    traders[1].set_supplies(1);
    traders[1].set_demands(2);
    traders[2].set_supplies(2);
    traders[2].set_demands(1);
    traders[3].set_supplies(3);
    traders[3].set_demands(1);
    traders[4].set_supplies(1);
    traders[4].set_demands(3);
    traders[5].set_supplies(2);
    traders[5].set_demands(3);
    traders[6].set_supplies(3);
    traders[6].set_demands(2);

    // Provide shop pointers to traders
    for (int i = 1; i <= 6; ++i)
        traders[i].set_shops(&shops);

    // --- Shops: provide various goods and price configurations ---
    // Each shop has g[0] = "supply side counter-good" and g[1] = "demand side counter-good"
    // and P[0]/P[1] their corresponding prices. We set them to create both win/lose cases.
    shops[1].idx = 1;
    shops[1].active = 1;
    shops[1].g[0] = 1;
    shops[1].g[1] = 2;
    shops[1].P[0] = 1.20;
    shops[1].P[1] = 1.00;

    shops[2].idx = 2;
    shops[2].active = 1;
    shops[2].g[0] = 2;
    shops[2].g[1] = 1;
    shops[2].P[0] = 1.50;
    shops[2].P[1] = 1.30;

    shops[3].idx = 3;
    shops[3].active = 1;
    shops[3].g[0] = 3;
    shops[3].g[1] = 1;
    shops[3].P[0] = 0.90;
    shops[3].P[1] = 1.60;

    shops[4].idx = 4;
    shops[4].active = 0; // inactive shop to ensure filtering works
    shops[4].g[0] = 2;
    shops[4].g[1] = 3;
    shops[4].P[0] = 2.00;
    shops[4].P[1] = 2.00;

    // --- Produces/Consumes adjacency (for any_comrade/soulmate & random_* lookups) ---
    // Producers by good
    produces[1] = {1, 4}; // traders producing good 1
    produces[2] = {2, 5}; // traders producing good 2
    produces[3] = {3, 6}; // traders producing good 3

    // Consumers by good (who demand that good)
    consumes[1] = {2, 3}; // need 1
    consumes[2] = {1, 6}; // need 2
    consumes[3] = {4, 5}; // need 3

    // --- TestResearch with parameters to exercise overhead & entry calculus ---
    // m = number of active trader slots (use 6), K = 4 shops, f1 & slope control overhead,
    // xMax bounds random target draws, C is fixed cost.
    int m = 6, K = 4;
    double f1 = 0.05, slope = 0.25, xMax = 10.0, C = 0.40;

    TestResearch test(traders, shops, produces, consumes, m, K, f1, slope, xMax, C);

    // --- Execute research across a suite of traders to hit diverse branches ---
    // Expectation:
    //  - t1 vs t2 exercises "comrade" branch
    //  - t4 vs t3 exercises "soulmate" branch
    //  - t5/t6 exercise random consumer/producer fallbacks
    //  - inactive shop(4) ensures compatibility filter
    for (int tid : {1, 2, 3, 4, 5, 6})
    {
        ResearchResults res = test.research(traders[tid]);
        std::cout << "[t" << tid << "] enter=" << res.enter
                  << " targ0=" << res.targ0
                  << " targ1=" << res.targ1 << "\n";
    }

    return 0;
}
