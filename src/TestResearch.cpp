#include <TestResearch.h>
#include <algorithm>

// Constructor
TestResearch::TestResearch(std::vector<Trader>& traders, std::vector<Shop>& shops, std::vector<std::vector<int>>& produces, std::vector<std::vector<int>>& consumes, int m, int K, double f1, double slope, double xMax, double C)
    : traders(traders), shops(shops), produces(produces), consumes(consumes), m(m), K(K), f1(f1), slope(slope), xMax(xMax), C(C) {}


std::function<double(int)> make_overhead(double f1, double slope) {
    return [f1, slope](int i) -> double {
        return f1 + (i - 1) * slope;
    };
}
// Research process for a prospective owner r
// Should return an object ResearchResults, with targ0 and targ1 set as the local variables, and enter set as the return value of the function
ResearchResults TestResearch::research(Trader& trader) {
    
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
        Since these are conditions for trading, TestResearch is considering
        non-matching roles as potential trade partners.
    */
    // Test with a comrade
    Trader& partner = this->traders[trader.any_comrade(this->produces)];
    double U = 0.0;
    double Ucomp = partner.utility(this->shops);
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
        Trader& partner = this->traders[trader.soulmate(this->consumes)];
        // partner = traders[fr];
        Ucomp = partner.utility(this->shops);
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
        Ucomp = partner.utility(this->shops);
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
            Ucomp = partner.utility(this->shops);
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

Trader& TestResearch::random_consumer(int good) {
    int k = 1 + rng.uniform_int(std::max(1, (int)this->consumes[good].size()) - 1);
    return this->traders[this->consumes[good][k - 1]];
}

Trader& TestResearch::random_producer(int good) {
    int k = 1 + rng.uniform_int(std::max(1, (int)this->produces[good].size()) - 1);
    return this->traders[this->produces[good][k - 1]];
}

Shop& TestResearch::random_shop() {
    return this->shops[rng.uniform_int(this->K) + 1];
}

void TestResearch::addshop(const Trader* trader, Shop* shop, std::vector<int>& cand) {
    if (shop && shop->active && trader->is_compatible_with(shop)) {
        if (std::find(cand.begin(), cand.end(), shop->idx) == cand.end()) {
            if (!(shop->idx != trader->get_buyer_idx() && shop->idx != trader->get_seller_idx())) {
                printf("hit!\n");
            }
            cand.push_back(shop->idx);
        }
    }
}

void TestResearch::lineup() {
    // random permutation of 1..m
    int line[this->m + 1];
    int start[this->m + 1];
    for (int i = 1; i <= this->m; ++i)
        start[i] = i;
    for (int j = 1; j <= this->m; ++j) {
        int k = rng.uniform_int(this->m - j + 1) + 1;
        line[j] = start[k];
        for (int i = k; i <= this->m - j; ++i) {
            start[i] = start[i + 1];
        }
    }
    this->trader_line.clear();
    for (int i = 0; i < this->m; ++i) {
        this->trader_line.push_back(&this->traders[line[i + 1]]);
    }
}

void TestResearch::try_barter(Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    // iterate candidates from index 2 onward
    for (size_t idx = 2; idx < c.size(); ++idx) {
        Shop& shop = this->shops[c[idx]];
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

void TestResearch::try_one(const Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    int s = trader.get_supplies();
    int d = trader.get_demands();
    // Track which side matches for current c[0] and c[1]
    for (size_t idx = 2; idx < c.size(); ++idx) {
        Shop& shop = this->shops[c[idx]];
        Shop& candidate_0 = this->shops[c[0]];
        Shop& candidate_1 = this->shops[c[1]];
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

void TestResearch::try_two(const Trader& trader, std::vector<int>& c, struct MatchEvaluation& eval) {
    for (size_t ia = 2; ia < c.size(); ++ia) {
        int a = c[ia];
        if (this->shops[a].provides(trader.get_supplies())) {
            for (size_t ib = 2; ib < c.size(); ++ib) {
                if (ia == ib) continue;
                int b = c[ib];
                if (this->shops[b].provides(trader.get_demands())) {
                    // common intermediary condition
                    if (this->shops[a].get_good(trader.get_supplies()) == this->shops[b].get_good(trader.get_demands())) {
                        double val = this->shops[a].get_price_supply(trader.get_supplies()) * this->shops[b].get_price(trader.get_demands());
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
