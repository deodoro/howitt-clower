#include "PCG32.h"
#include "TestMatch.h"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <stack>
#include <vector>

PCG32 rng;
static constexpr double lambda = 0.05;

// Constructor
TestMatch::TestMatch(std::vector<Trader>& traders, std::vector<Shop>& shops,
	std::vector<std::vector<int>>& produces,
	std::vector<std::vector<int>>& consumes, int m, int K,
	double f1, double slope, double xMax, double C)
	: traders(traders), shops(shops), produces(produces), consumes(consumes),
	m(m), K(K), f1(f1), slope(slope), xMax(xMax), C(C) {
	::rng.seed(12345, 67890);
	this->rng.seed(12345, 67890);
}

std::function<double(int)> make_overhead(double f1, double slope) {
	return [f1, slope](int i) -> double { return f1 + (i - 1) * slope; };
}

// Matching: agents sample a small set of shops and adopt best links
std::vector<MatchEvaluation>*
TestMatch::weekly_matching(std::vector<Trader*> trader_line,
	std::vector<Shop>& shops) {
	std::vector<MatchEvaluation>* response = new std::vector<MatchEvaluation>();

	for (Trader* trader_p : trader_line) {
		Trader& trader = *trader_p;
		double U = trader.utility(shops);
		double psearch = (U > 0.0 ? lambda : 1.0);
		// Skip condition: random or already owns a shop
		if (rng.uniform01_inclusive() < psearch &&
			trader.get_familyshop() == nullptr) {
			struct MatchEvaluation eval;
			eval.Ucomp = U;
			eval.candidate_seller = trader.get_outlet();
			eval.candidate_buyer = trader.get_source();
			// candidate initialization with current links
			std::stack<Shop*> cand_stack;
			cand_stack.push(trader.get_outlet());
			cand_stack.push(trader.get_source());

			std::vector<int> cand;
			cand.reserve(8);
			cand.push_back(trader.get_outlet_idx()); // c[0]
			cand.push_back(trader.get_source_idx());  // c[1]

			// add friend outlets/sources and one random shop
			int comrade_idx = trader.trade_comrade(produces);
			Trader& comrade_ = traders[comrade_idx];
			/*
			NOTE: I added comrade_ as a friend outlet by mistake, and it worked. Why?
			*/
			// if (comrade_idx != 0)
			addshop(&trader, comrade_.get_outlet(), cand);

			Trader& soulmate_ = traders[trader.soulmate(consumes)];
			addshop(&trader, soulmate_.get_source(), cand);

			// Shop* random = &random_shop();
			Shop* random = &shops[1];
			addshop(&trader, random, cand);

			if (cand.size() > 2) {
				try_barter(trader, cand, eval);
				if ((shops[cand[0]].get_the_other_good(trader.get_supplied_good()) !=
					trader.get_supplied_good() ||
					shops[cand[0]].get_price_supply(trader.get_supplied_good()) ==
					0.0)) {
					try_one(trader, cand, eval);
				}
				try_two(trader, cand, eval);

				if (eval.Ucomp < eval.Ubarter && eval.barter != nullptr) {
					trader.set_source(nullptr);
					trader.set_outlet(eval.barter);
				}
				else {
					// adopt c[0], c[1] as improved chain if any
					trader.set_outlet(eval.candidate_seller);
					trader.set_source(eval.candidate_buyer);
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
	}
	return response;
}

Trader& TestMatch::random_consumer(int good) {
	auto& vec = this->consumes[good];
	int idx = rng.uniform_int(vec.size());
	return this->traders[vec[idx]];
}

Trader& TestMatch::random_producer(int good) {
	auto& vec = this->produces[good];
	int idx = rng.uniform_int(vec.size());
	return this->traders[vec[idx]];
}

Shop& TestMatch::random_shop() {
	return this->shops[rng.uniform_int(this->K) + 1];
}

void TestMatch::addshop(const Trader* trader, Shop* shop,
	std::vector<int>& cand) {
	if (shop && shop->active && trader->is_compatible_with(shop)) {
		if (std::find(cand.begin(), cand.end(), shop->idx) == cand.end()) {
			if (!(shop->idx != trader->get_source_idx() &&
				shop->idx != trader->get_outlet_idx())) {
				printf("hit!\n");
			}
			cand.push_back(shop->idx);
		}
	}
}

void TestMatch::lineup() {
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

void TestMatch::try_barter(Trader& trader, std::vector<int>& c,
	struct MatchEvaluation& eval) {
	// iterate candidates from index 2 onward
	for (size_t idx = 2; idx < c.size(); ++idx) {
		Shop& shop = this->shops[c[idx]];
		// Does the shop provides both supply and demand good?
		// Shouldn't it be a real match of demand against supply?
		if (trader.allows_barter_with(shop)) {
			double val = shop.get_price_supply(trader.get_supplied_good());
			if (val > std::max(eval.Ucomp, eval.Ubarter)) {
				eval.barter = &shop;
				eval.Ubarter = val;
				c.erase(c.begin() + idx);
				--idx;
			}
		}
	}
}

void TestMatch::try_one(const Trader& trader, std::vector<int>& c,
	struct MatchEvaluation& eval) {
    int s = trader.get_supplied_good();
    int d = trader.get_demand_good();
    // Track which side matches for current c[0] and c[1]
    for (size_t idx = 2; idx < c.size(); ++idx) {
        Shop& shop = shops[c[idx]];
        Shop& candidate_0 = shops[c[0]];
        Shop& candidate_1 = shops[c[1]];
        // improve outlet (sell s)
        if (shop.provides(s)) {
            if ((shop.get_the_other_good(s) == candidate_1.get_the_other_good(d)) || (candidate_0.get_price_supply(s) == 0.0)) {
                if (candidate_0.get_price_supply(s) < shop.get_price_supply(s)) {
                    if (shop.get_the_other_good(s) == candidate_1.get_the_other_good(d)) {
                        eval.Ucomp = shop.get_price_supply(s) * candidate_1.get_price_supply(s);
                    }
                    c[0] = c[idx];
                    eval.candidate_seller = &shop;
                    c.erase(c.begin() + idx);
                    --idx;
                }
            }
        }
        else {
            // improve source (buy d)
            if (shop.provides(d)) {
                if ((shop.get_the_other_good(d) == candidate_0.get_the_other_good(s)) || (candidate_1.get_price_demand(d) == 0.0)) {
                    if (candidate_1.get_price_demand(d) < shop.get_price_demand(d)) {
                        if (shop.get_the_other_good(d) == candidate_0.get_the_other_good(s)) {
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

void TestMatch::try_two(const Trader& trader, std::vector<int>& c,
	struct MatchEvaluation& eval) {
	for (size_t ia = 2; ia < c.size(); ++ia) {
		int a = c[ia];
		if (this->shops[a].provides(trader.get_supplied_good())) {
			for (size_t ib = 2; ib < c.size(); ++ib) {
				if (ia == ib)
					continue;
				int b = c[ib];
				if (this->shops[b].provides(trader.get_supplied_good())) {
					// common intermediary condition
					if (this->shops[a].get_the_other_good(trader.get_supplied_good()) ==
						this->shops[b].get_the_other_good(trader.get_supplied_good())) {
						double val =
							this->shops[a].get_price_supply(trader.get_supplied_good()) *
							this->shops[b].get_price_demand(trader.get_supplied_good());
						if (eval.Ucomp < val) {
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

int main() {
	std::vector<Trader> traders(7);
	std::vector<Shop> shops(5);

	// Map of producers/consumers by good id (size = max_good_id + 1)
	std::vector<std::vector<int>> produces(4); // 0..3
	std::vector<std::vector<int>> consumes(4); // 0..3

	traders[1].set_supplies(1);
	traders[1].set_demands(2);
	traders[2].set_supplies(1);
	traders[2].set_demands(3);
	traders[3].q = 1;
	traders[3].set_supplies(2);
	traders[3].set_demands(1);
	traders[4].set_supplies(2);
	traders[4].set_demands(3);
	traders[5].q = 1;
	traders[5].set_supplies(3);
	traders[5].set_demands(1);
	traders[6].q = 1;
	traders[6].set_supplies(3);
	traders[6].set_demands(2);

	// Trader 1: q = 0, s = 1, d = 2
	// Trader 2: q = 0, s = 1, d = 3
	// Trader 3: q = 1, s = 2, d = 1
	// Trader 4: q = 0, s = 2, d = 3
	// Trader 5: q = 1, s = 3, d = 1
	// Trader 6: q = 1, s = 3, d = 2
	// Shop 1: owner = 3,  g = [1,2], P = [.54, 1.26], tr = [100, 75]

	// Provide shop pointers to traders
	for (int i = 1; i <= 6; ++i) {
		traders[i].idx = i;
		traders[i].set_traders(&traders);
		traders[i].set_shops(&shops);
	}

	// --- Shops: provide various goods and price configurations ---
	// Each shop has g[0] = "supply side counter-good" and g[1] = "demand side
	// counter-good" and P[0]/P[1] their corresponding prices. We set them to
	// create both win/lose cases.
	shops[1].idx = 1;
	shops[1].active = 1;
	shops[1].g[0] = 1;
	shops[1].g[1] = 2;
	shops[1].P[0] = 1.20;
	shops[1].P[1] = 1.00;

	for (int i = 2; i <= 4; ++i) {
		shops[i].idx = i;
	}

	// --- Produces/Consumes adjacency (for any_comrade/soulmate & random_*
	// lookups) --- Producers by good
	produces[1] = { 1, 2 }; // traders producing good 1
	produces[2] = { 3, 4 }; // traders producing good 2
	produces[3] = { 5, 6 }; // traders producing good 3

	// Consumers by good (who demand that good)
	consumes[1] = { 3, 5 }; // need 1
	consumes[2] = { 1, 6 }; // need 2
	consumes[3] = { 2, 4 }; // need 3

	// --- TestMatch with parameters to exercise overhead & entry calculus ---
	// m = number of active trader slots (use 6), K = 4 shops, f1 & slope control
	// overhead, xMax bounds random target draws, C is fixed cost.
	int m = 6, K = 4;
	double f1 = 0.00, slope = 0.25, xMax = 200.0, C = 5.00;

	TestMatch test(traders, shops, produces, consumes, m, K, f1, slope, xMax, C);

	// --- Execute research across a suite of traders to hit diverse branches ---
	// Expectation:
	//  - t1 vs t2 exercises "comrade" branch
	//  - t4 vs t3 exercises "soulmate" branch
	//  - t5/t6 exercise random consumer/producer fallbacks
	//  - inactive shop(4) ensures compatibility filter
	std::vector<Trader*> trader_line{ &traders[1], &traders[2], &traders[3],
																		&traders[4], &traders[5], &traders[6] };
	std::vector<MatchEvaluation>* res = test.weekly_matching(trader_line, shops);

	int i = 0;
	for (const auto& eval : *res) {
		std::cout << "[t" << i++ << "].Ucomp=" << eval.Ucomp << "\n";
	}

	return 0;
}

//
// Round 1: seller_shop =1, buyer_shop=0
// Round 1: seller_shop =1, buyer_shop=0
// Round 1: barter =1, buyer_shop=0

// Trader 1: q = 0, s = 1, d = 2, idx = 0
// Trader 2: q = 0, s = 1, d = 3, idx = 0
// Trader 3: q = 1, s = 2, d = 1, idx = 0
// Trader 4: q = 0, s = 2, d = 3, idx = 0
// Trader 5: q = 1, s = 3, d = 1, idx = 0
// Trader 6: q = 1, s = 3, d = 2, idx = 0
// Shop 1: owner = 3,  g = [1,2], P = [.54, 1.26], tr = [100, 75]
