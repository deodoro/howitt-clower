#ifndef TESTRESEARCH_H
#define TESTRESEARCH_H

#include <vector>
#include <ResearchResults.h>
#include <Trader.h>
#include <Shop.h>
#include <PCG32.h>

extern PCG32 rng;

struct MatchEvaluation {
    double Ucomp = 0.0;
    double Ubarter = 0.0;
    Shop* barter = nullptr;
    Shop* candidate_seller = nullptr;
    Shop* candidate_buyer = nullptr;
};

class TestMatch {
public:
    TestMatch(std::vector<Trader>& traders, std::vector<Shop>& shops, std::vector<std::vector<int>>& produces, std::vector<std::vector<int>>& consumes, int m, int K, double f1, double slope, double xMax, double C);

    std::vector <MatchEvaluation> *weekly_matching(std::vector<Trader*> trader_line, std::vector<Shop> &shops);

    Trader& random_consumer(int good);
    Trader& random_producer(int good);
    Shop& random_shop();
    void addshop(const Trader* trader, Shop* shop, std::vector<int>& cand);
    void lineup();
    void try_barter(Trader& trader, std::vector<int>& c, MatchEvaluation& eval);
    void try_one(const Trader& trader, std::vector<int>& c, MatchEvaluation& eval);
    void try_two(const Trader& trader, std::vector<int>& c, MatchEvaluation& eval);

private:
    std::vector<Trader>& traders;
    std::vector<Shop>& shops;
    std::vector<std::vector<int>>& produces;
    std::vector<std::vector<int>>& consumes;
    std::vector<Trader*> trader_line;
    int m;
    int K;
    double f1;
    double slope;
    double xMax;
    double C;
    PCG32 rng;
};

#endif // TESTRESEARCH_H
