#ifndef TRADER_H
#define TRADER_H

#include <string>
#include <vector>

// Forward declarations
class Shop;
extern class PCG32 rng;

class Trader {
public:
    int idx{0};             // PROVISIONAL: to make the code compatible while it's refactored
    int supplies{0};        // produced good
    int demands{0};         // desired good
    int q{0};               // 0 if s<=d else 1 (orders shop's good pair)
    int seller_idx{0};
    int buyer_idx{0};
    int familyshop{0};      // owned shop index (0 if none)

    std::string to_string() const;
    void sever_links(Shop& shop);
    int comrade(const std::vector<std::vector<int>>& produces) const;
    int soulmate(const std::vector<std::vector<int>>& consumes) const;
    double utility(std::vector<Shop>& shops) const;
};

#endif // TRADER_H
