#include "Trader.h"
#include "Shop.h"
#include "PCG32.h"
#include <algorithm>

extern PCG32 rng;

std::string Trader::to_string() const {
    return "Trader{s=" + std::to_string(supplies) + ", d=" + std::to_string(demands) + ", q=" + std::to_string(q) +
           ", sh=[" + std::to_string(seller_idx) + "," + std::to_string(buyer_idx) + "], familyshop=" + std::to_string(familyshop) + "}";
}

void Trader::sever_links(Shop& shop) {
    if (seller_idx == shop.idx) seller_idx = 0;
    if (buyer_idx == shop.idx) buyer_idx = 0;
}

/*
 NOTE: ChatGPT says:
 Checking fr >= r is wrong: it skips any trader with id >= r (not just self) and can bias selection.
 The current index math (k+1) can go out of range when k is the last index.
 Simpler, correct approach: pick a random index in the full vector; if it equals r and size>1, pick the next (or wrap) — this guarantees a different trader without bias and avoids bounds errors.
*/
int Trader::comrade(const std::vector<std::vector<int>>& produces) const {
    // someone else producing s[r] (the same production good)
    int k = rng.uniform_int(std::max(0, std::max(1, (int)produces[supplies].size()) - 1));
    int fr = produces[supplies][k];
    if (fr >= idx) fr = produces[supplies][k+1]; // skip self
    return fr;
}

int Trader::soulmate(const std::vector<std::vector<int>>& consumes) const {
    // someone else consuming d[r] (the  same consumption good)
    int k = rng.uniform_int(std::max(0, std::max(1, (int)consumes[demands].size()) - 1));
    int fr = consumes[demands][k];
    if (fr >= idx) fr = consumes[demands][k+1]; // skip self
    return fr;
}

double Trader::utility(std::vector<Shop>& shops) const {
    // attainable consumption for r via current links
    double X = 0.0;
    if (seller_idx > 0) {
        Shop& sell_shop = shops[seller_idx];
        if (sell_shop.get_good(supplies) == demands) {
            X = sell_shop.get_price(supplies, true);
        } else {
            if (buyer_idx > 0) {
                Shop& buy_shop = shops[buyer_idx];
                if (sell_shop.get_good(supplies) == buy_shop.get_good(demands)) {
                    X = sell_shop.get_price(supplies, true) * buy_shop.get_price(demands);
                }
            }
        }
    }
    return X;
}
