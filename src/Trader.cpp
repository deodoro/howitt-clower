/*
 * Trader.cpp
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

#include "Trader.h"
#include "Shop.h"
#include "PCG32.h"
#include <algorithm>

extern PCG32 rng;

std::string Trader::to_string() const {
    return "Trader{s=" + std::to_string(get_supplies()) + ", d=" + std::to_string(get_demands()) + ", q=" + std::to_string(q) +
           ", sh=[" + std::to_string(get_seller_idx()) + "," + std::to_string(get_buyer_idx()) + "], familyshop=" + std::to_string(get_familyshop()) + "}";
}

void Trader::sever_links(Shop& shop) {
    if (familyshop &&  get_familyshop() == shop.idx) set_familyshop(0);
    if (seller_idx && get_seller_idx() == shop.idx) set_seller_idx(0);
    if (buyer_idx && get_buyer_idx() == shop.idx) set_buyer_idx(0);
}

/*
 NOTE: ChatGPT says:
 Checking fr >= r is wrong: it skips any trader with id >= r (not just self) and can bias selection.
 The current index math (k+1) can go out of range when k is the last index.
 Simpler, correct approach: pick a random index in the full vector; if it equals r and size>1, pick the next (or wrap) — this guarantees a different trader without bias and avoids bounds errors.
*/
int Trader::any_comrade(const std::vector<std::vector<int>>& produces) const {
    // someone else producing s[r] (the same production good)
    int k = rng.uniform_int(std::max(0, std::max(1, (int)produces[supplies].size()) - 1));
    int fr = produces[supplies][k];
    if (fr >= idx) fr = produces[supplies][k+1]; // skip self
    return fr;
}

int Trader::trade_comrade(const std::vector<std::vector<int>>& produces) const {
    int fr = any_comrade(produces);
    // if (fr == buyer_idx || fr == seller_idx) return 0;
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
    if (get_seller_idx() > 0) {
        Shop& sell_shop = shops[get_seller_idx()];
        if (sell_shop.get_good(supplies) == demands) {
            X = sell_shop.get_price_supply(supplies);
        } else {
            if (get_buyer_idx() > 0) {
                Shop& buy_shop = shops[get_buyer_idx()];
                if (sell_shop.get_good(supplies) == buy_shop.get_good(demands)) {
                    X = sell_shop.get_price_supply(supplies) * buy_shop.get_price(demands);
                }
            }
        }
    }
    return X;
}

bool Trader::open_shop(Shop &shop)
{
    if (shop.active) {
        return false;
    }
    else {
        shop.active = 1;
        shop.g[1 - q] = demands;
        shop.g[q]     = supplies;
        // owner links
        set_seller_idx(shop.idx);
        set_buyer_idx(0);
        set_familyshop(shop.idx);
        shop.owner = idx;
        return true;
    }
}

void Trader::set_supplies(int s) {
    supplies = s;
    q = (supplies > demands);
}

void Trader::set_demands(int d) {
    demands = d;
    q = (supplies > demands);
}

bool Trader::is_compatible_with(const Shop& shop) const {
    return shop.is_active() && (shop.provides(supplies) || shop.provides(demands));
}

bool Trader::is_compatible_with(const Shop* shop) const {
    return shop->is_active() && (shop->provides(supplies) || shop->provides(demands));
}

bool Trader::allows_barter_with(const Shop& shop) const {
    return shop.is_active() && shop.provides(supplies) && shop.provides(demands);
}

bool Trader::wants_to_trade_in(int good) {
    return demands == good;
}

bool Trader::wants_to_trade_out(int good) {
    return supplies == good;
}

void Trader::set_buyer_shop(Shop* shop) { 
    buyer_shop = shop; buyer_idx = shop ? shop->idx : 0; 
}

void Trader::set_seller_shop(Shop* shop) { 
    seller_shop = shop; seller_idx = shop ? shop->idx : 0; 
}
