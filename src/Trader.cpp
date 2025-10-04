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
    // Returns a string representation of the trader's state for debugging.
    // Simulation rule: Each trader is defined by their supply and demand goods, shop links, and family shop status.
    return "Trader{s=" + std::to_string(get_supplied_good()) + ", d=" + std::to_string(get_supplied_good()) + ", q=" + std::to_string(q) +
           ", sh=[" + std::to_string(get_seller_idx()) + "," + std::to_string(get_buyer_idx()) + "], familyshop=" + std::to_string(get_familyshop()) + "}";
}

void Trader::sever_links(Shop& shop) {
    // Simulation rule: Traders sever links to shops that become inactive or unprofitable.
    if (familyshop &&  get_familyshop() == shop.idx) set_familyshop(0);
    if (seller_idx && get_seller_idx() == shop.idx) set_seller_idx(0);
    if (buyer_idx && get_buyer_idx() == shop.idx) set_buyer_idx(0);
}

/*
 NOTE: Check this algorithm. ChatGPT says:
 Checking fr >= r is wrong: it skips any trader with id >= r (not just self) and can bias selection.
 The current index math (k+1) can go out of range when k is the last index.
 Simpler, correct approach: pick a random index in the full vector; if it equals r and size>1, pick the next (or wrap) — this guarantees a different trader without bias and avoids bounds errors.
*/
int Trader::any_comrade(const std::vector<std::vector<int>>& produces) const {
    // Simulation rule: Selects another trader who produces the same good (potential direct competitor or trading partner).
    int k = rng.uniform_int(std::max(0, std::max(1, (int)produces[supplies].size()) - 1));
    int fr = produces[supplies][k];
    if (fr >= idx) fr = produces[supplies][k+1]; // skip self
    return fr;
}

// NOTE: should it consider comrades that have the same supply? They would have the same goods in both sides. Same applies to soulmates (in reverse)
int Trader::trade_comrade(const std::vector<std::vector<int>>& produces) const {
    // Simulation rule: Finds a comrade for potential trade, using the same supply good.
    int fr = any_comrade(produces);
    return fr;
}

// NOTE: Check this, same issue as comrades
int Trader::soulmate(const std::vector<std::vector<int>>& consumes) const {
    // Simulation rule: Selects another trader who consumes the same good (potential indirect trading partner).
    int k = rng.uniform_int(std::max(0, std::max(1, (int)consumes[demands].size()) - 1));
    int fr = consumes[demands][k];
    if (fr >= idx) fr = consumes[demands][k+1]; // skip self
    return fr;
}

double Trader::utility() const {
    // Simulation rule: Calculates the attainable consumption for this trader via current shop links.
    // Direct barter: If the seller shop offers the demanded good, utility is its price.
    // Indirect trade: If a buyer shop is linked and goods match, utility is the product of prices along the trade path.
    double X = 0.0;
    if (get_seller_idx() > 0) {
        Shop& sell_shop = shops_ref->at(get_seller_idx());
        if (sell_shop.get_the_other_good(supplies) == demands) {
            X = sell_shop.get_price_supply(supplies);
        } else {
            if (get_buyer_idx() > 0) {
                Shop& buy_shop = shops_ref->at(get_buyer_idx());
                if (sell_shop.get_the_other_good(supplies) == buy_shop.get_the_other_good(demands)) {
                    X = sell_shop.get_price_supply(supplies) * buy_shop.get_price_demand(demands);
                }
            }
        }
    }
    return X;
}

bool Trader::open_shop(Shop &shop)
{
    // Simulation rule: Trader opens a shop if it is inactive, sets goods and ownership, and links to self.
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
    // Simulation rule: Sets the supply good for the trader and updates the orientation (q) for shop assignment.
    supplies = s;
    q = (supplies > demands);
}

void Trader::set_demands(int d) {
    // Simulation rule: Sets the demand good for the trader and updates the orientation (q) for shop assignment.
    demands = d;
    q = (supplies > demands);
}

bool Trader::is_compatible_with(const Shop& shop) const {
    // Simulation rule: Trader is compatible with a shop if it is active and trades either the supply or demand good.
    return shop.is_active() && (shop.provides(supplies) || shop.provides(demands));
}

bool Trader::is_compatible_with(const Shop* shop) const {
    // Simulation rule: Trader is compatible with a shop pointer if it is active and trades either the supply or demand good.
    return shop->is_active() && (shop->provides(supplies) || shop->provides(demands));
}

bool Trader::allows_barter_with(const Shop& shop) const {
    // Simulation rule: Trader allows barter with a shop if it is active and trades both the supply and demand goods.
    return shop.is_active() && shop.provides(supplies) && shop.provides(demands);
}

bool Trader::wants_to_trade_in(int good) {
    // Simulation rule: Trader wants to trade in if the good matches their demand.
    return demands == good;
}

bool Trader::wants_to_trade_out(int good) {
    // Simulation rule: Trader wants to trade out if the good matches their supply.
    return supplies == good;
}

void Trader::set_buyer_shop(Shop* shop) {
    // Simulation rule: Sets the buyer shop pointer and index for the trader.
    buyer_shop = shop; buyer_idx = shop ? shop->idx : 0;
}

void Trader::set_seller_shop(Shop* shop) {
    // Simulation rule: Sets the seller shop pointer and index for the trader.
    seller_shop = shop; seller_idx = shop ? shop->idx : 0;
}
