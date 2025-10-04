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
           ", sh=[" + std::to_string(get_outlet_idx()) + "," + std::to_string(get_source_idx()) + "], familyshop=" + std::to_string(get_familyshop_idx()) + "}";
}

void Trader::sever_links(Shop& shop) {
    // Simulation rule: Traders sever links to shops that become inactive or unprofitable.
    if (familyshop &&  get_familyshop_idx() == shop.idx) set_familyshop_idx(0);
    if (outlet_idx && get_outlet_idx() == shop.idx) set_outlet_idx(0);
    if (source_idx && get_source_idx() == shop.idx) set_source_idx(0);
}

/*
 NOTE: Check this algorithm. ChatGPT says:
 Checking fr >= r is wrong: it skips any trader with id >= r (not just self) and can bias selection.
 The current index math (k+1) can go out of range when k is the last index.
 Simpler, correct approach: pick a random index in the full vector; if it equals r and size>1, pick the next (or wrap) — this guarantees a different trader without bias and avoids bounds errors.
*/
int Trader::any_comrade(const std::vector<std::vector<int>>& produces) const {
    // Simulation rule: Selects another trader who produces the same good (potential direct competitor or trading partner).
    int k = rng.uniform_int(std::max(0, std::max(1, (int)produces[supplied_good].size()) - 1));
    int fr = produces[supplied_good][k];
    if (fr >= idx) fr = produces[supplied_good][k+1]; // skip self
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
    int k = rng.uniform_int(std::max(0, std::max(1, (int)consumes[demanded_good].size()) - 1));
    int fr = consumes[demanded_good][k];
    if (fr >= idx) fr = consumes[demanded_good][k+1]; // skip self
    return fr;
}

double Trader::utility() const {
    // Simulation rule: Calculates the attainable consumption for this trader via current shop links.
    // Direct barter: If the seller shop offers the demanded good, utility is its price.
    // Indirect trade: If a buyer shop is linked and goods match, utility is the product of prices along the trade path.
    double X = 0.0;
    if (get_outlet() != nullptr) {
        Shop* sell_shop = get_outlet();
        if (sell_shop->get_the_other_good(supplied_good) == demanded_good) {
            X = sell_shop->get_price_supply(supplied_good);
        } else {
            if (get_source() != nullptr) {
                Shop* buy_shop = get_source();
                if (sell_shop->get_the_other_good(supplied_good) == buy_shop->get_the_other_good(demanded_good)) {
                    X = sell_shop->get_price_supply(supplied_good) * buy_shop->get_price_demand(demanded_good);
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
        shop.g[1 - q] = demanded_good;
        shop.g[q]     = supplied_good;
        // owner links
        set_outlet_idx(shop.idx);
        set_source_idx(0);
        set_familyshop_idx(shop.idx);
        shop.owner = idx;
        return true;
    }
}

void Trader::clear_shops() {
    set_outlet(nullptr);
    set_source(nullptr);
    set_familyshop(nullptr);
}

void Trader::set_supplies(int s) {
    // Simulation rule: Sets the supply good for the trader and updates the orientation (q) for shop assignment.
    supplied_good = s;
    q = (supplied_good > demanded_good);
}

void Trader::set_demands(int d) {
    // Simulation rule: Sets the demand good for the trader and updates the orientation (q) for shop assignment.
    demanded_good = d;
    q = (supplied_good > demanded_good);
}

bool Trader::is_compatible_with(const Shop& shop) const {
    // Simulation rule: Trader is compatible with a shop if it is active and trades either the supply or demand good.
    return shop.is_active() && (shop.provides(supplied_good) || shop.provides(demanded_good));
}

bool Trader::is_compatible_with(const Shop* shop) const {
    // Simulation rule: Trader is compatible with a shop pointer if it is active and trades either the supply or demand good.
    return shop->is_active() && (shop->provides(supplied_good) || shop->provides(demanded_good));
}

bool Trader::allows_barter_with(const Shop& shop) const {
    // Simulation rule: Trader allows barter with a shop if it is active and trades both the supply and demand goods.
    return shop.is_active() && shop.provides(supplied_good) && shop.provides(demanded_good);
}

bool Trader::wants_to_trade_in(int good) {
    // Simulation rule: Trader wants to trade in if the good matches their demand.
    return demanded_good == good;
}

bool Trader::wants_to_trade_out(int good) {
    // Simulation rule: Trader wants to trade out if the good matches their supply.
    return supplied_good == good;
}

void Trader::set_source(Shop* shop) {
    // Simulation rule: Sets the buyer shop pointer and index for the trader.
    source = shop; source_idx = shop ? shop->idx : 0;
}

void Trader::set_outlet(Shop* shop) {
    // Simulation rule: Sets the seller shop pointer and index for the trader.
    outlet = shop; outlet_idx = shop ? shop->idx : 0;
}

void Trader::set_familyshop(Shop* shop) {
    // Set the pointer to the family shop (owned shop).
    family_shop = shop;
    familyshop = shop ? shop->idx : 0;
}
