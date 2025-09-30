/*
 * Trader.h
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
    int q{0};               // 0 if s<=d else 1 (orders shop's good pair)

    // Getters and setters for buyer_idx and seller_idx
    int get_buyer_idx() const { return buyer_idx; }
    void set_buyer_idx(int idx) { 
        buyer_idx = idx; 
        buyer_shop = (idx > 0 && shops_ref) ? &(*shops_ref)[idx] : nullptr; 
    }
    int get_seller_idx() const { return seller_idx; }
    void set_seller_idx(int idx) { 
        seller_idx = idx; 
        seller_shop = (idx > 0 && shops_ref) ? &(*shops_ref)[idx] : nullptr; 
    }

    // Getters and setters for familyshop
    int get_familyshop() const { return familyshop; }
    void set_familyshop(int idx) { 
        familyshop = idx; 
        family_shop = (idx > 0 && shops_ref) ? &(*shops_ref)[idx] : nullptr; 
    }

    // Getters for shop pointers
    Shop* get_buyer_shop() const { return buyer_shop; }
    Shop* get_seller_shop() const { return seller_shop; }
    Shop* get_family_shop() const { return family_shop; }

    // Set shops reference
    void set_shops(std::vector<Shop>* s) { shops_ref = s; }

    // Getters
    int get_supplies() const { return supplies; }
    int get_demands() const { return demands; }

    // Setters
    void set_supplies(int s);
    void set_demands(int d);

    bool is_compatible_with(const Shop &shop) const;
    bool is_compatible_with(const Shop *shop) const;

    std::string to_string() const;
    void sever_links(Shop& shop);
    int comrade(const std::vector<std::vector<int>>& produces) const;
    int soulmate(const std::vector<std::vector<int>>& consumes) const;
    double utility(std::vector<Shop>& shops) const;
    bool allows_barter_with(const Shop& shop) const;
    bool wants_to_trade_in(int good);
    bool wants_to_trade_out(int good);
    bool open_shop(Shop &shop);

private:
    int supplies{0};        // produced good
    int demands{0};         // desired good
    int buyer_idx{0};
    int seller_idx{0};
    int familyshop{0};      // owned shop index (0 if none)
    std::vector<Shop>* shops_ref{nullptr};
    Shop* buyer_shop{nullptr};
    Shop* seller_shop{nullptr};
    Shop* family_shop{nullptr};
};

#endif // TRADER_H
