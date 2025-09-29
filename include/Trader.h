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
    int seller_idx{0};
    int buyer_idx{0};
    int familyshop{0};      // owned shop index (0 if none)

    // Getters
    int get_supplies() const { return supplies; }
    int get_demands() const { return demands; }

    // Setters
    void set_supplies(int s);
    void set_demands(int d);

    std::string to_string() const;
    void sever_links(Shop& shop);
    int comrade(const std::vector<std::vector<int>>& produces) const;
    int soulmate(const std::vector<std::vector<int>>& consumes) const;
    double utility(std::vector<Shop>& shops) const;
    bool open_shop(Shop& shop);

private:
    int supplies{0};        // produced good
    int demands{0};         // desired good
};

#endif // TRADER_H
