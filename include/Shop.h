/*
 * Shop.h
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

#ifndef SHOP_H
#define SHOP_H

#include <string>
#include <functional>

class Shop {
public:
    int idx{0};
    int active{0};         // 1 if active
    int g[2]{0,0};         // traded goods, ordered so g[0]<=g[1]
    int owner{0};          // trader id that owns this shop
    double P[2]{0.0,0.0};  // posted buying prices per side
    double y[2]{0.0,0.0};  // realized incomes this week per side
    double tr[2]{0.0,0.0}; // income targets per side

    std::string to_string() const;
    bool provides(int good) const;
    void clear();
    int get_good(int side, bool opposite = false) const;
    double get_income(int side, bool opposite = false) const;
    void add_income(int side, double val, bool opposite = false);
    double get_price(int side, bool opposite = false) const;
    void update_prices(double C, std::function<double(int)> overhead_f);
    void update_targets(double alpha);
    void update_targets(double targ0, double targ1);
    bool is_profitable(std::function<double(int)> overhead = nullptr) const;
    void reset_weekly_incomes();
};

#endif // SHOP_H
