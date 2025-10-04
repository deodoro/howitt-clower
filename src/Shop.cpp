/*
 * Shop.cpp
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

#include "Shop.h"

// Shop class models a trading post in the simulation.
// Each shop has two goods it can trade (g[0], g[1]), an owner, posted prices (P[0], P[1]),
// weekly incomes (y[0], y[1]), and target trading volumes (tr[0], tr[1]).
// The simulation rules encoded here govern how shops interact with traders, update prices,
// accumulate income, and determine profitability and survival.

std::string Shop::to_string() const {
    // Returns a string representation of the shop's state for debugging.
    return "Shop{active=" + std::to_string(active) + ", g=[" + std::to_string(g[0]) + "," + std::to_string(g[1]) +
           "], owner=" + std::to_string(owner) + ", P=[" + std::to_string(P[0]) + "," + std::to_string(P[1]) +
           "], y=[" + std::to_string(y[0]) + "," + std::to_string(y[1]) + "], tr=[" + std::to_string(tr[0]) + "," + std::to_string(tr[1]) + "]}";
}

bool Shop::provides(int good) const {
    // Checks if the shop trades the specified good.
    return g[0] == good || g[1] == good;
}

void Shop::clear() {
    // Resets the shop to an inactive state, clearing all goods, prices, incomes, and targets.
    active = 0;
    g[0] = g[1] = 0;
    owner = 0;
    P[0] = P[1] = 0.0;
    y[0] = y[1] = 0.0;
    tr[0] = tr[1] = 0.0;
}

int Shop::get_the_other_good(int side) const {
    // Returns the good traded on the opposite side of the shop.
    // Used to determine barter possibilities.
    return g[g[0] == side];
}

double Shop::get_income(int side, bool opposite) const {
    // Returns the weekly income for the specified side.
    // If opposite is true, returns income for the other good.
    if (opposite)
        return y[g[0] != side];
    else
        return y[g[0] == side];
}

void Shop::add_income(int side, double val, bool opposite) {
    // Adds income to the shop for the specified good and side.
    // Used when traders complete transactions.
    int idx;
    if (opposite)
        idx = g[0] != side;
    else
        idx = g[0] == side;
    y[idx] += val;
}

double Shop::get_price_demand(int side) const {
    // Returns the posted price for the good demanded by the trader.
    return P[g[0] == side];
}

double Shop::get_price_supply(int side) const {
    // Returns the posted price for the good supplied by the trader.
    return P[g[0] != side];
}

void Shop::update_prices(double C, std::function<double(int)> overhead_f) {
    // Updates posted prices for both goods based on current targets, costs, and overhead.
    // Simulation rule: Prices are set to maximize profit, but are zero if not profitable.
    auto priceF = [C](double tr0, double tr1, double f_other) {
        return (tr1 - f_other - C > 0.0) ? ((tr1 - C - f_other) / tr0) : 0.0;
    };
    P[0] = priceF(tr[0], tr[1], overhead_f(g[1]));
    P[1] = priceF(tr[1], tr[0], overhead_f(g[0]));
}

void Shop::update_targets(double alpha) {
    // Adjusts target trading volumes toward observed weekly incomes.
    // Simulation rule: Shops adapt their expectations based on recent performance.
    tr[0] += alpha * (y[0] - tr[0]);
    tr[1] += alpha * (y[1] - tr[1]);
}

void Shop::set_targets(double targ0, double targ1) {
    // Directly sets target trading volumes for both goods.
    tr[0] = targ0;
    tr[1] = targ1;
}

bool Shop::is_profitable(std::function<double(int)> overhead) const {
    // Determines if the shop is profitable for both goods, considering prices, incomes, and overhead costs.
    // Simulation rule: Shops exit if not profitable on both sides.
    return (y[0] - P[1] * y[1] - (overhead ? overhead(g[0]) : 0.0)) > 0 &&
           (y[1] - P[0] * y[0] - (overhead ? overhead(g[1]) : 0.0)) > 0;
}

void Shop::reset_weekly_incomes() {
    // Resets weekly incomes to zero at the start of each simulation week.
    if (active) {
        y[0] = y[1] = 0.0;
    }
}
