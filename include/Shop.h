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

class Trader;

/**
 * @brief The Shop class models a trading post in the simulation.
 *
 * Simulation rules:
 * - Each shop trades two goods (g[0], g[1]), and is owned by a trader.
 * - Shops post prices for each good, track weekly incomes, and set target trading volumes.
 * - Shops interact with traders, accumulate income from transactions, and adapt their prices and targets based on performance.
 * - Shops exit the market if not profitable, based on realized incomes, posted prices, and overhead costs.
 * - The math in methods encodes these rules, ensuring shops only survive if they contribute positively to the simulated economy.
 */
class Shop {
private:
    Trader* owner{nullptr};          // trader that owns this shop
public:
    int idx{0};
    int active{0};         // 1 if active
    int g[2]{0,0};         // traded goods, ordered so g[0]<=g[1]
    double P[2]{0.0,0.0};  // posted buying prices per side
    double y[2]{0.0,0.0};  // realized incomes this week per side
    double tr[2]{0.0,0.0}; // income targets per side

    Shop() {};
    /**
     * Returns a string representation of the shop's state for debugging.
     * Includes active status, traded goods, owner, prices, incomes, and targets.
     */
    std::string to_string() const;
    /**
     * Checks if the shop trades the specified good.
     */
    bool provides(int good) const;
    /**
     * Resets the shop to an inactive state, clearing all goods, prices, incomes, and targets.
     */
    void clear();
    /**
     * Returns the good traded on the opposite side of the shop.
     * Used to determine barter possibilities.
     */
    int get_the_other_good(int side) const;
    /**
     * Returns the weekly income for the specified side.
     * If opposite is true, returns income for the other good.
     */
    double get_income(int side, bool opposite = false) const;
    /**
     * Adds income to the shop for the specified good and side.
     * Used when traders complete transactions.
     */
    void add_income(int side, double val, bool opposite = false);
    /**
     * Returns the posted price for the good demanded by the trader.
     */
    double get_price_demand(int side) const;
    /**
     * Returns the posted price for the good supplied by the trader.
     */
    double get_price_supply(int side) const;
    /**
     * Updates posted prices for both goods based on current targets, costs, and overhead.
     * Simulation rule: Prices are set to maximize profit, but are zero if not profitable.
     */
    void update_prices(double C, std::function<double(int)> overhead_f);
    /**
     * Adjusts target trading volumes toward observed weekly incomes.
     * Simulation rule: Shops adapt their expectations based on recent performance.
     */
    void update_targets(double alpha);
    /**
     * Directly sets target trading volumes for both goods.
     */
    void set_targets(double targ0, double targ1);
    /**
     * Determines if the shop is profitable for both goods, considering prices, incomes, and overhead costs.
     * Simulation rule: Shops exit if not profitable on both sides.
     */
    bool is_profitable(std::function<double(int)> overhead = nullptr) const;
    /**
     * Resets weekly incomes to zero at the start of each simulation week.
     */
    void reset_weekly_incomes();
    /**
     * Returns true if the shop is active.
     */
    bool is_active() const { return active == 1; }
    /**
     * Sets the owner of the shop.
     */
    void set_owner(Trader*);
    /**
     * Returns the owner of the shop.
     */
    Trader* get_owner() const { return owner; }
};

#endif // SHOP_H
