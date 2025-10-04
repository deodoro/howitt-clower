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

/**
 * @brief The Trader class models an agent in the simulation who produces and consumes goods.
 *
 * Simulation rules:
 * - Each trader is assigned a supply good and a demand good.
 * - Traders interact with shops to buy and sell goods, and may open their own shop.
 * - Traders form links to shops for trading, and sever links if shops become unprofitable or inactive.
 * - Utility is calculated based on attainable consumption via current shop links, considering direct barter and indirect trade.
 * - The math in methods encodes these rules, ensuring traders adapt their strategies to maximize utility and survive in the simulated economy.
 */
class Trader {
public:
    int idx{0};             // PROVISIONAL: to make the code compatible while it's refactored
    int q{0};               // 0 if s<=d else 1 (orders shop's good pair)

    /**
     * Get the pointer to the buyer shop.
     */
    Shop* get_source() const { return source; }
    /**
     * Set the pointer to the buyer shop and update the index.
     */
    void set_source(Shop* shop);
    /**
     * Get the pointer to the seller shop.
     */
    Shop* get_outlet() const { return outlet; }
    /**
     * Set the pointer to the seller shop and update the index.
     */
    void set_outlet(Shop* shop);
    /**
     * Get the pointer to the family shop (owned shop).
     */
    Shop* get_familyshop() const { return family_shop; }
    /**
     * Set the pointer to the family shop (owned shop).
     */
    void set_familyshop(Shop* shop);
    /**
     * Set the reference to the shops vector.
     */
    void set_shops(std::vector<Shop>* s) { shops_ref = s; }

    /**
     * Get the supply good for the trader.
     */
    int get_supplied_good() const { return supplied_good; }
    /**
     * Get the demand good for the trader.
     */
    int get_demand_good() const { return demanded_good; }

    /**
     * Set the supply good for the trader and update orientation.
     */
    void set_supplies(int s);
    /**
     * Set the demand good for the trader and update orientation.
     */
    void set_demands(int d);

    /**
     * Check if the trader is compatible with a shop (active and trades supply or demand good).
     */
    bool is_compatible_with(const Shop &shop) const;
    /**
     * Check if the trader is compatible with a shop pointer (active and trades supply or demand good).
     */
    bool is_compatible_with(const Shop *shop) const;

    /**
     * Returns a string representation of the trader's state for debugging.
     * Includes supply/demand goods, shop links, and family shop status.
     */
    std::string to_string() const;
    /**
     * Sever links to a shop if it becomes inactive or unprofitable. Null = sever to all.
     */
    void sever_links(Shop* shop = nullptr);
    /**
     * Select another trader who produces the same good (potential direct competitor or trading partner).
     */
    int any_comrade(const std::vector<std::vector<int>>& produces) const;
    /**
     * Finds a comrade for potential trade, using the same supply good.
     */
    int trade_comrade(const std::vector<std::vector<int>>& produces) const;
    /**
     * Select another trader who consumes the same good (potential indirect trading partner).
     */
    int soulmate(const std::vector<std::vector<int>>& consumes) const;
    /**
     * Calculates the attainable consumption for this trader via current shop links.
     * Considers direct barter and indirect trade.
     */
    double utility() const;
    /**
     * Check if the trader allows barter with a shop (active and trades both supply and demand goods).
     */
    bool allows_barter_with(const Shop& shop) const;
    /**
     * Check if the trader wants to trade in a good (matches demand).
     */
    bool wants_to_trade_in(int good);
    /**
     * Check if the trader wants to trade out a good (matches supply).
     */
    bool wants_to_trade_out(int good);
    /**
     * Open a shop if inactive, set goods and ownership, and link to self.
     */
    bool open_shop(Shop &shop);

private:
    int supplied_good{0};        // produced good
    int demanded_good{0};         // desired good
    int source_idx{0};
    int outlet_idx{0};
    int familyshop_idx{0};      // owned shop index (0 if none)
    std::vector<Shop>* shops_ref{nullptr};
    // TODO: confirm if outlet/source are the in same direction as in the paper
    Shop* source{nullptr};
    Shop* outlet{nullptr};
    Shop* family_shop{nullptr};
};

#endif // TRADER_H
