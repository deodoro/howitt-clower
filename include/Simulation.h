/*
 * Simulation.h
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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "PCG32.h"
#include "Shop.h"
#include "Trader.h"
#include "ResearchResults.h"
#include "RunInfo.h"
#include "SimulationInfo.h"
#include <vector>
#include <array>
#include <ctime>
#include <cstdio>
#include <functional>

#include "SimulationSerializer.h"

struct MatchEvaluation;

/**
 * @brief The Simulation class orchestrates the agent-based trading economy.
 *
 * Simulation rules:
 * - Initializes a population of traders and shops, assigning goods and locations.
 * - Runs the simulation over multiple time steps (weeks), allowing agents to interact, trade, and adapt.
 * - Implements entry, matching, trading, price updates, and exit processes for shops and traders.
 * - Tracks aggregate statistics, monetary emergence, and equilibrium conditions.
 * - The math in methods encodes these rules, ensuring the simulated market evolves according to agent decisions and economic constraints.
 */
class Simulation {
    friend class SimulationSerializer;
public:
    /**
     * Constructor builds the static economy (agents and lists).
     */
    Simulation(SimulationInfo info = SimulationInfo());

    /**
     * Runs all simulation sweeps and time steps, collecting results and statistics.
     */
    void run_all();

    /**
     * Get the collected run info for all slopes and runs.
     */
    const std::vector<std::vector<RunInfo*>>& get_runs_per_slope() const { return runs_per_slope; }

private:
    SimulationInfo info;

    double slope{16.0};

    std::vector<std::vector<RunInfo*>> runs_per_slope; // stores all runs for each slope

    // State
    std::vector<Trader> traders;
    std::vector<Shop> shops;
    std::vector<Trader*> trader_line; // Weekly lineup - permutation of traders

    std::vector<std::vector<int>> produces; // by good i: trader ids producing i
    std::vector<std::vector<int>> consumes; // by good j: trader ids desiring j

    // Helper methods
    /**
     * Computes the overhead cost function for shops based on good index and slope parameter.
     * Simulation rule: Overhead increases with good index and slope, affecting shop profitability.
     */
    std::function<double(int)> make_overhead(double f1, double slope);
    /**
     * Initializes the static structure of the economy: traders, shops, and goods assignments.
     * Simulation rule: Sets up agent populations and their supply/demand relationships.
     */
    void init_static();
    /**
     * Initializes the state for a new simulation run, resetting agents and shops.
     * Simulation rule: Prepares the market for a fresh round of trading and adaptation.
     */
    void init_run(RunInfo& runInfo);
    /**
     * Randomizes the weekly lineup of traders for matching and trading.
     * Simulation rule: Ensures agent interactions are permuted each week.
     */
    void lineup();
    /**
     * Handles entry of new shops by prospective entrepreneurs.
     * Simulation rule: Traders may open shops if conditions are favorable.
     */
    void weekly_entry(RunInfo& runInfo);
    /**
     * Matches traders to shops based on utility and compatibility, updating links.
     * Simulation rule: Agents seek optimal trading partners each week.
     */
    std::vector <MatchEvaluation>* weekly_matching();
    /**
     * Reports the state of a trader for debugging and analysis.
     */
    void report_trader(Trader const* trader);
    /**
     * Executes weekly trading, income accounting, and stochastic exit of unprofitable shops.
     * Simulation rule: Shops may exit if not profitable, severing links with traders.
     */
    void weekly_trade_and_exit(RunInfo& runInfo);
    /**
     * Updates shop targets and posted prices adaptively based on recent performance.
     * Simulation rule: Shops adjust strategies to maximize profit and survive.
     */
    void weekly_update_prices();
    /**
     * Simulates the research process for a prospective shop owner, determining entry viability.
     * Simulation rule: Entry depends on utility comparisons with existing trading options.
     */
    ResearchResults research(Trader& trader);
    /**
     * Selects a random consumer of a given good from the trader population.
     */
    Trader &random_consumer(int good);
    /**
     * Selects a random producer of a given good from the trader population.
     */
    Trader &random_producer(int good);
    /**
     * Selects a random shop from the shop population.
     */
    Shop &random_shop();
    /**
     * Attempts to match a trader with shops for barter opportunities.
     * Simulation rule: Barter is considered if compatible shops are available.
     */
    void try_barter(const Trader* trader, std::vector<int>& c, struct MatchEvaluation& eval);
    /**
     * Attempts to improve a trader's outlet or source by matching with candidate shops.
     * Simulation rule: Agents seek better trading terms through shop selection.
     */
    void try_one(const Trader* trader, std::vector<int>& c, struct MatchEvaluation& eval);
    /**
     * Attempts to match a trader with two shops for indirect trade via a common intermediary.
     * Simulation rule: Indirect trade is considered for maximizing utility.
     */
    void try_two(const Trader* trader, std::vector<int>& c, struct MatchEvaluation& eval);
    /**
     * Calculates monetary equilibrium and tracks the emergence of money in the simulation.
     * Simulation rule: Aggregates statistics to detect monetary phases and equilibrium.
     */
    int calc1(RunInfo& runInfo);
    /**
     * Calculates final statistics for the simulation run, including surplus and shop counts.
     * Simulation rule: Summarizes market outcomes for analysis.
     */
    void calc2(RunInfo& runInfo);
    /**
     * Computes root mean square error for price and volume statistics.
     * Simulation rule: Measures market efficiency and price dispersion.
     */
    void rmse(RunInfo& runInfo, std::vector<double> &Pinv);
    /**
     * Reports simulation progress and statistics at specified intervals.
     */
    void report(RunInfo& runInfo);
    /**
     * Prints debug information about the simulation state.
     */
    void print_debug(std::string title) const;
    /**
     * Prints the state of all traders for debugging.
     */
    void print_traders() const;
    /**
     * Prints the state of all shops for debugging.
     */
    void print_shops() const;
};

#endif // SIMULATION_H
