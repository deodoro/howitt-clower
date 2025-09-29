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

std::string Shop::to_string() const {
    return "Shop{active=" + std::to_string(active) + ", g=[" + std::to_string(g[0]) + "," + std::to_string(g[1]) +
           "], owner=" + std::to_string(owner) + ", P=[" + std::to_string(P[0]) + "," + std::to_string(P[1]) +
           "], y=[" + std::to_string(y[0]) + "," + std::to_string(y[1]) + "], tr=[" + std::to_string(tr[0]) + "," + std::to_string(tr[1]) + "]}";
}

bool Shop::provides(int good) const {
    return g[0] == good || g[1] == good;
}

void Shop::clear() {
    active = 0;
    g[0] = g[1] = 0;
    owner = 0;
    P[0] = P[1] = 0.0;
    y[0] = y[1] = 0.0;
    tr[0] = tr[1] = 0.0;
}

int Shop::get_good(int side, bool opposite) const {
    if (opposite)
        return g[g[0] != side];
    else
        return g[g[0] == side];
}

double Shop::get_income(int side, bool opposite) const {
    if (opposite)
        return y[g[0] != side];
    else
        return y[g[0] == side];
}

void Shop::add_income(int side, double val, bool opposite) {
    int idx;
    if (opposite)
        idx = g[0] != side;
    else
        idx = g[0] == side;
    y[idx] += val;
}

double Shop::get_price(int side, bool opposite) const {
    if (opposite)
        return P[g[0] != side];
    else
        return P[g[0] == side];
}

void Shop::update_prices(double C, std::function<double(int)> overhead_f) {
    auto priceF = [C](double tr0, double tr1, double f_other) {
        return (tr1 - f_other - C > 0.0) ? ((tr1 - C - f_other) / tr0) : 0.0;
    };
    P[0] = priceF(tr[0], tr[1], overhead_f(g[1]));
    P[1] = priceF(tr[1], tr[0], overhead_f(g[0]));
}

void Shop::update_targets(double alpha) {
    tr[0] += alpha * (y[0] - tr[0]);
    tr[1] += alpha * (y[1] - tr[1]);
}

void Shop::set_targets(double targ0, double targ1) {
    tr[0] = targ0;
    tr[1] = targ1;
}

bool Shop::is_profitable(std::function<double(int)> overhead) const {
    return (y[0] - P[1] * y[1] - (overhead ? overhead(g[0]) : 0.0)) > 0 && 
           (y[1] - P[0] * y[0] - (overhead ? overhead(g[1]) : 0.0)) > 0;
}

void Shop::reset_weekly_incomes() {
    if (active) {
        y[0] = y[1] = 0.0;
    }
}
