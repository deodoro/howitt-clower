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
    void update_targets(double targ0, double targ1);
    bool is_profitable(std::function<double(int)> overhead = nullptr) const;
    void reset_weekly_incomes();
};

#endif // SHOP_H
