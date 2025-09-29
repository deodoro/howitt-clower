import numpy as np
import random
from collections import defaultdict, Counter

# ===============================
# Parameters (close to Table 1)
# ===============================
T = 5000                 # number of weeks (lower than paper default for speed)
n = 10                   # number of commodities
b = 24                   # per (i,j)-type -> m = 2160 for n=10
m = b * n * (n - 1)      # number of agents
K = 200                  # shop locations
xMax = 200               # max initial target "animal spirits"
search_intensity = 0.05  # weekly probability an agent searches
adaptation_speed = 0.25  # targets adjust toward realized incomes
exit_rate = 0.01         # weekly exit prob if operating surplus < 0
C_setup = 5              # one-time setup cost (charged on entry only)
s_overhead = 8           # slope of overhead cost by commodity index

# Reproducibility
RANDOM_SEED = 0
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)

# ===============================
# Overhead schedule
# ===============================
def overhead(g_idx: int) -> float:
    # Increasing fixed weekly overhead per traded commodity
    return s_overhead * (g_idx + 1)

# ===============================
# Data structures
# ===============================
class Transactor:
    __slots__ = ("id","prod","cons","owned_shop_id","outlet","source")
    def __init__(self, idx, prod, cons):
        self.id = idx
        self.prod = prod
        self.cons = cons
        self.owned_shop_id = None
        self.outlet = None       # shop where they SELL prod
        self.source = None       # shop where they BUY cons

class Shop:
    """
    Trades goods (g0, g1). We keep two posted prices:
      p01 = price of g0 in units of g1  (g1 per 1 unit g0)
      p10 = price of g1 in units of g0  (g0 per 1 unit g1)
    Targets (tr0, tr1) are income targets in commodities g0 and g1, respectively.
    We maintain last week's realized volumes (q0_hat, q1_hat) to set prices.
    """
    __slots__ = ("id","g0","g1","owner","tr0","tr1","p01","p10",
                 "q0_hat","q1_hat","y0","y1","customers_sell_g0","customers_buy_g0",
                 "customers_sell_g1","customers_buy_g1","setup_paid")
    def __init__(self, idx, g0, g1, owner, tr0, tr1):
        self.id = idx
        self.g0, self.g1 = g0, g1
        self.owner = owner
        self.tr0, self.tr1 = float(tr0), float(tr1)
        # initialize prices with something positive
        self.p01 = 1.0
        self.p10 = 1.0
        # last week's volumes (expectations for pricing)
        self.q0_hat = 1.0
        self.q1_hat = 1.0
        # last week's realized incomes
        self.y0 = 0.0
        self.y1 = 0.0
        # client sets
        self.customers_sell_g0 = set()
        self.customers_buy_g0 = set()
        self.customers_sell_g1 = set()
        self.customers_buy_g1 = set()
        # charge setup once
        self.setup_paid = False

    def set_prices_full_cost(self):
        """
        'Full-cost' toward targets: set prices high enough
        (given expected volume) to cover targets and overhead on the *other* side.
        Intuition: income in g0 is primarily from selling g1 (at price p10); to
        deliver target tr0 and cover overhead on g0-side flows, raise p10 when
        expected q1_hat is low, and analogously for p01.
        """
        f0 = overhead(self.g0)
        f1 = overhead(self.g1)
        q0h = max(self.q0_hat, 1.0)
        q1h = max(self.q1_hat, 1.0)
        # keep prices nonnegative and finite
        self.p01 = max((self.tr1 + f1) / q0h, 0.0)  # g1 per g0 unit
        self.p10 = max((self.tr0 + f0) / q1h, 0.0)  # g0 per g1 unit

    def update_targets(self, y0_real, y1_real):
        self.y0, self.y1 = y0_real, y1_real
        self.tr0 = (1 - adaptation_speed) * self.tr0 + adaptation_speed * y0_real
        self.tr1 = (1 - adaptation_speed) * self.tr1 + adaptation_speed * y1_real

    def weekly_reset(self):
        self.customers_sell_g0.clear()
        self.customers_buy_g0.clear()
        self.customers_sell_g1.clear()
        self.customers_buy_g1.clear()

# ===============================
# Population & initial state
# ===============================
def create_transactors():
    T_list = []
    idx = 0
    for i in range(n):
        for j in range(n):
            if i == j: 
                continue
            for _ in range(b):
                T_list.append(Transactor(idx, i, j))
                idx += 1
    return T_list

# ===============================
# Helper: expected throughput proxy
# ===============================
def expected_q_for_pair(shops, pair):
    """Return (q0_hat, q1_hat) for the shop trading this pair if it exists,
       else small default expectations to let pricing remain finite."""
    sid = pair_to_shop_id.get(pair)
    if sid is None:
        return (1.0, 1.0)
    sh = shops[sid]
    return (max(sh.q0_hat, 1.0), max(sh.q1_hat, 1.0))

# ===============================
# Agent search: best path (direct or via m)
# ===============================
def best_path_for_agent(agent, shops):
    """
    Return (best_outlet, best_source, expected_consumption, intermediary_or_None).
    Consumption measured in units of 'cons' good per week (proxy).
    We compare:
      - Direct: one shop (prod, cons) if present (outlet and source may be same or distinct)
      - Indirect via m in {0..n-1}, m != prod != cons:
          outlet uses (prod, m), source uses (m, cons)
    Expected flow proxy is the min of the relevant shop q-hats on the constrained sides.
    """
    prod, cons = agent.prod, agent.cons
    best_tuple = (agent.outlet, agent.source, 0.0, None)

    # DIRECT
    sid_direct = pair_to_shop_id.get((prod, cons))
    if sid_direct is not None:
        # proxy: available flow limited by both sides at that shop
        sh = shops[sid_direct]
        # for direct, seller side is g0=prod to g1=cons, buyer side for cons is the same shop
        q_proxy = min(max(sh.q0_hat, 0.0), max(sh.q1_hat, 0.0))
        if q_proxy > best_tuple[2]:
            best_tuple = (sid_direct, sid_direct, q_proxy, None)

    # INDIRECT via m
    for m_good in range(n):
        if m_good == prod or m_good == cons:
            continue
        sid_out = pair_to_shop_id.get((prod, m_good))
        sid_src = pair_to_shop_id.get((m_good, cons))
        if sid_out is None or sid_src is None:
            continue
        sh_out = shops[sid_out]
        sh_src = shops[sid_src]
        # flow limited by the weaker link's q-hat
        q_proxy = min(max(sh_out.q0_hat, 0.0), max(sh_src.q0_hat, 0.0),
                      max(sh_out.q1_hat, 0.0), max(sh_src.q1_hat, 0.0))
        if q_proxy > best_tuple[2]:
            best_tuple = (sid_out, sid_src, q_proxy, m_good)

    return best_tuple

# ===============================
# Entry (market research) & Exit
# ===============================
def consider_entry(transactors, shops, available_shop_ids):
    """Entrepreneur proposes shop (prod, cons) using owner's types; 
       opens only if expected surplus >= 0 under tentative prices."""
    if len(shops) >= K or not available_shop_ids:
        return
    candidates = [r for r in transactors if r.owned_shop_id is None]
    if not candidates:
        return
    r = random.choice(candidates)
    g0, g1 = r.prod, r.cons
    if (g0, g1) in pair_to_shop_id:
        return  # shop already exists

    new_id = available_shop_ids.pop(0)
    tr0 = random.randint(1, xMax)
    tr1 = random.randint(1, xMax)
    sh = Shop(new_id, g0, g1, r.id, tr0, tr1)

    # crude "market research": look at currently unmatched types on both sides
    # and set initial q-hats as a function of potential clientele size
    sellers_g0 = sum(1 for a in transactors if a.prod == g0 and a.outlet in (None, -1))
    buyers_g0  = sum(1 for a in transactors if a.cons == g0 and a.source in (None, -1))
    sellers_g1 = sum(1 for a in transactors if a.prod == g1 and a.outlet in (None, -1))
    buyers_g1  = sum(1 for a in transactors if a.cons == g1 and a.source in (None, -1))
    sh.q0_hat = max(1.0, min(sellers_g0, buyers_g0))
    sh.q1_hat = max(1.0, min(sellers_g1, buyers_g1))
    sh.set_prices_full_cost()

    # Expected incomes under these q-hats
    # Interpret incomes as receipts in opposite commodity from the paired market:
    #   y0 ≈ p10 * q1_hat   (receipts in g0 from selling g1)
    #   y1 ≈ p01 * q0_hat   (receipts in g1 from selling g0)
    y0_exp = sh.p10 * sh.q1_hat
    y1_exp = sh.p01 * sh.q0_hat
    f0, f1 = overhead(g0), overhead(g1)
    surplus_exp = y0_exp + y1_exp - (f0 + f1)

    if surplus_exp >= 0.0:
        shops[new_id] = sh
        pair_to_shop_id[(g0, g1)] = new_id
        r.owned_shop_id = new_id
        sh.setup_paid = True  # one-time
    else:
        available_shop_ids.insert(0, new_id)

def consider_exit(shops, transactors, available_shop_ids):
    to_close = []
    for sid, sh in shops.items():
        f0, f1 = overhead(sh.g0), overhead(sh.g1)
        surplus = (sh.y0 + sh.y1) - (f0 + f1)
        if surplus < 0.0 and random.random() < exit_rate:
            to_close.append(sid)
    for sid in to_close:
        owner = shops[sid].owner
        transactors[owner].owned_shop_id = None
        pair_to_shop_id.pop((shops[sid].g0, shops[sid].g1), None)
        del shops[sid]
        available_shop_ids.append(sid)

# ===============================
# Weekly clearing with rationing
# ===============================
def weekly_clear(shops):
    """
    Each shop clears *two* markets in parallel:
      - g0-for-g1 (sellers of g0 vs buyers of g0)
      - g1-for-g0 (sellers of g1 vs buyers of g1)
    Each agent supplies/demands 1 unit; realized quantities are rationed
    proportionally by the short side.
    We update:
      - realized volumes q0_real, q1_real (throughput)
      - realized incomes: y0 = p10 * q1_real ; y1 = p01 * q0_real
      - q-hats for next period (exponential smoothing by simply assigning)
    """
    total_y0 = total_y1 = 0.0
    for sh in shops.values():
        # collect side counts
        S0 = len([a for a in sh.customers_sell_g0])
        B0 = len([a for a in sh.customers_buy_g0])
        S1 = len([a for a in sh.customers_sell_g1])
        B1 = len([a for a in sh.customers_buy_g1])

        # volumes per side = min(sellers, buyers)
        q0_real = float(min(S0, B0))  # g0 traded
        q1_real = float(min(S1, B1))  # g1 traded

        # realized incomes in g0 and g1
        y0 = sh.p10 * q1_real  # receipts in g0 from selling g1
        y1 = sh.p01 * q0_real  # receipts in g1 from selling g0
        sh.update_targets(y0, y1)

        # expectations for pricing next week
        sh.q0_hat = max(q0_real, 1.0)
        sh.q1_hat = max(q1_real, 1.0)

        total_y0 += y0
        total_y1 += y1
    return total_y0, total_y1

# ===============================
# Money detection (intermediary share)
# ===============================
def update_money_stats(intermediary_used_counts, chosen_intermediary):
    if chosen_intermediary is not None:
        intermediary_used_counts[chosen_intermediary] += 1

def compute_money_share(intermediary_used_counts):
    if not intermediary_used_counts:
        return None, 0.0
    total = sum(intermediary_used_counts.values())
    if total == 0:
        return None, 0.0
    g, cnt = max(intermediary_used_counts.items(), key=lambda kv: kv[1])
    return g, cnt / total

# ===============================
# Simulation
# ===============================
def simulate():
    global pair_to_shop_id
    transactors = create_transactors()
    shops = {}
    pair_to_shop_id = {}  # map (g0,g1) -> shop_id
    available_shop_ids = list(range(K))

    # seed a few shops to avoid cold start
    for _ in range(min(5, K)):
        r = random.choice(transactors)
        if (r.prod, r.cons) in pair_to_shop_id:
            continue
        sid = available_shop_ids.pop(0)
        tr0 = random.randint(1, xMax)
        tr1 = random.randint(1, xMax)
        sh = Shop(sid, r.prod, r.cons, r.id, tr0, tr1)
        sh.q0_hat = 1.0
        sh.q1_hat = 1.0
        sh.set_prices_full_cost()
        sh.setup_paid = True
        shops[sid] = sh
        pair_to_shop_id[(r.prod, r.cons)] = sid
        r.owned_shop_id = sid

    # metrics
    intermediary_used_counts = Counter()
    gdp_series = []
    shop_counts_series = []
    price_dispersion_series = []

    for week in range(T):
        # ------- ENTRY -------
        consider_entry(transactors, shops, available_shop_ids)

        # ------- SHOPPING / SEARCH -------
        # Reset client sets
        for sh in shops.values():
            sh.weekly_reset()

        # random order
        order = np.random.permutation(m)
        for idx in order:
            a = transactors[idx]
            if a.owned_shop_id is not None:
                # owner doesn't shop as a customer
                continue
            # search with probability
            if random.random() < search_intensity or (a.outlet is None or a.source is None):
                best_outlet, best_source, qstar, intermediary = best_path_for_agent(a, shops)
                # switch only if improves expected consumption proxy
                # current expected
                curr_q = 0.0
                if a.outlet is not None and a.source is not None:
                    # approximate current expected cons by min of both shop's q-hats
                    sh_o = shops.get(a.outlet)
                    sh_s = shops.get(a.source)
                    if sh_o and sh_s:
                        curr_q = min(sh_o.q0_hat, sh_o.q1_hat, sh_s.q0_hat, sh_s.q1_hat)
                if best_outlet is not None and best_source is not None and qstar > curr_q:
                    a.outlet = best_outlet
                    a.source = best_source
                    update_money_stats(intermediary_used_counts, intermediary)

            # After decisions, enroll the agent with the chosen shops
            if a.outlet in shops:
                sh = shops[a.outlet]
                # they SELL their production good at outlet
                if sh.g0 == a.prod and sh.g1 != a.prod:
                    sh.customers_sell_g0.add(a.id)
                elif sh.g1 == a.prod and sh.g0 != a.prod:
                    sh.customers_sell_g1.add(a.id)
            if a.source in shops:
                sh = shops[a.source]
                # they BUY their consumption good at source
                if sh.g0 == a.cons:
                    sh.customers_buy_g0.add(a.id)
                elif sh.g1 == a.cons:
                    sh.customers_buy_g1.add(a.id)

        # ------- CLEARING -------
        total_y0, total_y1 = weekly_clear(shops)

        # ------- EXIT -------
        consider_exit(shops, transactors, available_shop_ids)

        # ------- PRICING UPDATE -------
        for sh in shops.values():
            sh.set_prices_full_cost()

        # ------- METRICS -------
        # GDP proxy: sum of (operating surplus)+ across shops
        gdp = 0.0
        prices = []
        for sh in shops.values():
            f0, f1 = overhead(sh.g0), overhead(sh.g1)
            surplus = (sh.y0 + sh.y1) - (f0 + f1)
            gdp += max(surplus, 0.0)
            prices.append(sh.p01)
            prices.append(sh.p10)
        shop_counts_series.append(len(shops))
        gdp_series.append(gdp)
        if prices:
            price_dispersion_series.append(np.std(prices))
        else:
            price_dispersion_series.append(0.0)

        # progress
        if week % 500 == 0:
            money_good, money_share = compute_money_share(intermediary_used_counts)
            print(
                f"Week {week:5d} | Shops: {len(shops):3d} | "
                f"GDP: {gdp:8.2f} | Price σ: {price_dispersion_series[-1]:.3f} | "
                f"Money: {money_good if money_good is not None else '-'} "
                f"({money_share:.2%})"
            )

    # ------- FINAL REPORT -------
    money_good, money_share = compute_money_share(intermediary_used_counts)
    print("\n=== Simulation complete ===")
    print(f"Weeks: {T}")
    print(f"Final #shops: {len(shops)}")
    if money_good is not None:
        print(f"Emergent money candidate: good {money_good} with share {money_share:.2%}")
    else:
        print("No clear monetary intermediary emerged.")
    print(f"Final GDP (weekly proxy): {gdp_series[-1]:.2f}")
    print(f"Avg GDP (last 10%): {np.mean(gdp_series[int(0.9*T):]) if gdp_series else 0.0:.2f}")
    print(f"Avg shops (last 10%): {np.mean(shop_counts_series[int(0.9*T):]) if shop_counts_series else 0.0:.1f}")
    print(f"Price dispersion σ (last): {price_dispersion_series[-1]:.4f}")

if __name__ == "__main__":
    simulate()
