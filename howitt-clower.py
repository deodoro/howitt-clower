#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import time
from collections import Counter, defaultdict

# =========================================
# Defaults (Table 1 in Howitt–Clower 2000)
# =========================================
DEF_T = 1000              # weeks per run
DEF_n = 10                 # number of commodities
DEF_b = 24                 # agents per (i,j)-type; m=b*n*(n-1)=2160 for n=10
DEF_K = 200                # potential shop sites
DEF_xMax = 200             # max initial "animal spirits" (income targets)
DEF_search = 0.05          # search intensity
DEF_lambda = 0.25          # adaptation speed toward realized incomes
DEF_exit = 0.01            # exit probability if operating surplus < 0
DEF_Csetup = 5             # one-time setup cost (bookkeeping only)
DEF_s_overhead = 8         # overhead slope; varies in experiments
SEED_BASE = 1729           # base RNG seed for reproducibility across runs
RAND = 1 

class PCG32:
    def __init__(self, initstate: int, initseq: int):
        self.state = 0
        self.inc   = ((initseq & ((1<<64)-1)) << 1) | 1
        self._next()
        self.state = (self.state + (initstate & ((1<<64)-1))) & ((1<<64)-1)
        self._next()

    def _next(self) -> int:
        oldstate = self.state
        self.state = (oldstate * 6364136223846793005 + self.inc) & ((1<<64)-1)
        xorshifted = (((oldstate >> 18) ^ oldstate) >> 27) & 0xFFFFFFFF
        rot = (oldstate >> 59) & 31
        return ((xorshifted >> rot) | ((xorshifted << ((-rot) & 31)) & 0xFFFFFFFF)) & 0xFFFFFFFF

    def next_u32(self) -> int:
        return self._next()

random_generator = None

def RANDOM_N(n):
    return random_generator.next_u32() % n

def rnd01():
    return (1 + RANDOM_N(1000)) / 1000.0

# ------------- CLI -------------
def parse_args():
    ap = argparse.ArgumentParser(description="Howitt–Clower emergence simulation (extended).")
    ap.add_argument("--T", type=int, default=DEF_T, help="weeks per run")
    ap.add_argument("--n", type=int, default=DEF_n, help="number of goods")
    ap.add_argument("--b", type=int, default=DEF_b, help="agents per (i,j)")
    ap.add_argument("--K", type=int, default=DEF_K, help="shop sites")
    ap.add_argument("--xMax", type=int, default=DEF_xMax, help="max initial target")
    ap.add_argument("--search", type=float, default=DEF_search, help="search intensity")
    ap.add_argument("--lam", type=float, default=DEF_lambda, help="target adaptation speed")
    ap.add_argument("--exit", type=float, default=DEF_exit, help="exit prob when surplus<0")
    ap.add_argument("--s", type=float, default=DEF_s_overhead, help="overhead slope")
    ap.add_argument("--runs", type=int, default=1, help="number of independent runs")
    ap.add_argument("--quiet", action="store_true", help="suppress per-week logs")
    ap.add_argument("--warm_shops", type=int, default=5, help="seed shops at t=0")
    ap.add_argument("--entry_sample", type=int, default=150, help="market-research sample size per side")
    ap.add_argument("--log_every", type=int, default=500, help="progress print interval")
    return ap.parse_args()

# =========================================
# Core model components
# =========================================
class Params:
    def __init__(self, T, n, b, K, xMax, search, lam, exit_rate, Csetup, s_overhead,
                 warm_shops, entry_sample, log_every, quiet):
        self.T=T; self.n=n; self.b=b; self.K=K; self.xMax=xMax
        self.search=search; self.lam=lam; self.exit_rate=exit_rate
        self.Csetup=Csetup; self.s_overhead=s_overhead
        self.warm_shops=warm_shops
        self.entry_sample=entry_sample
        self.log_every=log_every
        self.quiet=quiet
        self.m = b*n*(n-1)

def overhead(g_idx, s_overhead):
    return s_overhead * g_idx  # C version: f(i) = f1 + (i-1)*slope, with f1=0 and i 1-indexed

class Transactor:
    __slots__ = ("id","prod","cons","owned_shop_id","outlet","source")
    def __init__(self, idx, prod, cons):
        self.id = idx
        self.prod = prod
        self.cons = cons
        self.owned_shop_id = None
        self.outlet = None
        self.source = None

class Shop:
    """
    Trades (g0, g1) with posted prices:
      p01 = price of g0 in units of g1 (g1 per 1 unit g0)
      p10 = price of g1 in units of g0 (g0 per 1 unit g1)
    Targets tr0 (income in g0), tr1 (income in g1) adapt to realized incomes.
    """
    __slots__ = ("id","g0","g1","owner","tr0","tr1","p01","p10",
                 "q0_hat","q1_hat","y0","y1",
                 "customers_sell_g0","customers_buy_g0",
                 "customers_sell_g1","customers_buy_g1",
                 "setup_paid","lam","s_overhead")

    def __init__(self, idx, g0, g1, owner, tr0, tr1, lam, s_overhead):
        self.id=idx; self.g0=g0; self.g1=g1; self.owner=owner
        self.tr0=float(tr0); self.tr1=float(tr1)
        self.p01=1.0; self.p10=1.0
        self.q0_hat=1.0; self.q1_hat=1.0
        self.y0=0.0; self.y1=0.0
        self.customers_sell_g0=set(); self.customers_buy_g0=set()
        self.customers_sell_g1=set(); self.customers_buy_g1=set()
        self.setup_paid=False
        self.lam=lam
        self.s_overhead=s_overhead

    def set_prices_full_cost(self):
        f0 = overhead(self.g0, self.s_overhead)
        f1 = overhead(self.g1, self.s_overhead)
        q0h = max(self.q0_hat, 1.0)
        q1h = max(self.q1_hat, 1.0)
        # aim to hit incomes (tr) and cover overheads on opposite side
        self.p01 = max((self.tr1 + f1) / q0h, 0.0)  # g1 / g0
        self.p10 = max((self.tr0 + f0) / q1h, 0.0)  # g0 / g1

    def update_targets(self, y0_real, y1_real):
        lam = self.lam
        self.y0, self.y1 = y0_real, y1_real
        self.tr0 = (1 - lam) * self.tr0 + lam * y0_real
        self.tr1 = (1 - lam) * self.tr1 + lam * y1_real

    def weekly_reset(self):
        self.customers_sell_g0.clear(); self.customers_buy_g0.clear()
        self.customers_sell_g1.clear(); self.customers_buy_g1.clear()

# =========================================
# Simulation utilities
# =========================================
def create_transactors(n, b):
    arr=[]
    idx=0
    for i in range(n):
        for j in range(n):
            if i==j: continue
            for _ in range(b):
                arr.append(Transactor(idx, i, j))
                idx+=1
    return arr

def seed_randomness(seed_val):
    global random_generator
    random_generator = PCG32(seed_val, seed_val)
    return seed_val

def expected_q_for_pair(shops, pair, default=1.0):
    sid = pair_to_shop_id.get(pair)
    if sid is None: return (default, default)
    sh = shops[sid]
    return (max(sh.q0_hat, default), max(sh.q1_hat, default))

def best_path_for_agent(agent, shops, n):
    """Return (outlet_id, source_id, q_proxy, intermediary)."""
    prod, cons = agent.prod, agent.cons
    best = (agent.outlet, agent.source, 0.0, None)

    # direct
    sid = pair_to_shop_id.get((prod, cons))
    if sid is not None:
        sh = shops[sid]
        q_proxy = min(sh.q0_hat, sh.q1_hat)
        if q_proxy > best[2]:
            best = (sid, sid, q_proxy, None)

    # indirect via m
    for m_good in range(n):
        if m_good==prod or m_good==cons: continue
        sid_out = pair_to_shop_id.get((prod, m_good))
        sid_src = pair_to_shop_id.get((m_good, cons))
        if sid_out is None or sid_src is None: continue
        sh_o = shops[sid_out]; sh_s = shops[sid_src]
        q_proxy = min(sh_o.q0_hat, sh_o.q1_hat, sh_s.q0_hat, sh_s.q1_hat)
        if q_proxy > best[2]:
            best = (sid_out, sid_src, q_proxy, m_good)
    return best

def market_research_survey(transactors, g, role, sample_size):
    """
    role: 'sell' or 'buy' for good g.
    Return rough available counts among a random sample.
    """
    pool = []
    if role=='sell':
        pool = [a for a in transactors if a.prod==g]
    else:
        pool = [a for a in transactors if a.cons==g]
    if not pool: return 0
    k = min(sample_size, len(pool))
    # Use PCG32 for sampling instead of random.sample
    sample = []
    pool_copy = pool.copy()
    for _ in range(k):
        idx = RANDOM_N(len(pool_copy))
        sample.append(pool_copy.pop(idx))
    
    # Count those currently *not* matched on that role as a proxy for latent demand/supply
    if role=='sell':
        cnt = sum(1 for a in sample if a.outlet is None)
    else:
        cnt = sum(1 for a in sample if a.source is None)
    # scale back up from sample to population
    return int(round(cnt * (len(pool)/k)))

def consider_entry(params, transactors, shops, available_shop_ids):
    """Propose (prod, cons) using owner's types; open only if expected surplus >= 0."""
    if len(shops) >= params.K or not available_shop_ids: return
    candidates = [r for r in transactors if r.owned_shop_id is None]
    if not candidates: return
    r_idx = RANDOM_N(len(candidates))
    r = candidates[r_idx]
    g0, g1 = r.prod, r.cons
    if (g0, g1) in pair_to_shop_id: return

    new_id = available_shop_ids.pop(0)
    tr0 = RANDOM_N(params.xMax) + 1
    tr1 = RANDOM_N(params.xMax) + 1
    sh = Shop(new_id, g0, g1, r.id, tr0, tr1, params.lam, params.s_overhead)

    # survey both sides to initialize expectations
    sellers_g0 = market_research_survey(transactors, g0, 'sell', params.entry_sample)
    buyers_g0  = market_research_survey(transactors, g0, 'buy',  params.entry_sample)
    sellers_g1 = market_research_survey(transactors, g1, 'sell', params.entry_sample)
    buyers_g1  = market_research_survey(transactors, g1, 'buy',  params.entry_sample)
    sh.q0_hat = max(1.0, float(min(sellers_g0, buyers_g0)))
    sh.q1_hat = max(1.0, float(min(sellers_g1, buyers_g1)))
    sh.set_prices_full_cost()

    y0_exp = sh.p10 * sh.q1_hat
    y1_exp = sh.p01 * sh.q0_hat
    f0 = overhead(g0, params.s_overhead)
    f1 = overhead(g1, params.s_overhead)
    surplus_exp = (y0_exp + y1_exp) - (f0 + f1)
    if surplus_exp >= 0.0:
        shops[new_id] = sh
        pair_to_shop_id[(g0, g1)] = new_id
        r.owned_shop_id = new_id
        sh.setup_paid = True
    else:
        available_shop_ids.insert(0, new_id)

def consider_exit(params, shops, transactors, available_shop_ids):
    to_close=[]
    for sid, sh in shops.items():
        f0 = overhead(sh.g0, params.s_overhead)
        f1 = overhead(sh.g1, params.s_overhead)
        surplus = (sh.y0 + sh.y1) - (f0 + f1)
        if surplus < 0.0 and rnd01() <= params.exit_rate:
            to_close.append(sid)
    for sid in to_close:
        owner = shops[sid].owner
        transactors[owner].owned_shop_id = None
        pair_to_shop_id.pop((shops[sid].g0, shops[sid].g1), None)
        del shops[sid]
        available_shop_ids.append(sid)

def weekly_clear(shops):
    total_y0 = total_y1 = 0.0
    for sh in shops.values():
        S0 = len(sh.customers_sell_g0); B0 = len(sh.customers_buy_g0)
        S1 = len(sh.customers_sell_g1); B1 = len(sh.customers_buy_g1)
        q0_real = float(min(S0, B0))
        q1_real = float(min(S1, B1))
        y0 = sh.p10 * q1_real
        y1 = sh.p01 * q0_real
        sh.update_targets(y0, y1)
        sh.q0_hat = max(q0_real, 1.0)
        sh.q1_hat = max(q1_real, 1.0)
        total_y0 += y0; total_y1 += y1
    return total_y0, total_y1

def compute_money_share(intermediary_counts):
    total = sum(intermediary_counts.values())
    if total<=0: return (None, 0.0)
    g, cnt = max(intermediary_counts.items(), key=lambda kv: kv[1])
    return (g, cnt/total)

def compute_c_stats(transactors, shops, n):
    part = 0
    moneytraders = 0
    usingmoney = [0] * (n+1)
    for a in transactors:
        if a.outlet is not None:
            sh_a = shops.get(a.outlet)
            sh_b = shops.get(a.source) if a.source is not None else None
            if sh_a is None: continue
            ma = 1 if sh_a.g0 == a.prod else 0
            mb = 1 if sh_b and sh_b.g0 == a.cons else 0
            g_ma_a = sh_a.g1 if ma else sh_a.g0
            g_mb_b = sh_b.g1 if mb else sh_b.g0 if sh_b else None
            if g_ma_a == a.cons or (g_mb_b is not None and g_ma_a == g_mb_b):
                part += 1
            if g_mb_b is not None and g_ma_a == g_mb_b:
                moneytraders += 1
                usingmoney[g_mb_b] += 1
    return part, moneytraders, usingmoney

# =========================================
# Stationary benchmark (compact heuristic)
# =========================================
def stationary_benchmark(params, last_window_stats):
    """
    A compact comparator for the stationary commodity-money case:
      - choose money m* as the lowest-overhead good
      - shops: (i, m*) for all i != m*
      - use last-window average throughput per active shop as q-hat
      - full-cost prices that exactly cover overheads and steady incomes
    This is a *pragmatic* benchmark to measure how close the run is to a
    single-money structure; it is not a closed-form reproduction from the paper.
    """
    n = params.n
    m_star = min(range(n), key=lambda g: overhead(g, params.s_overhead))
    # use recent averages
    avg_q0 = max(last_window_stats.get("avg_q0_hat", 1.0), 1.0)
    avg_q1 = max(last_window_stats.get("avg_q1_hat", 1.0), 1.0)
    f_money = overhead(m_star, params.s_overhead)

    # For each non-money good i, there is one shop (i, m*)
    # Prices that would cover overheads at those throughputs:
    # p01 ~ f(m*)/avg_q0 ; p10 ~ f(i)/avg_q1
    # GDP_bench ~ sum positive surpluses; here ~ 0 by construction (covering overheads)
    prices = []
    gdp_bench = 0.0
    for i in range(n):
        if i==m_star: continue
        f_i = overhead(i, params.s_overhead)
        p01 = f_money / avg_q0  # price of i in m*
        p10 = f_i     / avg_q1  # price of m* in i
        prices.extend([p01, p10])

    price_sigma = np.std(prices) if prices else 0.0
    return {
        "money_star": m_star,
        "price_sigma": price_sigma,
        "gdp_bench": gdp_bench
    }

# =========================================
# Single run
# =========================================
def run_once(run_idx, params):
    global pair_to_shop_id
    seed_used = seed_randomness(RAND + run_idx - 1)
    start_time = time.time()

    transactors = create_transactors(params.n, params.b)
    shops = {}
    pair_to_shop_id = {}
    available_shop_ids = list(range(params.K))

    # warm start with a few shops
    tries = 0
    while sum(1 for _ in shops) < min(params.warm_shops, params.K) and tries < params.K*2:
        r_idx = RANDOM_N(len(transactors))
        r = transactors[r_idx]
        if (r.prod, r.cons) in pair_to_shop_id:
            tries += 1; continue
        sid = available_shop_ids.pop(0)
        tr0 = RANDOM_N(params.xMax) + 1
        tr1 = RANDOM_N(params.xMax) + 1
        sh = Shop(sid, r.prod, r.cons, r.id, tr0, tr1, params.lam, params.s_overhead)
        sh.q0_hat = 1.0; sh.q1_hat = 1.0
        sh.set_prices_full_cost()
        sh.setup_paid = True
        shops[sid] = sh
        pair_to_shop_id[(r.prod, r.cons)] = sid
        r.owned_shop_id = sid

    intermediary_counts = Counter()
    gdp_series=[]; shops_series=[]; price_sigma_series=[]
    q0h_series=[]; q1h_series=[]
    connect_both_series=[]

    if not params.quiet:
        print("Number  Using  Using  Using  Using  Using  Using                ")
        print("Active  Money  good1  good2  good3  good4  good5   Year   NS   ")

    def log_progress(week):
        if params.quiet: return
        # if (week+1) % 50 == 0:
        if True:
            part, moneytraders, usingmoney = compute_c_stats(transactors, shops, params.n)
            year = (week+1) // 50
            ns = len(shops)
            print(f"{part:6.0f} {moneytraders:6.0f} {usingmoney[1]:6.0f} {usingmoney[2]:6.0f} {usingmoney[3]:6.0f} {usingmoney[4]:6.0f} {usingmoney[5]:6.0f} {year:6d} {ns:4d}")
            exit(0)

    # main loop
    for t in range(params.T):
        # Entry
        consider_entry(params, transactors, shops, available_shop_ids)

        # Reset weekly customer sets
        for sh in shops.values(): sh.weekly_reset()

        # Shopping/search
        # Create a random permutation like the C version's lineup()
        line = list(range(params.m))
        # Fisher-Yates shuffle using PCG32 to match C version
        for j in range(params.m):
            k = RANDOM_N(params.m - j) + j
            if k < len(line):
                line[j], line[k] = line[k], line[j]
        
        for idx in line:
            a = transactors[idx]
            if a.owned_shop_id is not None:  # owners don't shop as customers
                continue
            needs_link = (a.outlet is None or a.source is None)
            if needs_link or (rnd01() < params.search):
                best_outlet, best_source, qstar, interm = best_path_for_agent(a, shops, params.n)
                curr_q = 0.0
                if a.outlet in shops and a.source in shops:
                    sh_o = shops[a.outlet]; sh_s = shops[a.source]
                    curr_q = min(sh_o.q0_hat, sh_o.q1_hat, sh_s.q0_hat, sh_s.q1_hat)
                if best_outlet is not None and best_source is not None and qstar > curr_q:
                    a.outlet = best_outlet; a.source = best_source
                    if interm is not None:
                        intermediary_counts[interm] += 1

            # enroll
            if a.outlet in shops:
                sh = shops[a.outlet]
                if sh.g0 == a.prod and sh.g1 != a.prod:
                    sh.customers_sell_g0.add(a.id)
                elif sh.g1 == a.prod and sh.g0 != a.prod:
                    sh.customers_sell_g1.add(a.id)
            if a.source in shops:
                sh = shops[a.source]
                if sh.g0 == a.cons:
                    sh.customers_buy_g0.add(a.id)
                elif sh.g1 == a.cons:
                    sh.customers_buy_g1.add(a.id)

        # Clear markets
        total_y0, total_y1 = weekly_clear(shops)

        # Exit
        consider_exit(params, shops, transactors, available_shop_ids)

        # Repost prices
        for sh in shops.values():
            sh.set_prices_full_cost()

        # Metrics
        prices=[]
        gdp=0.0
        q0h_acc=0.0; q1h_acc=0.0; shcount=0
        for sh in shops.values():
            f0 = overhead(sh.g0, params.s_overhead)
            f1 = overhead(sh.g1, params.s_overhead)
            surplus=(sh.y0+sh.y1)-(f0+f1)
            gdp += max(surplus, 0.0)
            prices.extend([sh.p01, sh.p10])
            q0h_acc += sh.q0_hat; q1h_acc += sh.q1_hat; shcount+=1
        sigma = (np.std(prices) if prices else 0.0)
        gdp_series.append(gdp)
        shops_series.append(len(shops))
        price_sigma_series.append(sigma)
        if shcount>0:
            q0h_series.append(q0h_acc/shcount)
            q1h_series.append(q1h_acc/shcount)
        else:
            q0h_series.append(1.0); q1h_series.append(1.0)

        # connectivity: share of agents with both links
        have_both = sum(1 for a in transactors if a.outlet in shops and a.source in shops)
        connect_both_series.append(have_both/params.m)

        log_progress(t)

    # Last-window averages for benchmark
    W = max(1, params.T//10)  # last 10%
    last = slice(-W, None)
    last_stats = {
        "avg_q0_hat": float(np.mean(q0h_series[last])) if q0h_series else 1.0,
        "avg_q1_hat": float(np.mean(q1h_series[last])) if q1h_series else 1.0,
        "avg_price_sigma": float(np.mean(price_sigma_series[last])) if price_sigma_series else 0.0,
        "avg_shops": float(np.mean(shops_series[last])) if shops_series else 0.0,
        "avg_gdp": float(np.mean(gdp_series[last])) if gdp_series else 0.0,
        "avg_connect": float(np.mean(connect_both_series[last])) if connect_both_series else 0.0,
    }
    bench = stationary_benchmark(params, {
        "avg_q0_hat": last_stats["avg_q0_hat"],
        "avg_q1_hat": last_stats["avg_q1_hat"]
    })

    # money winner (by chosen intermediary share)
    money_good, money_share = compute_money_share(intermediary_counts)

    # A few summary flags (paper-style)
    market_developed = (last_stats["avg_connect"] > 0.80 and last_stats["avg_shops"] >= params.n-2)
    monetary_exchange = (money_good is not None and money_share >= 0.60)

    summary = {
        "seed": seed_used,
        "T": params.T,
        "n": params.n,
        "s_overhead": params.s_overhead,
        "shops_final": shops_series[-1] if shops_series else 0,
        "gdp_final": gdp_series[-1] if gdp_series else 0.0,
        "avg_gdp_last10p": last_stats["avg_gdp"],
        "avg_shops_last10p": last_stats["avg_shops"],
        "avg_price_sigma_last10p": last_stats["avg_price_sigma"],
        "avg_connect_last10p": last_stats["avg_connect"],
        "money_good": money_good,
        "money_share": money_share,
        "market_developed": market_developed,
        "monetary_exchange": monetary_exchange,
        "bench_money_star": bench["money_star"],
        "bench_price_sigma": bench["price_sigma"],
        "bench_gdp": bench["gdp_bench"],
    }
    elapsed = time.time() - start_time
    print(f"Run number {run_idx+1}. Time elapsed: {elapsed:.2f} seconds.")
    print(f"Slope equals {params.s_overhead:.0f}, xMax equals {params.xMax}")
    print()
    return summary

# =========================================
# Multi-run harness
# =========================================
def run_many(args):
    params = Params(
        T=args.T, n=args.n, b=args.b, K=args.K, xMax=args.xMax,
        search=args.search, lam=args.lam, exit_rate=args.exit,
        Csetup=DEF_Csetup, s_overhead=args.s,
        warm_shops=args.warm_shops, entry_sample=args.entry_sample,
        log_every=args.log_every, quiet=args.quiet,
    )
    results=[]
    for r in range(args.runs):
        if not args.quiet:
            print(f"\n=== RUN {r+1}/{args.runs} (T={params.T}, s={params.s_overhead}) ===")
        summary = run_once(r + 1, params)  # C uses 1-based run indexing
        results.append(summary)

    # Aggregate like a paper table
    dev_rate = np.mean([1.0 if s["market_developed"] else 0.0 for s in results])
    mon_rate = np.mean([1.0 if s["monetary_exchange"] else 0.0 for s in results])
    money_hist = Counter(s["money_good"] for s in results if s["money_good"] is not None)
    avg_shops = np.mean([s["avg_shops_last10p"] for s in results])
    avg_gdp   = np.mean([s["avg_gdp_last10p"] for s in results])
    avg_sigma = np.mean([s["avg_price_sigma_last10p"] for s in results])

    print("\n====== SUMMARY (multi-run) ======")
    print(f"Runs: {args.runs} | T: {params.T} | n: {params.n} | s: {params.s_overhead}")
    print(f"Market development rate: {dev_rate:.1%}")
    print(f"Monetary-exchange rate: {mon_rate:.1%}")
    if money_hist:
        top = money_hist.most_common(3)
        pretty = ", ".join([f"g{g}:{c}" for g,c in top])
        print(f"Money winners (count): {pretty}")
    else:
        print("Money winners: none")
    print(f"Avg shops (last 10%): {avg_shops:.2f}")
    print(f"Avg GDP   (last 10%): {avg_gdp:.2f}")
    print(f"Avg price σ (last 10%): {avg_sigma:.4f}")

    # Print each run succinctly
    print("\nPer-run:")
    for i,s in enumerate(results):
        mg = "-" if s["money_good"] is None else f"g{s['money_good']} ({s['money_share']:.0%})"
        print(f"  #{i+1}: dev={int(s['market_developed'])} mon={int(s['monetary_exchange'])} "
              f"shops={s['avg_shops_last10p']:.1f} gdp={s['avg_gdp_last10p']:.1f} "
              f"money={mg} bench_m*=g{s['bench_money_star']}")

if __name__ == "__main__":
    args = parse_args()
    run_many(args)
