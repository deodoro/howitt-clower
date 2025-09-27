"""
Agent‑based simulation approximating Howitt & Clower (2000),
"The emergence of economic organization" (Journal of Economic Behavior & Organization 41:55–84).

This script follows the paper’s key mechanisms:
- Transactors (type (i,j)) produce 1 unit of i and want to consume j each week.
- Trading only via specialist shops; each shop trades exactly two commodities (g0, g1).
- Transactors can keep at most two trading relationships: outlet (for their production) and source (for their consumption) per the paper.
- Weekly schedule: entry → shopping (search) → exchange/income realization → exit → price/targets update.
- Shop costs: overhead f(i)=s*(i-1), setup/normal return C.
- Pricing: full-cost/target pricing with adaptive targets (α) and posted bid prices p0, p1 (paper §5.4).
- Exit: when any operating surplus turns negative, shop exits with probability u (paper §5.3).
- Search: local sampling via one random shop + a "comrade" (same production) outlet + a "soulmate" (same consumption) source (paper §7.2).
- Entry: rare entrepreneurship with market research on both sides (paper §7.1).

Outputs (stdout):
- Periodic status lines.
- Final summary: GDP, gaps (vs. potential monetary equilibrium using commodity 1), monetary‑exchange emergence stats, price distance metrics.
- CSV snapshots (optional) to /mnt/data/ when SAVE_SNAPSHOTS=True.
"""
from __future__ import annotations
import random
import math
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple
from collections import defaultdict, Counter

# -------------------- Parameters (tune freely) --------------------
SEED = 42
T_WEEKS = 2000           # Fewer than paper's 20,000 for quick runs; increase for richer dynamics
N = 8                    # number of commodities (paper often uses 10)
B = 6                    # transactors per ordered pair type (i,j); m = B*N*(N-1)
K = 120                  # available shop locations (paper uses 200)
xMax = 200               # max animal spirits (target draws)
C = 5.0                  # setup/normal return per commodity
s = 8.0                  # overhead slope: f(i)=s*(i-1). try s in {0,2,4,...,22}
alpha = 0.25             # target adaptation speed
u = 0.01                 # weekly exit probability if any surplus negative
S_search = 0.05          # weekly search probability if already profitable
P_ENTRY = 1.0 / (B*N*(N-1))  # probability a random agent gets an entry idea per week
SAVE_SNAPSHOTS = False
SNAPSHOT_EVERY = 200

random.seed(SEED)

# -------------------- Helpers --------------------

def f_cost(i:int)->float:
    return s * (i-1)

@dataclass
class Shop:
    k: int                     # location id (0..K-1)
    g0: int
    g1: int
    owner: int                 # transactor id
    tr0: float
    tr1: float
    p0: float
    p1: float
    active: bool = True
    # updated weekly
    y0: float = 0.0
    y1: float = 0.0

    def price_update(self):
        # targets adapt towards realized incomes
        self.tr0 += alpha * (self.y0 - self.tr0)
        self.tr1 += alpha * (self.y1 - self.tr1)
        # pricing per paper (section 5.4):
        # p0 = max( (tr1 - f(g1) - C)/tr0, 0 ), p1 symmetric
        self.p0 = 0.0 if self.tr0 <= 1e-12 else max((self.tr1 - f_cost(self.g1) - C) / self.tr0, 0.0)
        self.p1 = 0.0 if self.tr1 <= 1e-12 else max((self.tr0 - f_cost(self.g0) - C) / self.tr1, 0.0)

    def operating_surpluses(self) -> Tuple[float,float]:
        # π0 = y0 - p1*y1 - f(g0); π1 = y1 - p0*y0 - f(g1)
        pi0 = self.y0 - self.p1 * self.y1 - f_cost(self.g0)
        pi1 = self.y1 - self.p0 * self.y0 - f_cost(self.g1)
        return pi0, pi1

@dataclass
class Agent:
    id: int
    prod: int  # i
    cons: int  # j
    outlet: int = -1  # shop id or -1
    source: int = -1  # shop id or -1
    is_owner: bool = False

    def current_consumption(self, shops: List[Optional[Shop]]) -> Tuple[float, Optional[int]]:
        """Return (quantity of desired cons good j, exchange_intermediary or None)
        Based on current outlet/source relationships.
        """
        # Direct: one shop trades both i and j
        if self.outlet >= 0:
            s_out = shops[self.outlet]
            if s_out and s_out.active and ((s_out.g0==self.prod and s_out.g1==self.cons) or (s_out.g1==self.prod and s_out.g0==self.cons)):
                # sells 1 unit of prod; buys at the posted price for prod
                if s_out.g0 == self.prod:
                    p = s_out.p0  # offer price for g0 (prod) is quantity of g1 paid
                else:
                    p = s_out.p1
                return (p, None)
        # Indirect: outlet trades (i,c) and source trades (c,j)
        if self.outlet >= 0 and self.source >= 0:
            s_out = shops[self.outlet]
            s_src = shops[self.source]
            if s_out and s_out.active and s_src and s_src.active:
                # find c common to both shops where outlet trades (prod,c) and source trades (c,cons)
                cands_out = {s_out.g0, s_out.g1}
                cands_src = {s_src.g0, s_src.g1}
                common = list(cands_out.intersection(cands_src))
                # common contains possibly prod or cons; need c != prod and c != cons
                for c in common:
                    if c!=self.prod and c!=self.cons:
                        # compute p (prod→c at outlet) and p' (c→cons at source)
                        if s_out.g0==self.prod and s_out.g1==c:
                            p = s_out.p0
                        elif s_out.g1==self.prod and s_out.g0==c:
                            p = s_out.p1
                        else:
                            continue
                        if s_src.g0==c and s_src.g1==self.cons:
                            pprime = s_src.p0
                        elif s_src.g1==c and s_src.g0==self.cons:
                            pprime = s_src.p1
                        else:
                            continue
                        return (p * pprime, c)
        return (0.0, None)

# -------------------- Economy --------------------
class Economy:
    def __init__(self):
        self.m = B*N*(N-1)
        # Build transactors evenly across ordered pairs (i,j), i!=j
        agents: List[Agent] = []
        aid = 0
        for i in range(1, N+1):
            for j in range(1, N+1):
                if i==j: continue
                for _ in range(B):
                    agents.append(Agent(aid, i, j))
                    aid += 1
        self.agents = agents
        # Shops array of size K; None when vacant
        self.shops: List[Optional[Shop]] = [None]*K
        self.free_slots = set(range(K))
        # fixed weekly permutation for shopping order (per run)
        self.search_order = list(range(self.m))
        random.shuffle(self.search_order)
        self.week = 0
        # Stats
        self.history = {
            'direct_traders': [],
            'indirect_traders': [],
            'money_users': [],
            'num_shops': [],
            'gdp': [],
        }

    # --------------- Utility funcs ---------------
    def vacant_slot(self) -> Optional[int]:
        return random.choice(tuple(self.free_slots)) if self.free_slots else None

    def post_prices_for_pair(self, g0:int, g1:int, tr0:float, tr1:float) -> Tuple[float,float]:
        p0 = 0.0 if tr0<=1e-12 else max((tr1 - f_cost(g1) - C)/tr0, 0.0)
        p1 = 0.0 if tr1<=1e-12 else max((tr0 - f_cost(g0) - C)/tr1, 0.0)
        return p0, p1

    # --------------- Weekly stages ---------------
    def stage_entry(self):
        # One random potential entrepreneur per week with probability P_ENTRY
        if random.random() > P_ENTRY:
            return
        r = random.randrange(self.m)
        ag = self.agents[r]
        if ag.is_owner: return
        if not self.free_slots: return
        # propose shop trading (prod,cons)
        g0, g1 = ag.prod, ag.cons
        # draw initial targets
        tr0 = random.randint(1, xMax)
        tr1 = random.randint(1, xMax)
        p0, p1 = self.post_prices_for_pair(g0,g1,tr0,tr1)
        # market research: sample four prospects (paper §7.1)
        # (i) supplier of g0; (ii) consumer of g1; (iii) supplier of g1; (iv) consumer of g0
        def pick_with(cond):
            cand = [a for a in self.agents if cond(a)]
            return random.choice(cand) if cand else None
        cand_i_prod = pick_with(lambda a: a.prod==g0 and not a.is_owner)
        cand_j_cons = pick_with(lambda a: a.cons==g1 and not a.is_owner)
        cand_j_prod = pick_with(lambda a: a.prod==g1 and not a.is_owner)
        cand_i_cons = pick_with(lambda a: a.cons==g0 and not a.is_owner)
        # build a fake shop instance (not yet placed) to test adoption in samples
        fake = Shop(-1,g0,g1,r,tr0,tr1,p0,p1,active=True)

        def would_adopt(agent:Agent, role:str) -> bool:
            # role is 'supplier' or 'consumer' side; but adoption logic just checks if fake yields higher consumption
            cur_q,_ = agent.current_consumption(self.shops)
            # Build a search sample akin to weekly search: include fake + agent's outlet/source if any
            sample_shop_ids = set()
            if agent.outlet>=0: sample_shop_ids.add(agent.outlet)
            if agent.source>=0: sample_shop_ids.add(agent.source)
            # We'll test candidate sets that include fake as either outlet or source as appropriate
            best_q = cur_q
            # option 1: direct if fake trades (prod,cons)
            if (fake.g0==agent.prod and fake.g1==agent.cons) or (fake.g1==agent.prod and fake.g0==agent.cons):
                # price for selling prod at fake
                p = fake.p0 if fake.g0==agent.prod else fake.p1
                best_q = max(best_q, p)
            # option 2: indirect via common c with existing real shop (outlet or source)
            # try fake as outlet
            for sid in list(sample_shop_ids):
                s = self.shops[sid]
                if not s or not s.active: continue
                # fake (agent.prod, c) and s (c, agent.cons)
                for c in {fake.g0,fake.g1}.intersection({s.g0,s.g1}):
                    if c!=agent.prod and c!=agent.cons:
                        # get p at fake
                        if fake.g0==agent.prod and fake.g1==c:
                            p = fake.p0
                        elif fake.g1==agent.prod and fake.g0==c:
                            p = fake.p1
                        else:
                            continue
                        # get p' at s
                        if s.g0==c and s.g1==agent.cons:
                            pprime = s.p0
                        elif s.g1==c and s.g0==agent.cons:
                            pprime = s.p1
                        else:
                            continue
                        best_q = max(best_q, p*pprime)
            # try fake as source with existing outlet
            for sid in list(sample_shop_ids):
                s = self.shops[sid]
                if not s or not s.active: continue
                for c in {s.g0,s.g1}.intersection({fake.g0,fake.g1}):
                    if c!=agent.prod and c!=agent.cons:
                        if s.g0==agent.prod and s.g1==c:
                            p = s.p0
                        elif s.g1==agent.prod and s.g0==c:
                            p = s.p1
                        else:
                            continue
                        if fake.g0==c and fake.g1==agent.cons:
                            pprime = fake.p0
                        elif fake.g1==c and fake.g0==agent.cons:
                            pprime = fake.p1
                        else:
                            continue
                        best_q = max(best_q, p*pprime)
            return best_q > cur_q + 1e-12

        left_ok = any(would_adopt(c, 'left') for c in [cand_i_prod, cand_j_cons] if c)
        right_ok = any(would_adopt(c, 'right') for c in [cand_j_prod, cand_i_cons] if c)
        if not (left_ok and right_ok):
            return
        # Open the shop
        slot = self.vacant_slot()
        if slot is None: return
        new_shop = Shop(slot,g0,g1,r,tr0=float(tr0),tr1=float(tr1),p0=p0,p1=p1,active=True)
        self.shops[slot] = new_shop
        self.free_slots.remove(slot)
        ag.is_owner = True
        # entrepreneur becomes first customer (paper §7.1)
        # choose best role for them given their type
        # set as direct if possible, else as outlet by default
        if (g0==ag.prod and g1==ag.cons) or (g1==ag.prod and g0==ag.cons):
            ag.outlet = slot
            ag.source = -1
        else:
            # default: make it their outlet
            ag.outlet = slot

    def stage_shopping(self):
        for r in self.search_order:
            ag = self.agents[r]
            if ag.is_owner:  # owners do not search as customers this week (paper §7.2)
                continue
            cur_q, _ = ag.current_consumption(self.shops)
            will_search = (cur_q <= 1e-12) or (random.random() < S_search)
            if not will_search:
                continue
            # build sample
            sample_ids = set()
            # one random shop
            existing = [k for k,s in enumerate(self.shops) if s and s.active]
            if existing:
                sample_ids.add(random.choice(existing))
            # comrade (same production) outlet
            comrades = [a for a in self.agents if a.prod==ag.prod and a.outlet>=0]
            if comrades:
                sample_ids.add(random.choice(comrades).outlet)
            # soulmate (same consumption) source
            soulmates = [a for a in self.agents if a.cons==ag.cons and a.source>=0]
            if soulmates:
                sample_ids.add(random.choice(soulmates).source)
            # evaluate best option from sample
            best_q = cur_q
            best_outlet = ag.outlet
            best_source = ag.source
            # consider direct options
            for sid in list(sample_ids):
                s = self.shops[sid]
                if not s or not s.active: continue
                if (s.g0==ag.prod and s.g1==ag.cons):
                    if s.p0 > best_q + 1e-12:
                        best_q = s.p0; best_outlet = sid; best_source = -1
                if (s.g1==ag.prod and s.g0==ag.cons):
                    if s.p1 > best_q + 1e-12:
                        best_q = s.p1; best_outlet = sid; best_source = -1
            # consider indirect via common c
            for sid_o in list(sample_ids):
                so = self.shops[sid_o]
                if not so or not so.active: continue
                for sid_s in list(sample_ids):
                    ss = self.shops[sid_s]
                    if not ss or not ss.active: continue
                    commons = {so.g0,so.g1}.intersection({ss.g0,ss.g1})
                    for c in commons:
                        if c==ag.prod or c==ag.cons: continue
                        # outlet path prod->c
                        if so.g0==ag.prod and so.g1==c:
                            p = so.p0
                        elif so.g1==ag.prod and so.g0==c:
                            p = so.p1
                        else:
                            continue
                        # source path c->cons
                        if ss.g0==c and ss.g1==ag.cons:
                            pprime = ss.p0
                        elif ss.g1==c and ss.g0==ag.cons:
                            pprime = ss.p1
                        else:
                            continue
                        q = p*pprime
                        if q > best_q + 1e-12:
                            best_q = q; best_outlet = sid_o; best_source = sid_s
            # adopt
            if best_q > cur_q + 1e-12:
                ag.outlet, ag.source = best_outlet, best_source
            elif cur_q <= 1e-12:
                # emergency: try to latch onto any outlet offering positive price for prod
                # in the sample
                adopted = False
                for sid in list(sample_ids):
                    s = self.shops[sid]
                    if not s or not s.active: continue
                    price = None
                    if s.g0==ag.prod: price = s.p0
                    if s.g1==ag.prod: price = s.p1 if price is None else max(price, s.p1)
                    if price and price>1e-12:
                        ag.outlet = sid
                        adopted = True
                        break
                if not adopted:
                    # leave as is (may search next week)
                    pass

    def stage_exchange(self):
        # reset incomes
        for s in self.shops:
            if s: s.y0 = s.y1 = 0.0
        # each agent sells 1 unit at outlet if possible, and (if indirect) buys at source; but incomes y are just inflows per shop
        for ag in self.agents:
            # direct trade at outlet if it trades (prod,cons)
            if ag.outlet>=0:
                s = self.shops[ag.outlet]
                if s and s.active:
                    if s.g0==ag.prod:
                        s.y0 += 1.0
                    elif s.g1==ag.prod:
                        s.y1 += 1.0
        # Note: We do not need to track physical flows for GDP; incomes y are enough for pricing and surplus.

    def stage_exit(self):
        for k,s in enumerate(self.shops):
            if not s or not s.active: continue
            pi0, pi1 = s.operating_surpluses()
            if (pi0 <= 0.0 or pi1 <= 0.0) and random.random() < u:
                # close
                s.active = False
                self.free_slots.add(k)
                # sever relationships
                for ag in self.agents:
                    if ag.outlet == k: ag.outlet = -1
                    if ag.source == k: ag.source = -1

    def stage_pricing(self):
        for s in self.shops:
            if s and s.active:
                s.price_update()

    # --------------- Metrics ---------------
    def agent_consumption_and_money(self) -> Tuple[int,int,Counter]:
        direct = 0
        indirect = 0
        money_counter = Counter()
        for ag in self.agents:
            q,c = ag.current_consumption(self.shops)
            if q>1e-12:
                if c is None:
                    direct += 1
                else:
                    indirect += 1
                    money_counter[c]+=1
        return direct, indirect, money_counter

    def gdp_now(self) -> float:
        # GDP = sum over shops of (y0 + y1 - f(g0) - f(g1))
        gdp = 0.0
        for s in self.shops:
            if s and s.active:
                gdp += s.y0 + s.y1 - f_cost(s.g0) - f_cost(s.g1)
        return gdp

    def potential_gdp_money1(self) -> float:
        # potential GDP in stationary money eq. with commodity 1 as money: m - sum_{i=2..N} f(i)
        return self.m - sum(f_cost(i) for i in range(2,N+1))

    def price_distance_from_stationary(self, money:int) -> Tuple[float,float]:
        # The paper computes closed-form stationary prices; here we approximate distance by RMS across active shops
        # Dist0: wholesale (offer) price distance for non-money goods; Dist1: retail (reciprocal) distances
        # We will compute average p0/p1 for shops trading (j,money) and compare across them.
        ps = []
        invps = []
        for s in self.shops:
            if not s or not s.active: continue
            # if s trades (j,m)
            if money in (s.g0, s.g1):
                if s.g0==money:
                    # p1 is price of j when seller brings j (since p1 is offer price for g1 buyers of g0?)
                    # Keep simple: collect both posted offer prices
                    ps.append(s.p0)
                    ps.append(s.p1)
                    if s.p0>1e-12: invps.append(1.0/s.p0)
                    if s.p1>1e-12: invps.append(1.0/s.p1)
                else:
                    ps.append(s.p0)
                    ps.append(s.p1)
                    if s.p0>1e-12: invps.append(1.0/s.p0)
                    if s.p1>1e-12: invps.append(1.0/s.p1)
        if not ps:
            return 0.0, 0.0
        mean_p = sum(ps)/len(ps)
        mean_inv = sum(invps)/len(invps) if invps else 0.0
        rms_p = math.sqrt(sum((x-mean_p)**2 for x in ps)/len(ps)) if len(ps)>1 else 0.0
        rms_inv = math.sqrt(sum((x-mean_inv)**2 for x in invps)/len(invps)) if len(invps)>1 else 0.0
        return rms_p, rms_inv

    # --------------- Run ---------------
    def run(self):
        pot = self.potential_gdp_money1()
        monetary_emerged = False
        monetary_commodity = None
        max_indirect_theoretical = int(self.m * (N-2)/N)
        for t in range(1, T_WEEKS+1):
            self.week = t
            self.stage_entry()
            self.stage_shopping()
            self.stage_exchange()
            self.stage_exit()
            self.stage_pricing()
            # stats
            direct, indirect, mc = self.agent_consumption_and_money()
            gdp = self.gdp_now()
            active_shops = sum(1 for s in self.shops if s and s.active)
            self.history['direct_traders'].append(direct)
            self.history['indirect_traders'].append(indirect)
            self.history['num_shops'].append(active_shops)
            self.history['gdp'].append(gdp)
            if mc:
                money, users = mc.most_common(1)[0]
                self.history['money_users'].append(users)
                # detect emergence (within 1% of max indirect users for 50 consecutive weeks)
                if users >= 0.99*max_indirect_theoretical:
                    # require persistence over last 50 weeks
                    window = self.history['money_users'][-50:]
                    if len(window)==50 and min(window) >= 0.99*max_indirect_theoretical:
                        monetary_emerged = True
                        monetary_commodity = money
                        # We may continue running to let shops consolidate.
            if SAVE_SNAPSHOTS and t % SNAPSHOT_EVERY == 0:
                try:
                    import csv, os
                    with open(f"/mnt/data/snap_{t}.csv","w",newline="") as f:
                        w = csv.writer(f)
                        w.writerow(["week","direct","indirect","money_users","shops","gdp"])
                        w.writerow([t,direct,indirect,self.history['money_users'][-1] if self.history['money_users'] else 0,active_shops,gdp])
                except Exception:
                    pass
            if t % 200 == 0:
                print(f"Week {t:5d} | shops={active_shops:3d} direct={direct:4d} indirect={indirect:4d} gdp={gdp:7.2f}")
        # final metrics
        direct, indirect, mc = self.agent_consumption_and_money()
        gdp = self.gdp_now()
        gdp_gap = (pot - gdp)/pot if pot>0 else float('nan')
        if mc:
            money, users = mc.most_common(1)[0]
        else:
            money, users = None, 0
        dist0, dist1 = self.price_distance_from_stationary(money) if money else (0.0,0.0)
        print("\n----- FINAL SUMMARY -----")
        print(f"Weeks run: {T_WEEKS}")
        print(f"Active shops: {sum(1 for s in self.shops if s and s.active)} / {K}")
        print(f"Direct traders: {direct}  Indirect traders: {indirect}")
        print(f"Most common intermediary: {money} with users={users} (theoretical max {int(self.m*(N-2)/N)})")
        print(f"Monetary exchange emerged? {monetary_emerged} (money={monetary_commodity})")
        print(f"GDP (this run): {gdp:.2f}  Potential (money=1): {pot:.2f}  Gap: {gdp_gap*100:.2f}%")
        print(f"Price distances (RMS) approx → Dist0={dist0:.4f} Dist1={dist1:.4f}")

if __name__ == "__main__":
    econ = Economy()
    econ.run()
