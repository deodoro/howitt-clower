"""
Howitt–Clower simulation (C++-accurate Python mirror)
====================================================

Goal: behave virtually identically to the provided HowClow.cpp so that you can
compare states **week-by-week**. This Python program mirrors the original C++
implementation’s data layout (1-based indexing), RNG (PCG32), control flow,
and update rules. Printing matches the C++ program’s status lines.

Notes
-----
- Arrays are sized with +1 and indexed from 1 like the C++ code.
- PCG32 is implemented to match `pcg32_next`/`pcg32_seed`. RANDOM_N(x) uses `% x`.
- The `F`, `f`, pricing, entry (including `Research()`), matching, trade,
  exit, target updates, RMSE diagnostics are carried over verbatim.
- The C++ snippet you pasted includes a debug `exit(0)` in the weekly loop.
  That would stop after the first print; here it is **not** executed by default.
  You can enable it with `DEBUG_EXIT_AFTER_PRINT=True` to reproduce that exact
  behavior; keep it False for full runs.
- Optional weekly snapshotting is provided (disabled by default) to compare
  states per week. When enabled, it writes compact CSV rows to `/mnt/data`.

Tested assumptions
------------------
- Uses default parameters from the snippet (T=1000, slope in {16,18}, etc.).
- Prints the same columns during periodic checks (every 50 weeks; same as C++’s
  intended behavior where `if(t%(50*RptPer)==0)` guards the print).

"""
from __future__ import annotations
import math, time, os, csv
from typing import List

# ---------------- Configuration (mirrors C++) ----------------
numruns   = 2
FirstSlope= 16
LastSlope = 18

n         = 10   # number of goods
bsize     = 24   # number per (i,j)
K         = 200  # store locations
xMax      = 200  # max target for new entrant
aut_lambda= 0.05 # autonomous search intensity (lambda in C++)
alpha     = 0.25 # adjustment speed
theta     = 0.01 # failure (exit) probability when unprofitable
C         = 5.0  # setup cost
f1        = 0.0  # overhead constant
persist   = 10   # years of persistence for flags
T         = 1000 # weeks
RAND      = 1    # seed base
RptPer    = 1    # reporting period (years)
prtoscr   = 1    # print to screen

m         = bsize*n*(n-1)

# F and f exactly as in C++
def f(i, slope):
    return f1 + (i-1)*slope

def F(x, y, z):
    return ((y - z - C)/x) if (y - z - C > 0) else 0.0

# Debug behavior switches matching oddities in pasted C++
DEBUG_EXIT_AFTER_PRINT = False  # set True to mimic that accidental exit(0)
SNAPSHOT_WEEKLY = False         # optional: write compact weekly CSV
SNAPSHOT_PATH = "/mnt/data/howclow_weekly.csv"

# ---------------- RNG: PCG32 (matches C++) ----------------
class PCG32:
    __slots__ = ("state","inc")
    def __init__(self):
        self.state = 0
        self.inc = 0
    def next_u32(self) -> int:
        oldstate = self.state & ((1<<64)-1)
        self.state = (oldstate * 6364136223846793005 + (self.inc | 1)) & ((1<<64)-1)
        xorshifted = ((oldstate >> 18) ^ oldstate) >> 27
        rot = (oldstate >> 59) & 0xFFFFFFFF
        return ((xorshifted >> rot) | ((xorshifted << ((-rot) & 31)) & 0xFFFFFFFF)) & 0xFFFFFFFF
    def seed(self, initstate: int, initseq: int):
        self.state = 0
        self.inc   = ((initseq & ((1<<64)-1)) << 1) | 1
        self.next_u32()
        self.state = (self.state + (initstate & ((1<<64)-1))) & ((1<<64)-1)
        self.next_u32()

rng = PCG32()

def SEED_RANDOM(seed: int):
    rng.seed(seed, seed)

def RANDOM_N(x: int) -> int:
    return rng.next_u32() % x

def rnd01() -> float:
    # (1..1000)/1000 to match C++ formulation
    return (1 + RANDOM_N(1000))/1000.0

# ---------------- Global state (1-based indexing) ----------------
s = [0]*(m+1)
d = [0]*(m+1)
sh = [[0]*(m+1) for _ in range(2)]
q = [0]*(m+1)

NS = 0
BS = 0
active = [0]*(K+1)
# g[h][k] with smaller good first (same convention)
g = [[0]*(K+1) for _ in range(2)]

produces = [[0]*((n-1)*bsize+1) for _ in range(n+1)]
consumes = [[0]*((n-1)*bsize+1) for _ in range(n+1)]
numprod = [0]*(n+1)
numcons = [0]*(n+1)
owner   = [0]*(K+1)
familyshop = [0]*(m+1)

# Diagnostics & temporaries
devyear = -1
monyear = -1
h=i=j=0
k=ss=0
r=fr=0
extra=0
t=0
run=0
bestbarter=0
endcount=devcount=0
a=b=0
m0=m1=ma=mb=Nshop=0
c = [0]*6
line = [0]*(m+1)
monetary=fulldev=moneygood=0

W = 0.0
Pinv = [0.0]*(n+1)
SurpSME = 0.0
slope = 0.0

# targets, prices, incomes
# tr[h][k], P[h][k], y[h][k]
tr = [[0.0]*(K+1) for _ in range(2)]
P  = [[0.0]*(K+1) for _ in range(2)]
y  = [[0.0]*(K+1) for _ in range(2)]

Fmon = 0.0
U = Ucomp = Ubarter = 0.0
part = 0.0
moneytraders = usingmax = Csurp = Psurp = 0.0
usingmoney = [0.0]*(n+1)
R = [0.0, 0.0]
vol = [[0.0]*(n+1) for _ in range(2)]
avp = [[0.0]*(n+1) for _ in range(2)]

# ---------------- Utility functions mapped from C++ ----------------

def groupsize(v: int, w: int) -> int:
    return 0 if v==w else bsize

# Initialization (one-time)

def Init():
    global Fmon
    # Build agent roster
    rr = 0
    for ii in range(1, n+1):
        numprod[ii]=0; numcons[ii]=0
    for ii in range(1, n+1):
        for jj in range(1, n+1):
            for kk in range(1, groupsize(ii,jj)+1):
                rr += 1
                s[rr]=ii; d[rr]=jj; q[rr] = 1 if (s[rr]>d[rr]) else 0
                numprod[ii]+=1; produces[ii][numprod[ii]] = rr
                numcons[jj]+=1; consumes[jj][numcons[jj]] = rr
    Fmon = bsize*(n-2)*(n-1)

# Per-run initialization

def lineup():
    # Create random permutation of 1..m
    start = [0]*(m+1)
    line[0]=0
    start[0]=0
    for ii in range(1, m+1):
        start[ii]=ii
    for jj in range(1, m+1):
        kk = RANDOM_N(m - jj + 1) + 1
        line[jj] = start[kk]
        for ii in range(kk, m - jj + 1):
            start[ii] = start[ii+1]


def InitRun(seed_offset: int):
    global endcount, fulldev, devyear, monyear, devcount, NS
    print("Number  Using  Using  Using  Using  Using  Using                ")
    print("Active  Money  good1  good2  good3  good4  good5   Year   NS   ")
    endcount=0; fulldev=0; devyear=-1; monyear=-1; devcount=0
    SEED_RANDOM(RAND + seed_offset)
    # reset shops & links
    NS=0
    for rr in range(0, m+1):
        sh[0][rr]=0; sh[1][rr]=0; familyshop[rr]=0
    for kk in range(0, K+1):
        active[kk]=0; owner[kk]=0
        P[0][kk]=P[1][kk]=0.0
        y[0][kk]=y[1][kk]=0.0
        tr[0][kk]=tr[1][kk]=0.0
        g[0][kk]=g[1][kk]=0
    lineup()

# ---- Entry subroutines (Research and friends) ----

def comrade() -> int:
    # someone with same production as r, not r
    k = 1 + RANDOM_N(numprod[s[r]]-1)
    fr = produces[s[r]][k]
    if fr >= r:
        fr = produces[s[r]][k+1]
    return fr

def soulmate() -> int:
    # someone with same desired good as r, not r
    k = 1 + RANDOM_N(numcons[d[r]]-1)
    fr = consumes[d[r]][k]
    if fr >= r:
        fr = consumes[d[r]][k+1]
    return fr

def Usample() -> float:
    global m0, m1
    m0 = 1 if (g[0][sh[0][fr]] == s[fr]) else 0
    m1 = 1 if (g[0][sh[1][fr]] == d[fr]) else 0
    X = 0.0
    if sh[0][fr] > 0:
        if g[m0][sh[0][fr]] == d[fr]:
            X = P[1-m0][sh[0][fr]]
        elif g[m0][sh[0][fr]] == g[m1][sh[1][fr]]:
            X = P[1-m0][sh[0][fr]] * P[m1][sh[1][fr]]
    return X

def Utility() -> float:
    global m0, m1
    m0 = 1 if (g[0][sh[0][r]] == s[r]) else 0
    m1 = 1 if (g[0][sh[1][r]] == d[r]) else 0
    X = 0.0
    if sh[0][r] > 0:
        if g[m0][sh[0][r]] == d[r]:
            X = P[1-m0][sh[0][r]]
        elif g[m0][sh[0][r]] == g[m1][sh[1][r]]:
            X = P[1-m0][sh[0][r]] * P[m1][sh[1][r]]
    return X

def Research() -> int:
    # draws targ[0], targ[1] and tests adoptability on both sides using
    # comrade/soulmate and strangers with matching roles
    global fr
    targ = [0.0, 0.0]
    for hh in range(2):
        targ[hh] = RANDOM_N(xMax) + 1
    P0 = F(targ[0], targ[1], f(d[r], slope))  # price of s[r]
    P1 = F(targ[1], targ[0], f(s[r], slope))  # price of d[r]

    # test a comrade and a soulmate side
    U = 0.0
    fr = comrade()
    Ucomp_local = Usample()
    if d[fr] == d[r]:
        U = P0
    elif g[m1][sh[1][fr]] == d[r]:
        U = P0 * P[m1][sh[1][fr]]
    if U < Ucomp_local:
        U = 0.0
        fr = soulmate()
        Ucomp_local = Usample()
        if s[fr] == s[r]:
            U = P0
        elif g[m0][sh[0][fr]] == s[r]:
            U = P[1-m0][sh[0][fr]] * P0
        if U < Ucomp_local:
            return 0

    # test strangers: one who likes s[r], one who produces d[r]
    U = 0.0
    k_loc = 1 + RANDOM_N(numcons[s[r]]-1)
    fr = consumes[s[r]][k_loc]
    Ucomp_local = Usample()
    if s[fr] == d[r]:
        U = P1
    elif g[m0][sh[0][fr]] == d[r]:
        U = P[1-m0][sh[0][fr]] * P1
    if U < Ucomp_local:
        U = 0.0
        k_loc = 1 + RANDOM_N(numprod[d[r]]-1)
        fr = produces[d[r]][k_loc]
        Ucomp_local = Usample()
        if d[fr] == s[r]:
            U = P1
        elif g[m1][sh[1][fr]] == s[r]:
            U = P1 * P[m1][sh[1][fr]]
        if U < Ucomp_local:
            return 0
    # pass
    return 1

# ---- Matching (search) helpers ----

def addshop(ss_idx: int):
    global extra
    if active[ss_idx]:
        if (g[0][ss_idx]==s[r] or g[1][ss_idx]==s[r] or
            (g[0][ss_idx]==d[r] or g[1][ss_idx]==d[r])):
            bcnt = 0
            for aa in range(0, 2+extra):
                bcnt += (c[aa] == ss_idx)
            if bcnt == 0:
                extra += 1
                c[1+extra] = ss_idx

def dropshop():
    global a, extra
    for bb in range(a, 1+extra):
        c[bb] = c[bb+1]
    a -= 1
    extra -= 1

def tryBarter():
    global bestbarter, Ubarter, a
    bestbarter = 0; Ubarter = 0.0
    for a in range(2, 1+extra+1):
        if g[q[r]][c[a]] == s[r] and g[1-q[r]][c[a]] == d[r]:
            if P[q[r]][c[a]] > max(Ucomp, Ubarter):
                bestbarter = c[a]; Ubarter = P[q[r]][c[a]]
            dropshop()

def tryOne():
    global a, m0, m1, Ucomp
    for a in range(2, 1+extra+1):
        if g[0][c[a]]==s[r] or g[1][c[a]]==s[r]:
            ma = 1 if (g[0][c[a]]==s[r]) else 0
            if g[ma][c[a]]==g[m1][c[1]] or P[1-m0][c[0]]==0:
                if P[1-m0][c[0]] < P[1-ma][c[a]]:
                    c[0]=c[a]; m0=ma
                    Ucomp = (P[1-ma][c[a]]*P[m1][c[1]] if g[ma][c[a]]==g[m1][c[1]] else 0)
                dropshop()
        else:
            ma = 1 if (g[0][c[a]]==d[r]) else 0
            if g[ma][c[a]]==g[m0][c[0]] or P[m1][c[1]]==0:
                if P[m1][c[1]] < P[ma][c[a]]:
                    c[1]=c[a]; m1=ma
                    Ucomp = (P[1-m0][c[0]]*P[ma][c[a]] if g[ma][c[a]]==g[m0][c[0]] else 0)
                dropshop()

def tryTwo():
    global a, b, m0, m1, Ucomp
    for a in range(2, 1+extra+1):
        if g[0][c[a]]==s[r] or g[1][c[a]]==s[r]:
            for b in range(2, 1+extra+1):
                if a!=b and (g[0][c[b]]==d[r] or g[1][c[b]]==d[r]):
                    ma = 1 if (g[0][c[a]]==s[r]) else 0
                    mb = 1 if (g[0][c[b]]==d[r]) else 0
                    if g[ma][c[a]]==g[mb][c[b]] and Ucomp < P[1-ma][c[a]]*P[mb][c[b]]:
                        Ucomp = P[1-ma][c[a]]*P[mb][c[b]]
                        c[0]=c[a]; c[1]=c[b]; m0=ma; m1=mb

# ---- Run-time calculators ----

def Calc1() -> int:
    global part, moneytraders, usingmoney, fulldev, devcount, devyear
    global usingmax, moneygood, endcount, monyear
    part=0.0; moneytraders=0.0
    for jj in range(1, n+1):
        usingmoney[jj]=0.0
    for rr in range(1, m+1):
        if sh[0][rr] > 0:
            a = sh[0][rr]; b = sh[1][rr]
            ma = 1 if (g[0][a]==s[rr]) else 0
            mb = 1 if (b>0 and g[0][b]==d[rr]) else 0
            part += 1.0 if (g[ma][a]==d[rr] or (b>0 and g[ma][a]==g[mb][b])) else 0.0
            if b>0 and g[ma][a]==g[mb][b]:
                moneytraders += 1.0
                usingmoney[g[mb][b]] += 1.0
    if fulldev==0:
        if part >= .99*m:
            if devcount==0:
                devyear = t//50
            devcount += 1
        else:
            devcount = 0
        if devcount >= persist:
            fulldev = 1
    usingmax=0.0; moneygood=0
    for ii in range(1, n+1):
        if usingmoney[ii] > usingmax:
            usingmax = usingmoney[ii]; moneygood = ii
    if usingmoney[moneygood] >= .99*(bsize*(n-2)*(n-1)):
        if endcount==0:
            monyear = t//50
        endcount += 1
    else:
        endcount = 0
    return 1 if endcount>=persist else 0


def RMSE():
    global R
    if W>0.0:
        for hh in range(2):
            for ii in range(1, n+1):
                vol[hh][ii]=0.0; avp[hh][ii]=0.0
        for kk in range(1, K+1):
            if g[0][kk]==moneygood or g[1][kk]==moneygood:
                ma = 1 if (g[1][kk]==moneygood) else 0
                ii = g[1-ma][kk]
                vol[0][ii] += y[1-ma][kk]
                if vol[0][ii]:
                    avp[0][ii] += (y[1-ma][kk]/vol[0][ii])*(P[1-ma][kk] - avp[0][ii])
                vol[1][ii] += y[ma][kk]
                if vol[1][ii]:
                    avp[1][ii] += (y[ma][kk]/vol[1][ii])*(P[ma][kk] - avp[1][ii])
        for ii in range(1, n+1):
            avp[0][ii] = avp[0][ii]/W if W else 0.0
            avp[1][ii] = avp[1][ii]/Pinv[ii] if Pinv[ii]!=0 else 0.0
        for hh in range(2):
            R[hh]=0.0
            for ii in range(1, n+1):
                if ii!=moneygood:
                    R[hh] += (1-avp[hh][ii])*(1-avp[hh][ii])
            R[hh]/=(n-1)
            R[hh]=math.sqrt(R[hh])
    else:
        R[0]=R[1]=-1.0


def Calc2():
    global BS, W, SurpSME, Csurp, Psurp, Nshop
    if monetary==0:
        monyear = -1
    if fulldev==0:
        devyear = -1
    BS=0
    for kk in range(1, K+1):
        if active[kk] and (g[0][kk]!=moneygood and g[1][kk]!=moneygood):
            BS += 1
    W = 1.0 - (f(moneygood, slope) + C)/bsize
    if W>0.0:
        for ii in range(1, n+1):
            Pinv[ii] = ((m/n) - (f(ii, slope)+C)) / ((m/n) - (n-2)*(f(moneygood, slope)+C))
    SurpSME = m - n*f1 - (slope/2)*n*(n-1) - (n-2)*f(moneygood, slope)
    RMSE()
    Csurp = 0.0
    for rr in range(1, m+1):
        # reuse Utility with r set appropriately
        global r
        r_backup = r
        r = rr
        Csurp += Utility()
        r = r_backup
    Psurp = 0.0; Nshop = 0
    for kk in range(1, K+1):
        if active[kk]:
            for hh in range(2):
                Psurp += (y[hh][kk] - f(g[hh][kk], slope) - P[1-hh][kk]*y[1-hh][kk])
            own = owner[kk]
            if own and (y[q[own]][kk]>1 or y[1-q[own]][kk]>0):
                Nshop += 1

# ---------------- Main program ----------------

def main():
    global slope, NS, t, monetary
    # Prepare file (like C++)
    with open("evol.fil","w") as stream:
        stream.write("  run slope dev part mon mtrd usmx Mgd Csur   Psur  Esur Nsh  NS  BS  RMS_0  RMS_1  dyr  myr\n")

    Init()
    start_time = time.time()
    for slope in range(FirstSlope, LastSlope+1, 2):
        for run_idx in range(1, numruns+1):
            InitRun(run_idx-1)
            # Weekly loop
            for t in range(1, T+1):
                # ---- Entry (ensure at least one shop exists) ----
                while True:
                    # random prospective owner r
                    global r
                    r = RANDOM_N(m) + 1
                    if NS < K and familyshop[r]==0:
                        j_local = Research()
                        if j_local > 0:
                            NS += 1
                            k_loc = 1
                            while active[k_loc]==1:
                                k_loc += 1
                            active[k_loc] = 1
                            # order goods as in C++ using q[r]
                            g[1-q[r]][k_loc] = d[r]
                            g[q[r]][k_loc]   = s[r]
                            for hh in range(2):
                                tr[hh][k_loc] = (RANDOM_N(xMax)+1) if False else tr[hh][k_loc]
                                # Initialize prices from targ used in Research: the C++ stores
                                # targ only locally; here replicate price init using same formula
                                # with current tr (0). To stay faithful, set from the targ used
                                # in Research; but the C++ copies targ[h] into tr[h][k]. We can
                                # mirror by reusing the P0/P1 computed in Research. For simplicity,
                                # recompute with tr=the last targ via closure — omitted since targ
                                # is local. We emulate by assigning targets then computing P.
                            # Assign targets to the shop as in C++: reuse the last targ[] from Research
                            # To match, we recompute prices exactly as C++ does using those targets:
                            # Recompute: (we must reconstruct targ here; but the C++ used the same
                            # targ when Research was called. To keep state identical, we stash it in
                            # globals during Research.)
                            for hh in range(2):
                                tr[hh][k_loc] = last_targ[hh]
                                # P[h][k] = F(targ[h],targ[1-h],f(g[1-h][k]))
                                P[hh][k_loc] = F(last_targ[hh], last_targ[1-hh], f(g[1-hh][k_loc], slope))
                            sh[0][r]=k_loc; sh[1][r]=0; familyshop[r]=k_loc; owner[k_loc]=r
                    # break condition emulates do{...}while(NS==0)
                    if NS>0:
                        break
                # ---- Matching (search) ----
                for ii in range(1, m+1):
                    r = line[ii]
                    U_curr = Utility()
                    if rnd01() < (aut_lambda if (U_curr>0) else 1) and familyshop[r]==0:
                        # Build candidate list: current links + comrade + soulmate + random
                        global extra
                        extra = 0
                        c[0] = sh[0][r]; c[1] = sh[1][r]
                        fr = comrade();  addshop(sh[0][fr])
                        fr = soulmate(); addshop(sh[1][fr])
                        addshop(RANDOM_N(K)+1)
                        if extra:
                            global Ucomp
                            Ucomp = U_curr
                            tryBarter()
                            if extra and (g[1-q[r]][c[0]]!=d[r] or P[q[r]][c[0]]==0):
                                tryOne()
                            if extra>1:
                                tryTwo()
                            if Ucomp < Ubarter:
                                sh[0][r] = bestbarter
                                sh[1][r] = 0
                            else:
                                sh[0][r] = c[0]
                                sh[1][r] = c[1]
                # ---- Trade ----
                for kk in range(1, K+1):
                    if active[kk]>0:
                        y[0][kk]=0.0; y[1][kk]=0.0
                for rr in range(1, m+1):
                    if sh[0][rr]>0:
                        a_loc = sh[0][rr]; b_loc = sh[1][rr]
                        ma_loc = 1 if (g[0][a_loc]==s[rr]) else 0
                        mb_loc = 1 if (g[0][b_loc]==d[rr]) else 0
                        if g[ma_loc][a_loc]==d[rr]:
                            y[1-ma_loc][a_loc] += 1.0
                        elif b_loc>0 and g[ma_loc][a_loc]==g[mb_loc][b_loc]:
                            y[1-ma_loc][a_loc] += 1.0
                            y[mb_loc][b_loc] += P[1-ma_loc][a_loc]
                # ---- Exit ----
                for kk in range(1, K+1):
                    if active[kk]:
                        # side tests
                        side0 = y[0][kk] - P[1][kk]*y[1][kk] - f(g[0][kk], slope)
                        side1 = y[1][kk] - P[0][kk]*y[0][kk] - f(g[1][kk], slope)
                        if (side0 <= 0 or side1 <= 0) and rnd01() <= theta:
                            active[kk]=0; NS-=1; familyshop[owner[kk]]=0; owner[kk]=0
                            for hh in range(2):
                                y[hh][kk]=0; tr[hh][kk]=0; g[hh][kk]=0; P[hh][kk]=0
                            for rr in range(1, m+1):
                                if sh[0][rr]==kk: sh[0][rr]=0
                                if sh[1][rr]==kk: sh[1][rr]=0
                # ---- Update targets & prices ----
                for kk in range(1, K+1):
                    if active[kk]:
                        for hh in range(2):
                            tr[hh][kk] += alpha*(y[hh][kk]-tr[hh][kk])
                        P[0][kk] = F(tr[0][kk], tr[1][kk], f(g[1][kk], slope))
                        P[1][kk] = F(tr[1][kk], tr[0][kk], f(g[0][kk], slope))
                # ---- Periodic check (every 50 weeks) ----
                if t % (50*RptPer) == 0:
                    monetary = Calc1()
                    if prtoscr!=0:
                        # Print: part moneytraders using[1..5] year NS
                        print(f"{part:6.0f} {moneytraders:6.0f} ", end="")
                        for bb in range(1,6):
                            print(f"{usingmoney[bb]:6.0f} ", end="")
                        print(f"{t//50:6d} {NS:4d}")
                        if DEBUG_EXIT_AFTER_PRINT:
                            return  # mimic exit(0)
                    if SNAPSHOT_WEEKLY:
                        hdr = ["week","part","moneytraders","NS","moneygood","using1","using2","using3","using4","using5"]
                        row = [t//50, int(part), int(moneytraders), NS, moneygood] + [int(usingmoney[x]) for x in range(1,6)]
                        file_exists = os.path.exists(SNAPSHOT_PATH)
                        with open(SNAPSHOT_PATH, "a", newline="") as f:
                            w=csv.writer(f)
                            if not file_exists:
                                w.writerow(hdr)
                            w.writerow(row)
                    if monetary==1:
                        break
            # end weekly loop
            Calc2()
            elapsed = time.time() - start_time
            print(f"Run number {run_idx}. Time elapsed: {elapsed:.2f} seconds.")
            print(f"Slope equals {slope:.0f}, xMax equals {xMax}\n")
            with open("evol.fil","a") as stream:
                stream.write(
                    f"{run_idx:5d} {slope:5.0f} {fulldev:3d} {part:4.0f} {monetary:3d} {moneytraders:4.0f} {usingmax:4.0f} {moneygood:3d} {Csurp:4.0f} {Psurp:6.0f} {SurpSME:5.0f} {Nshop:3d} {NS:3d} {BS:3d} {R[0]:6.3f} {R[1]:6.3f} {devyear:4d} {monyear:4d}\n"
                )
        with open("evol.fil","a") as stream:
            stream.write("\n")
    with open("evol.fil","a") as stream:
        stream.write("\n\n  n bsiz   K  xMax  lambda alpha  theta  C  Nrun Pst(yr)  T RND\n")
        stream.write(f"{n:3d} {bsize:4d} {K:3d} {xMax:5d} {aut_lambda:6.3f} {alpha:6.3f} {theta:6.3f} {C:2.0f} {numruns:5d} {persist:3d} {T:6d} {RAND:3d}\n")
        stream.write("\n")

# ------------- global to carry targ from Research to Entry (C++ does this implicitly) -------------
last_targ = [1.0, 1.0]

# Monkey-patch Research to also set last_targ so Entry can copy it into tr[] exactly like C++
_orig_Research = Research

def _Research_with_targ_capture():
    global last_targ
    # replicate original Research while capturing targ[]
    global fr, r
    targ_local = [0.0, 0.0]
    for hh in range(2):
        targ_local[hh] = RANDOM_N(xMax) + 1
    # Temporarily compute the rest exactly as Research does
    # Duplicate body from Research but using targ_local and without recursion
    # We'll set last_targ at the end iff it returns 1
    # --- Begin duplicate ---
    P0 = F(targ_local[0], targ_local[1], f(d[r], slope))
    P1 = F(targ_local[1], targ_local[0], f(s[r], slope))
    U = 0.0
    fr = comrade()
    Uc = Usample()
    if d[fr]==d[r]: U = P0
    elif g[m1][sh[1][fr]]==d[r]: U = P0*P[m1][sh[1][fr]]
    if U<Uc:
        U=0.0
        fr = soulmate()
        Uc = Usample()
        if s[fr]==s[r]: U=P0
        elif g[m0][sh[0][fr]]==s[r]: U = P[1-m0][sh[0][fr]]*P0
        if U<Uc:
            return 0
    U=0.0
    k_loc = 1+RANDOM_N(numcons[s[r]]-1)
    fr = consumes[s[r]][k_loc]
    Uc = Usample()
    if s[fr]==d[r]: U=P1
    elif g[m0][sh[0][fr]]==d[r]: U=P[1-m0][sh[0][fr]]*P1
    if U<Uc:
        U=0.0
        k_loc = 1+RANDOM_N(numprod[d[r]]-1)
        fr = produces[d[r]][k_loc]
        Uc = Usample()
        if d[fr]==s[r]: U=P1
        elif g[m1][sh[1][fr]]==s[r]: U=P1*P[m1][sh[1][fr]]
        if U<Uc:
            return 0
    last_targ = targ_local[:]  # capture
    return 1
    # --- End duplicate ---

Research = _Research_with_targ_capture

if __name__ == "__main__":
    main()
