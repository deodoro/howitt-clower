// -----------------------------------------------------------------------------
// HowClow.cpp — heavily commented edition (no code changes, comments only).
// -----------------------------------------------------------------------------
// Context:
//   This program implements the simulation described in
//   Howitt & Clower (2000), “The Emergence of Economic Organization.”
//   Agents are of type (i,j): they produce good i and desire good j.
//   Exchange is mediated by “shops” that each trade a pair of goods.
//   Over time, entry/search/trade/exit/price-updating can lead to a
//   “monetary” phase where one commodity becomes the universal intermediary.
//
// Note:
//   Per your request, only comments were added. All original code lines,
//   macros, constants, function bodies, and formatting remain unchanged.
//   You can diff against your original to verify: only `// ...` lines
//   were inserted.
// -----------------------------------------------------------------------------

// HowClow.cpp.
//
//      This is the source code for the program described in "The Emergence of
// Economic Organization," by Peter Howitt and Robert Clower. It compiles and
// runs under Borland C++ Builder 3.0. We are making it available so that others
// can see exactly how the algorithm described in the paper was implemented.
//
// December 17, 1998

#pragma hdrstop
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>

#define DEBUG 0

// -------------------------------
// Global configuration parameters
// -------------------------------
// numruns: how many independent runs for each slope configuration
// FirstSlope..LastSlope: iterate fixed overhead slope s in steps of 2
#define numruns 1     // Number of runs for each configuration
#define FirstSlope 16 // The program starts with this slope coeff
#define LastSlope 18  // The program ends with this slope coeff

// Core model sizes (match the paper’s large-scale experiments by setting
// n=10, bsize=24, K=200; each run is T weeks long)
#define n 10    // Number of goods
#define bsize 24 // Number of each type of person
#define K 200    // Number of store locations
// #define n 10     // Number of goods
// #define bsize 24 // Number of each type of person
// #define K 200    // Number of store locations

// Behavioral/tuning parameters
#define xMax 200    // Maximal target (animal spirits) for a new entrant
#define lambda 0.05 // Autonomous search intensity (higher: more search)
#define alpha 0.25  // Target adjustment speed (adaptive expectations)
#define theta 0.01  // Exit probability if a shop is unprofitable
#define C 5.0       // Setup/normal-return cost per good handled
#define f1 0.0      // Fixed overhead constant term
#define persist 10  // # of years criteria must persist for flags

// Horizon and randomness
// In the paper T≈20,000 weeks; here reduced for speed.
// Each “year” is 50 weeks; reporting uses those blocks.
#define T 20000
#define RAND 1    // Initial RNG seed base (per configuration)
#define RptPer 1  // Reporting period in years (multiplies 50-week blocks)
#define prtoscr 1 // 0 = quiet, 1 = print status to screen

// Derived constants/macros
#define m (bsize * n * (n - 1))     // Number of traders: bsize of each (i≠j) type
#define f(i) (f1 + (i - 1) * slope) // Overhead cost schedule for good i (depends on slope)
// Pricing rule F(tr0,tr1,f_other): posted buying price of good 0
//   Implements full-cost/target-pricing logic from the paper.
#define G(x, y, z) ((y - z - C > 0) ? ((y - C - z) / x) : 0) // F(tr0,tr1,f1) = buying price of good 0
#define MAX(a, b) ((a) > (b) ? (a) : (b))

double F(double tr0, double tr1, double f_other)
{
  if (DEBUG) {
    printf("priceF %.2f %.2f %.2f\n", tr0, tr1, f_other);
    printf("results -> %.2f\n", G(tr0, tr1, f_other));
  }
  return G(tr0, tr1, f_other);
}

#define PRINT_LOOP_N 6 // Number of loops between progress reports
// -----------------------------------
// Global state (agents, shops, stats)
// -----------------------------------
// Traders: r=1..m. s[r]=produced good, d[r]=desired good. q[r]=which side is supply in barter.
int s[m + 1], d[m + 1]; // s[r] d[r] are r's supply and demand goods
int sh[2][m + 1];       // Trader r sells to shop[0], buys from [1]
int q[m + 1];           // Index of the trader's supply good in barter

// Shops: k=1..K
int NS;            // Number of active stores
int BS;            // Number of active stores *not* using the money good
int active[K + 1]; // 0 if inactive shop location
int g[2][K + 1];   // Goods traded at shop k, ordered so g[0]<=g[1]

// Agent lists by production/consumption good (for sampling “friends”)
int produces[n + 1][(n - 1) * bsize + 1]; // Names of those endowed with good i
int consumes[n + 1][(n - 1) * bsize + 1]; // Names of those who like good i
int numprod[n + 1];                       // Count of producers of i
int numcons[n + 1];                       // Count of consumers of i

// Ownership/bookkeeping
int owner[K + 1];      // The individual that owns shop k
int familyshop[m + 1]; // The (unique) shop owned by r (0 if none)

// Timers/loop indices/temporaries
int devyear; // First year of “full development” (≥99% participants)
int monyear; // First year of “monetary exchange” (≥99% use same money)
int h, i, j;
int k, ss;                        // Used to index firms
int r, fr;                        // Used to index individuals and friends
int extra;                        // # of candidate shops sampled in search
int t;                            // Week index
int run;                          // Run index
int bestbarter;                   // Best direct (barter) shop found in search
int endcount, devcount;           // Persistence counters for stopping criteria
int a, b;                         // Temporary shop indices in search/utility
int m0, m1, ma, mb, Nshop;        // Side flags: which side matches s[r]/d[r], count shops
int c[6];                         // Sampled shop list (up to 5 extra + 1 base)
int line[m + 1];                  // Random weekly lineup of agents
int monetary, fulldev, moneygood; // Indicators and the dominant money commodity

// Stationary (SME) comparators and summary stats
double W;           // Wholesale price (SME) denominator
double Pinv[n + 1]; // Inverse retail prices in SME by good
double SurpSME;     // Surplus in SME with the same money good

// Pricing/targets/incomes
double slope;                                // Incremental fixed cost (varied across runs)
double targ[2];                              // Entry “animal spirits” targets (trials)
double tr[2][K + 1];                         // Shop income targets for the two goods
double Fmon;                                 // Theoretical max # money traders (for 99% test)
double U, Ucomp, Ubarter;                    // Utility indices: current, candidate, best barter
double part;                                 // # of participating traders this week
double P[2][K + 1];                          // Posted buying prices for the two goods
double y[2][K + 1];                          // Realized shop incomes this week (per good)
double moneytraders, usingmax, Csurp, Psurp; // Counts & welfare: agents’ and producers’ surplus
double usingmoney[n + 1];                    // Distribution of which good is used as money
time_t firstbegin, finish;
clock_t clock_begin, clock_finish;
FILE *stream;
double R[2];                         // Relative RMSEs for wholesale and inverse retail
double vol[2][n + 1], avp[2][n + 1]; // Volume-weighted average prices vs SME baselines
                                     // avp[0][j]: wholesale price relative to SME
                                     // avp[1][j]: inverse retail price relative to SME

// ------------------------
// Function declarations
// ------------------------
void Init(void);
int groupsize(int, int);
void InitRun(void);
void lineup(void);
int Research(void);
int comrade(void);
int soulmate(void);
double Usample(void);
double Utility(void);
void addshop(int);
void tryBarter(void);
void dropshop(void);
void tryOne(void);
void tryTwo(void);
int Calc1(void);
void Calc2(void);
void RMSE(void);
void report(int);
void print_debug(char *);
void print_shop(int);
void print_trader(int);

// ------------------------
// RNG: PCG32 implementation
// ------------------------
typedef struct
{
  uint64_t state;
  uint64_t inc; // must be odd
} pcg32_t;
pcg32_t rng;
uint32_t pcg32_next(pcg32_t *rng);
void pcg32_seed(pcg32_t *rng, uint64_t initstate, uint64_t initseq);
#define SEED_RANDOM(seed) pcg32_seed(&rng, seed, seed)
#define RANDOM_N(x) trap_rnd_gen(x)
#define rnd01 ((1 + RANDOM_N(1000)) / 1000.0) // Uniform(0,1] with 0.001 granularity
void report_trader(int r);
int trap_rnd_gen(int i);

int main()
{
  Init();
  for (int _slope = FirstSlope; _slope <= LastSlope; _slope += 2)
  { // Sweep overhead slope
    slope = _slope;
    for (run = 1; run <= numruns; run++)
    { // Repeat runs per slope
      InitRun();
      for (t = 1; t <= T; t++)
      { // Weekly loop

        // Entry: occasional entrepreneur tries to open a shop
        do
        {
          r = RANDOM_N(m) + 1; // Random prospective owner
          if (NS < K && familyshop[r] == 0)
          {                 // Must have an empty slot & not already own
            j = Research(); // Market research around r’s type
            if (j > 0)
            { // If expected to attract business both sides
              NS++;
              k = 1;
              while (active[k] == 1)
                k++; // Pick first free shop location
              active[k] = 1;
              g[1 - q[r]][k] = d[r];
              g[q[r]][k] = s[r]; // Trade pair ordered by q[r]
              for (h = 0; h < 2; h++)
              {
                tr[h][k] = targ[h];                                // Initialize targets
                P[h][k] = F(targ[h], targ[1 - h], f(g[1 - h][k])); // Post initial buying prices
              }
              sh[0][r] = k;
              sh[1][r] = 0;
              familyshop[r] = k;
              owner[k] = r; // Owner links to own shop
            }
          }
        } while (NS == 0); // Ensure at least one shop exists

        if (DEBUG)  {
          print_debug("Weekly entry completed");

          // Matching (search/adoption): agents sample a small set of shops
          printf("line: ");
          for (i = 1; i <= m; i++) {
            printf("%d ", line[i]);
          }
          printf("\n");          
        }
      
        for (i = 1; i <= m; i++)
        {
          r = line[i];   // Visit in randomized lineup order
          U = Utility(); // Current utility (0 if unmatched)
          double rr =  rnd01;
          if (rr < (U > 0 ? lambda : 1) && familyshop[r] == 0)
          {
            // Search: if not currently matched (U==0) search with prob 1, else with prob lambda
            extra = 0;
            for (h = 0; h <= 1; h++)
              c[h] = sh[h][r]; // Seed candidate list with current links
            fr = comrade();    // Friend with same production good → use their outlet
            addshop(sh[0][fr]);
            fr = soulmate(); // Friend with same consumption good → use their source
            addshop(sh[1][fr]);
            addshop(RANDOM_N(K) + 1); // One random shop (stranger)
            // NOTE: extra is set within addshop!!!
            if (extra)
            { // If any candidates were found, compare options
              Ucomp = U;
              if (extra)
                tryBarter(); // Check direct (barter) option
              if (extra && (g[1 - q[r]][c[0]] != d[r] || P[q[r]][c[0]] == 0))
                tryOne(); // One-side match
              if (extra > 1)
                tryTwo(); // Two-shop indirect (money) chain
              if (Ucomp < Ubarter)
              {
                sh[0][r] = bestbarter;
                sh[1][r] = 0;
              } // Adopt best
              else
              {
                sh[0][r] = c[0];
                sh[1][r] = c[1];
              }
            }
          }
          if (DEBUG)
            print_trader(r);
        }

        if (DEBUG)  {
          print_debug("Weekly trade begins");
        }

        // Trade (accounting): tally shop incomes from adopted relationships
        for (k = 1; k <= K; k++)
          if (active[k] > 0)
            for (h = 0; h <= 1; h++)
              y[h][k] = 0;
        for (r = 1; r <= m; r++)
          if (sh[0][r] > 0)
          {
            a = sh[0][r];
            b = sh[1][r];
            ma = (g[0][a] == s[r]);
            mb = (g[0][b] == d[r]);
            if (g[ma][a] == d[r])
              y[1 - ma][a]++; // Direct barter at outlet
            else if (g[ma][a] == g[mb][b])
            { // Indirect via common intermediary
              y[1 - ma][a]++;
              y[mb][b] += P[1 - ma][a]; // Transfer priced by outlet’s posted price
            }
          }
        // Exit: unprofitable shops close with probability theta
        for (k = 1; k <= K; k++)
          if (active[k])
            if ((y[0][k] - P[1][k] * y[1][k] - f(g[0][k]) <= 0 ||
                 y[1][k] - P[0][k] * y[0][k] - f(g[1][k]) <= 0) &&
                rnd01 <= theta)
            {
              active[k] = 0;
              NS--;
              familyshop[owner[k]] = 0;
              owner[k] = 0;
              for (h = 0; h <= 1; h++)
              {
                y[h][k] = 0;
                tr[h][k] = 0;
                g[h][k] = 0;
                P[h][k] = 0;
              }
              for (r = 1; r <= m; r++)
                for (h = 0; h <= 1; h++)
                  if (sh[h][r] == k)
                    sh[h][r] = 0; // Sever links
            }

        if (DEBUG)  {
          print_debug("Weekly trade and exit completed");
        }

        for (k = 1; k <= K; k++)
          if (active[k] && DEBUG)
            printf("Shop %d: tr0=%.2f tr1=%.2f\n", k, tr[0][k], tr[1][k]);

        // Update targets (adaptive) and recompute posted prices
        for (k = 1; k <= K; k++)
          if (active[k])
          {
            for (h = 0; h <= 1; h++)
            {
              tr[h][k] += alpha * (y[h][k] - tr[h][k]); // Partial adjustment toward realized y
            }
            P[0][k] = F(tr[0][k], tr[1][k], f(g[1][k])); // New posted prices from targets
            P[1][k] = F(tr[1][k], tr[0][k], f(g[0][k]));
          }

        // Periodic check during run (here printed every week with prtoscr)
        if(t%(PRINT_LOOP_N*RptPer)==0) {
          report(t);
          monetary = Calc1(); // Update participation & money stats          // if(t%(50*RptPer*20)==0) getch();
          if (monetary == 1)
            break;
        }
        // End of the week loop
      }
      Calc2(); // Final measurements for the run
      clock_finish = clock();
      printf("Run number %d. Time elapsed: %.2f seconds.\n", run, (clock_finish - clock_begin) / (double)CLOCKS_PER_SEC);
      printf("Slope equals %-.0f, xMax equals %d\n\n", slope, xMax);
      stream = fopen("evol.fil", "a");
      fprintf(stream, "%5d %5.0f %3d %4.0f %3d %4.0f %4.0f %3d %4.0f %6.0f %5.0f %3d %3d %3d %6.3f %6.3f %4d %4d\n",
              run, slope, fulldev, part, monetary, moneytraders, usingmax,
              moneygood, Csurp, Psurp, SurpSME, Nshop, NS, BS, R[0], R[1], devyear, monyear);
      fclose(stream);

      // End of the run loop
    }
    stream = fopen("evol.fil", "a");
    fprintf(stream, "\n");
    fclose(stream);
    // End of the fmean loop (slope sweep)
  }
  stream = fopen("evol.fil", "a");
  fprintf(stream, "\n\n  n bsiz   K  xMax  lambda alpha  theta  C  Nrun Pst(yr)  T RND\n");
  fprintf(stream, "%3d %4d %3d %5d %6.3f %6.3f %6.3f %2.0f %5d %3d %6d %3d\n",
          n, bsize, K, xMax, lambda, alpha, theta, C, numruns, persist, T, RAND);
  time(&finish);
  fprintf(stream, "\n    Started : %s", ctime(&firstbegin));
  fprintf(stream, "    Finished: %s\n", ctime(&finish));
  fclose(stream);
  return (0);
}

// ----------------------------------------
// Init(): one-time construction of economy
// ----------------------------------------
void Init(void)
{
  time(&firstbegin);
  stream = fopen("evol.fil", "w+");
  fprintf(stream, "  run slope dev part mon mtrd usmx Mgd Csur   Psur  Esur Nsh  NS  BS  RMS_0  RMS_1  dyr  myr\n");
  fclose(stream);
  r = 0;
  d[0] = 0;
  s[0] = 0;
  Pinv[0] = 0.0;
  // Clear producer/consumer rosters
  for (i = 0; i <= n; i++)
  {
    numprod[i] = 0;
    numcons[i] = 0;
    for (j = 0; j <= (n - 1) * bsize; j++)
    {
      produces[i][j] = 0;
      consumes[i][j] = 0;
    }
  }
  // Create all agent types (i!=j) with multiplicity bsize
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      for (k = 1; k <= groupsize(i, j); k++)
      {
        r++;
        s[r] = i;
        d[r] = j;
        q[r] = (s[r] > d[r]); // q[r]=1 if i>j else 0 (orders g pair)
        numprod[i]++;
        produces[i][numprod[i]] = r;
        numcons[j]++;
        consumes[j][numcons[j]] = r;
      }
  // Theoretical cap on # of “money traders”: bsize*(n-2)*(n-1)
  Fmon = bsize * (n - 2) * (n - 1);
  for (h = 0; h <= 1; h++)
  {
    vol[h][0] = 0;
    avp[h][0] = 0;
  }
}

// Size of each (i,j) group: bsize if i!=j, else 0
int groupsize(int v, int w)
{
  if (v == w)
    return (0);
  else
    return (bsize);
}

// ----------------------------------------
// InitRun(): reset state for a fresh run
// ----------------------------------------
void InitRun(void)
{
  printf("Number  Using  Using  Using  Using  Using  Using                \n");
  printf("Active  Money  good1  good2  good3  good4  good5   Year   NS   \n");
  endcount = 0;
  fulldev = 0;
  devyear = -1;
  monyear = -1;
  devcount = 0;
  SEED_RANDOM(RAND + run - 1);
  clock_begin = clock();
  t = 0;
  NS = 0;
  for (r = 0; r <= m; r++)
  {
    for (h = 0; h <= 1; h++)
      sh[h][r] = 0;
    familyshop[r] = 0;
  }
  for (k = 0; k <= K; k++)
  {
    active[k] = 0;
    owner[k] = 0;
    for (h = 0; h <= 1; h++)
    {
      P[h][k] = 0;
      y[h][k] = 0;
      tr[h][k] = 0;
      g[h][k] = 0;
    }
  }
  lineup(); // Randomize agent visitation order
}

// ----------------------------------------
// lineup(): generate random permutation
// ----------------------------------------
void lineup(void)
{
  int i, j, k;
  int start[m + 1];
  line[0] = 0;
  start[0] = 0;
  for (i = 1; i <= m; i++)
    start[i] = i;
  for (j = 1; j <= m; j++)
  {
    k = RANDOM_N(m - j + 1) + 1;
    line[j] = start[k];
    for (i = k; i <= m - j; i++)
      start[i] = start[i + 1];
  }
}

// -----------------------------------------------------------
// Research(): entry “market research” at prospective owner r.
// Returns 1 if opening the (s[r],d[r]) shop passes tests:
// - Can attract on both sides via friend/stranger samples.
// Else returns 0.
// -----------------------------------------------------------
int Research(void)
{ // Assume that we have tested U(r) first. Don't need to use m1
  for (h = 0; h <= 1; h++)
    targ[h] = RANDOM_N(xMax) + 1; // or m0 after that until we get to search.
  double U;
  double P0 = F(targ[0], targ[1], f(d[r])); // price of s[r]
  double P1 = F(targ[1], targ[0], f(s[r])); // price of d[r]
  U = 0;                                    // test a comrade and a friend
  fr = comrade();
  Ucomp = Usample();
  if (d[fr] == d[r])
    U = P0;
  else if (g[m1][sh[1][fr]] == d[r])
    U = P0 * P[m1][sh[1][fr]];
  if (U < Ucomp)
  {
    U = 0;
    fr = soulmate();
    Ucomp = Usample();
    if (s[fr] == s[r])
      U = P0;
    else if (g[m0][sh[0][fr]] == s[r])
      U = P[1 - m0][sh[0][fr]] * P0;
    if (U < Ucomp)
      return (0);
  }
  U = 0; // test a stranger who likes s[r] and one who produces d[r]
  int k = 1 + RANDOM_N(numcons[s[r]] - 1);
  fr = consumes[s[r]][k];
  Ucomp = Usample();
  if (s[fr] == d[r])
    U = P1;
  else if (g[m0][sh[0][fr]] == d[r])
    U = P[1 - m0][sh[0][fr]] * P1;
  if (U < Ucomp)
  {
    U = 0;
    k = 1 + RANDOM_N(numprod[d[r]] - 1);
    fr = produces[d[r]][k];
    Ucomp = Usample();
    if (d[fr] == s[r])
      U = P1;
    else if (g[m1][sh[1][fr]] == s[r])
      U = P1 * P[m1][sh[1][fr]];
    if (U < Ucomp)
      return (0);
  }
  return (1);
}

// Pick a “comrade”: someone with same production good s[r]
// Ensures not to pick r themself.
int comrade(void)
{
  int k = 1 + RANDOM_N(numprod[s[r]] - 1);
  int fr = produces[s[r]][k];
  if (fr >= r)
    fr = produces[s[r]][k + 1];
  return (fr);
}

// Pick a “soulmate”: someone with same desired good d[r]
int soulmate(void)
{
  int k = 1 + RANDOM_N(numcons[d[r]] - 1);
  int fr = consumes[d[r]][k];
  if (fr >= r)
    fr = consumes[d[r]][k + 1];
  return (fr);
}

// Usample(): compute fr’s attainable consumption via current links
double Usample(void)
{
  m0 = (g[0][sh[0][fr]] == s[fr]);
  m1 = (g[0][sh[1][fr]] == d[fr]);
  if (DEBUG) {
    printf("[1] SET m1 to %d\n", m1);
  }
  double X = 0;
  if (sh[0][fr] > 0)
  {
    if (g[m0][sh[0][fr]] == d[fr])
      X = P[1 - m0][sh[0][fr]]; // Direct
    else if (g[m0][sh[0][fr]] == g[m1][sh[1][fr]])
      X = P[1 - m0][sh[0][fr]] * P[m1][sh[1][fr]]; // Indirect via common good
  }
  return (X);
}

// Utility(): r’s present attainable consumption level via (sh[0],sh[1])
double Utility(void)
{
  m0 = (g[0][sh[0][r]] == s[r]);
  m1 = (g[0][sh[1][r]] == d[r]);
  if (DEBUG) {
    printf("[2] SET m1 to %d\n", m1);
  }
  double X = 0;
  if (sh[0][r] > 0)
  {
    if (g[m0][sh[0][r]] == d[r])
      X = P[1 - m0][sh[0][r]];
    else if (g[m0][sh[0][r]] == g[m1][sh[1][r]])
      X = P[1 - m0][sh[0][r]] * P[m1][sh[1][r]];
  }
  return (X);
}

// addshop(): add candidate shop ss to sample list c[*] if active and relevant to r
void addshop(int ss)
{
  if (active[ss])
    if (g[0][ss] == s[r] || g[1][ss] == s[r] || g[0][ss] == d[r] || g[1][ss] == d[r])
    {
      int b = 0;
      for (a = 0; a <= 1 + extra; a++)
        b += (c[a] == ss);
      if (b == 0)
      {
        extra++;
        c[1 + extra] = ss;
      }
    }
}

// tryBarter(): among candidates, pick best direct (barter) shop for r
void tryBarter(void)
{
  bestbarter = 0;
  Ubarter = 0;
  for (a = 2; a <= 1 + extra; a++)
    if (g[q[r]][c[a]] == s[r] && g[1 - q[r]][c[a]] == d[r])
    {
      if (P[q[r]][c[a]] > MAX(Ucomp, Ubarter))
      {
        bestbarter = c[a];
        Ubarter = P[q[r]][c[a]];
        dropshop();
      }
    }
}

// dropshop(): remove c[a] from the candidate list
void dropshop(void)
{
  for (int b = a; b < 1 + extra; b++)
    c[b] = c[b + 1];
  a--;
  extra--;
}

// tryOne(): try to improve either outlet or source when only one side matches
void tryOne(void)
{
  for (a = 2; a <= 1 + extra; a++)
    if ((g[0][c[a]] == s[r]) || (g[1][c[a]] == s[r]))
    {
      ma = (g[0][c[a]] == s[r]);
      if ((g[ma][c[a]] == g[m1][c[1]]) || (P[1 - m0][c[0]] == 0))
      {
        if (P[1 - m0][c[0]] < P[1 - ma][c[a]])
        {
          c[0] = c[a];
          m0 = ma;
          Ucomp = (g[ma][c[a]] == g[m1][c[1]] ? P[1 - ma][c[a]] * P[m1][c[1]] : 0);
          dropshop();
        }
      }
    }
    else
    {
      ma = (g[0][c[a]] == d[r]);
      if ((g[ma][c[a]] == g[m0][c[0]]) || (P[m1][c[1]] == 0))
      {
        if (P[m1][c[1]] < P[ma][c[a]])
        {
          c[1] = c[a];
          m1 = ma;
          if (DEBUG) {
            printf("[3] SET m1 to %d\n", m1);
          }

          Ucomp = (g[ma][c[a]] == g[m0][c[0]] ? P[1 - m0][c[0]] * P[ma][c[a]] : 0);
          dropshop();
        }
      }
    }
}

// tryTwo(): consider all two-shop chains (outlet, source) that form a valid intermediary
void tryTwo(void)
{
  for (a = 2; a <= 1 + extra; a++)
    if (g[0][c[a]] == s[r] || g[1][c[a]] == s[r])
    {
      for (b = 2; b <= 1 + extra; b++)
        if (a != b && (g[0][c[b]] == d[r] || g[1][c[b]] == d[r]))
        {
          ma = (g[0][c[a]] == s[r]);
          mb = (g[0][c[b]] == d[r]);
          if (g[ma][c[a]] == g[mb][c[b]] && Ucomp < P[1 - ma][c[a]] * P[mb][c[b]])
          {
            Ucomp = P[1 - ma][c[a]] * P[mb][c[b]];
            c[0] = c[a];
            c[1] = c[b];
            m0 = ma;
            m1 = mb;
            printf("[4] SET m1 to %d\n", m1);
          }
        }
    }
}

// --------------------------------------------------------
// Calc1(): compute participation & money dominance metrics
//   part: # agents who can trade (direct or indirect)
//   moneytraders: # using an intermediary (two-shop path)
//   usingmoney[i]: count of agents whose intermediary is good i
// Also updates persistence counters to mark full development
// and monetary exchange years (devyear, monyear).
// --------------------------------------------------------
int Calc1(void)
{
  part = 0;
  moneytraders = 0;
  for (j = 1; j <= n; j++)
    usingmoney[j] = 0;
  for (r = 1; r <= m; r++)
    if (sh[0][r] > 0)
    {
      a = sh[0][r];
      b = sh[1][r];
      ma = (g[0][a] == s[r]);
      mb = (g[0][b] == d[r]);
      part += ((g[ma][a] == d[r] || g[ma][a] == g[mb][b]));
      if (g[ma][a] == g[mb][b])
      {
        moneytraders++;
        usingmoney[g[mb][b]]++;
      }
    }
  if (fulldev == 0)
  {
    if (part >= .99 * m)
    {
      if (devcount == 0)
        devyear = t / 50;
      devcount++;
    }
    else
      devcount = 0;
    if (devcount >= persist)
      fulldev = 1;
  }
  usingmax = 0;
  moneygood = 0;
  for (i = 1; i <= n; i++)
    if (usingmoney[i] > usingmax)
    {
      usingmax = usingmoney[i];
      moneygood = i;
    }
  if (usingmoney[moneygood] >= .99 * Fmon)
  {
    if (endcount == 0)
      monyear = t / 50;
    endcount++;
  }
  else
    endcount = 0;
  if (endcount < persist)
    return (0);
  else
    return (1);
}

// --------------------------------------------------------
// Calc2(): end-of-run summary — flags, shop counts, SME
// comparators, welfare (consumer surplus Csurp, producer
// surplus Psurp), and diagnostic price RMSEs.
// --------------------------------------------------------
void Calc2(void)
{
  if (monetary == 0)
    monyear = -1;
  if (fulldev == 0)
    devyear = -1;
  BS = 0;
  for (k = 1; k <= K; k++)
    if (active[k] && g[0][k] != moneygood && g[1][k] != moneygood)
      BS++;
  W = 1.0 - (f(moneygood) + C) / bsize;
  if (W > 0.0)
    for (i = 1; i <= n; i++)
      Pinv[i] = ((m / n) - (f(i) + C)) / ((m / n) - (n - 2) * (f(moneygood) + C));
  SurpSME = m - n * f1 - (slope / 2) * n * (n - 1) - (n - 2) * f(moneygood);
  RMSE();
  Csurp = 0;
  for (r = 1; r <= m; r++)
    Csurp += Utility(); // Aggregate consumer surplus proxy
  Psurp = 0;
  Nshop = 0;
  for (k = 1; k <= K; k++)
    if (active[k])
    {
      for (h = 0; h <= 1; h++)
        Psurp += (y[h][k] - f(g[h][k]) - P[1 - h][k] * y[1 - h][k]); // Producer surplus proxy
      if (y[q[owner[k]]][k] > 1 || y[1 - q[owner[k]]][k] > 0)
        Nshop++; // Nonempty shops (activity)
    }
}

// --------------------------------------------------------
// RMSE(): compute price RMSEs relative to SME benchmarks.
// Uses volume-weighted averages (wholesale and inverse retail)
// by non-money goods paired with the chosen moneygood.
// --------------------------------------------------------
void RMSE(void)
{
  if (W > 0.0)
  {
    for (h = 0; h <= 1; h++)
      for (i = 1; i <= n; i++)
      {
        vol[h][i] = 0;
        avp[h][i] = 0;
      }
    for (k = 1; k <= K; k++)
      if (g[0][k] == moneygood || g[1][k] == moneygood)
      {
        ma = (g[1][k] == moneygood);
        i = g[1 - ma][k];
        vol[0][i] += y[1 - ma][k];
        if (vol[0][i])
          avp[0][i] += (y[1 - ma][k] / vol[0][i]) * (P[1 - ma][k] - avp[0][i]);
        vol[1][i] += y[ma][k];
        if (vol[1][i])
          avp[1][i] += (y[ma][k] / vol[1][i]) * (P[ma][k] - avp[1][i]);
      }
    for (i = 1; i <= n; i++)
    {
      avp[0][i] /= W;
      avp[1][i] /= Pinv[i];
    }
    for (h = 0; h <= 1; h++)
    {
      R[h] = 0;
      for (i = 1; i <= n; i++)
        if (i != moneygood)
          R[h] += (1 - avp[h][i]) * (1 - avp[h][i]);
      R[h] /= (n - 1);
      R[h] = sqrt(R[h]);
    }
  }
  else
  {
    for (h = 0; h <= 1; h++)
      R[h] = -1.0;
  }
}

// -----------------
// PCG32 RNG methods
// -----------------
uint32_t pcg32_next(pcg32_t *rng)
{
  uint64_t oldstate = rng->state;
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1ULL);
  uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  uint32_t rot = (uint32_t)(oldstate >> 59u);
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

void pcg32_seed(pcg32_t *rng, uint64_t initstate, uint64_t initseq)
{
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  pcg32_next(rng);
  rng->state += initstate;
  pcg32_next(rng);
}

void report(int t) 
{
  if (prtoscr != 0)
  {
    if (t == -1)
      printf("***");
    printf("%6.0f %6.0f ", part, moneytraders);
    for (b = 1; b <= 5; b++)
      printf("%6.0f ", usingmoney[b]);
    printf("%6d %4d\n", t / 50, NS);
  }
}

void print_debug(char *title) {
  printf("--------------------------------------------------\n");
  printf(title);
  printf("\nShops:\n");
  for (int i = 1; i <= K; i++) {
    print_shop(i);
  }
  printf("Traders:\n");
  for (int i = 1; i <= m; i++) {
    print_trader(i);
  }
  printf("--------------------------------------------------\n");
}

void print_shop(int i) {
  printf("Shop %d: Shop{active=%d, g=[%d,%d], owner=%d, P=[%f,%f], y=[%f,%f], tr=[%f,%f]}\n",
         i, active[i], g[0][i], g[1][i], owner[i], P[0][i], P[1][i], y[0][i], y[1][i], tr[0][i], tr[1][i]);
}

void print_trader(int i) 
{
  printf("Trader{s=%d, d=%d, q=%d, sh=[%d,%d], familyshop=%d}\n",
         s[i], d[i], q[i], sh[0][i], sh[1][i], familyshop[i]);
}

int trap_rnd_gen(int i)
{
  uint32_t u = pcg32_next(&rng);
  int next_n = (int)(u % i);
  if (DEBUG) {
    printf("RND: %u\n", u);
  }
  return next_n;
}
