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

#define numruns 2                               // Number of runs for each configuration
#define FirstSlope 16                           // The program starts with this slope coeff
#define LastSlope 18                            // The program ends with this slope coeff

#define n 10           				            // Number of goods
#define bsize 24                                // Number of each type of person
#define K 200                                   // Number of store locations
#define xMax 200                                // Maximal target of a new entrant
#define lambda 0.05   		    	            // Autonomous search intensity
#define alpha 0.25				                // Adjustment speed
#define theta 0.01      			            // Failure rate of unprofitable stores
#define C 5.0				                    // Setup cost in flow terms
#define f1 0.0                                  // The constant term in the overhead cost shdl.
#define persist 10                              // Number of years before "monetary"
#define T (20000 +50*(persist-1))               // Maximal number of weeks in each run
#define RAND 1                                  // Initial seed for each new configuration
#define RptPer 1                                // Reporting Period (in years)
#define prtoscr 1                               // Set this to 0 to kill printing to screen
#define m (bsize*n*(n-1))			            // Number of people - bsize of each type
#define f(i) (f1+(i-1)*slope)                   // Overhead cost for dealing in i
#define F(x,y,z) ((y-z-C>0)?((y-C-z)/x):0)      // F(tr0,tr1,f1) = buying price of good 0
#define MAX(a,b) ((a) > (b) ? (a) : (b))

int s[m+1], d[m+1];					            // s[r] d[r] are r's supply and demand goods
int sh[2][m+1];						            // Trader r sells to shop[0], buys from [1]
int q[m+1];                                     // Index of the trader's supply good in barter
int NS;									        // Number of stores
int BS;                                         // Number of stores not using moneygood
int active[K+1];						        // 0 if inactive shop location
int g[2][K+1];  						        // Goods traded at shop k, with smaller first.
int produces[n+1][(n-1)*bsize+1];               // Names of those endowed with good i
int consumes[n+1][(n-1)*bsize+1];               // Names of those who like good i
int numprod[n+1];                               // Number of those endowed with i
int numcons[n+1];                               // Number of those liking i
int owner[K+1];                                 // The individual that owns shop k
int familyshop[m+1];                            // The shop owned by r
int devyear;                                    // The first year of full development
int monyear;                                    // The first year of monetary exchange
int h, i, j;
int k, ss;                                      // Used to index firms
int r, fr;                                      // Used to index individuals and friends
int extra;                                      // Used to index the sample size of searcher
int t;                                          // Used for week
int run;                                        // Used to index the run
int bestbarter;                                 // Used to track barter shop in search
int endcount, devcount;                         // Used in the stopping rule for a run
int a, b;                                       // Used to index firms in search and Utility
int m0, m1, ma, mb, Nshop;                      // Used to index the money good of a firm
int c[6];                                       // The list of firms being sampled
int line[m+1];                                  // Random lineup before each week's search
int monetary, fulldev, moneygood;               // Summary statistics for each run
double W;                                       // Wholesale price in SME.
double Pinv[n+1];                               // Inverse retail prices in SME.
double SurpSME;                                 // Surplus in SME with same moneygood
double slope;      								// Incremental fixed cost
double targ[2];                                 // The target established in Research()
double tr[2][K+1];						        // Target values for y
double Fmon;                                    // Number of money traders in monetary SE
double U, Ucomp, Ubarter;			            // utility indices
double part;								    // number of participants
double P[2][K+1];						        // P[h][k] = buying price of good h
double y[2][K+1];						        // y[h][k] = shop income of good h
double moneytraders, usingmax, Csurp, Psurp;
double usingmoney[n+1];
time_t firstbegin, finish;
clock_t clock_begin, clock_finish;
FILE *stream;
double R[2];                                    // The relative RMSE on whole and inv prices.
double vol[2][n+1],avp[2][n+1];                 // avp[0][j]is the weighted average wholesale
                                                // price of [j] relative to the SME,
                                                // and avp[1][j] is the weighted avg
                                                // inverse retail price rel to SME
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

typedef struct {
  uint64_t state;
  uint64_t inc;   // must be odd
} pcg32_t;
pcg32_t rng;
uint32_t pcg32_next(pcg32_t *rng);
void pcg32_seed(pcg32_t *rng, uint64_t initstate, uint64_t initseq);
#define SEED_RANDOM(seed) pcg32_seed(&rng, seed, seed)
#define RANDOM_N(x) (pcg32_next(&rng) % (x))
#define rnd01 ((1+RANDOM_N(1000))/1000.0)

int main( ) {
  Init();
  for(slope=FirstSlope;slope<=LastSlope;slope+=2) {
    for(run=1;run<=numruns;run++)  {
      InitRun();
      for(t=1;t<=T;t++) {

                  // Entry
        do {
          r=RANDOM_N(m)+1;
          if( NS<K && familyshop[r]==0) {
            j=Research();
            if(j>0) {
              NS++; k=1; while(active[k]==1) k++;
              active[k]=1; g[1-q[r]][k]=d[r]; g[q[r]][k]=s[r];
              for(h=0;h<2;h++) {
                tr[h][k]=targ[h];
                P[h][k] = F(targ[h],targ[1-h],f(g[1-h][k]));
              }
              sh[0][r]=k; sh[1][r]=0; familyshop[r]=k; owner[k]=r;
            }
          }
        }
        while (NS==0);

                  // Matching
        for(i=1;i<=m;i++)  {
          r=line[i];
          U=Utility();
          if(rnd01<( U>0 ? lambda : 1 ) && familyshop[r]==0) {
            extra=0; for(h=0;h<=1;h++) c[h]=sh[h][r];
            fr=comrade();         // Check 2 friends' shops and find 1 on your own
            addshop(sh[0][fr]);
            fr=soulmate();
            addshop(sh[1][fr]);
            addshop(RANDOM_N(K)+1);
            if (extra) {          // Compare utilities
              Ucomp=U;
              if(extra) tryBarter();
              if(extra && ( g[1-q[r]][c[0]]!=d[r] || P[q[r]][c[0]]==0  ) ) tryOne();
              if(extra>1) tryTwo();
              if(Ucomp<Ubarter) {sh[0][r]=bestbarter; sh[1][r]=0;}
              else {sh[0][r]=c[0]; sh[1][r]=c[1];}
            }
          }
        }

                  // Trade
        for(k=1;k<=K;k++) if(active[k]>0) for(h=0;h<=1;h++) y[h][k]=0;
        for(r=1;r<=m;r++) if(sh[0][r]>0) {
          a=sh[0][r]; b=sh[1][r];ma=(g[0][a]==s[r]);mb=(g[0][b]==d[r]);
          if(g[ma][a]==d[r]) y[1-ma][a]++;
          else if(g[ma][a]==g[mb][b]) {
              y[1-ma][a]++;
              y[mb][b]+=P[1-ma][a];
          }
        }
                  // Exit
        for(k=1;k<=K;k++) if(active[k])
        if( ( y[0][k]-P[1][k]*y[1][k]-f(g[0][k]) <= 0 ||
              y[1][k]-P[0][k]*y[0][k]-f(g[1][k]) <= 0 ) && rnd01<=theta ){
          active[k]=0; NS--; familyshop[owner[k]]=0; owner[k]=0;
          for(h=0;h<=1;h++)	{
            y[h][k]=0; tr[h][k]=0; g[h][k]=0; P[h][k]=0;
          }
          for(r=1;r<=m;r++) for(h=0;h<=1;h++) if(sh[h][r]==k) sh[h][r]=0;
        }

                  // Update targets and prices
        for(k=1;k<=K;k++) if(active[k]) {
          for(h=0;h<=1;h++) {
            tr[h][k]+=alpha*(y[h][k]-tr[h][k]);
          }
          P[0][k]=F(tr[0][k],tr[1][k],f(g[1][k]));
          P[1][k]=F(tr[1][k],tr[0][k],f(g[0][k]));
        }


                  // Periodic check during run
        if(t%(50*RptPer)==0) {
          monetary=Calc1();
          i=prtoscr;
          if(i!=0) {
            printf("%6.0f %6.0f ", part, moneytraders );
            for(b=1;b<=5;b++) printf("%6.0f ", usingmoney[b]);
            printf("%6d %4d\n", t/50, NS);
        }
//if(t%(50*RptPer*20)==0) getch();
          if(monetary==1) break;
        }
              // End of the week loop
      }
      Calc2();
      clock_finish = clock();
      printf("Run number %d. Time elapsed: %.2f seconds.\n", run, (clock_finish - clock_begin)/(double)CLOCKS_PER_SEC);
      printf("Slope equals %-.0f, xMax equals %d\n\n", slope, xMax);
      stream=fopen("evol.fil", "a");
      fprintf(stream, "%5d %5.0f %3d %4.0f %3d %4.0f %4.0f %3d %4.0f %6.0f %5.0f %3d %3d %3d %6.3f %6.3f %4d %4d\n",
          run, slope, fulldev, part, monetary, moneytraders, usingmax,
          moneygood, Csurp, Psurp, SurpSME, Nshop, NS, BS, R[0], R[1], devyear, monyear);
      fclose(stream);

              // End of the run loop
    }
    stream=fopen("evol.fil", "a");
    fprintf(stream, "\n");
    fclose(stream);
                // End of the fmean loop
  }
  stream=fopen("evol.fil", "a");
  fprintf(stream, "\n\n  n bsiz   K  xMax  lambda alpha  theta  C  Nrun Pst(yr)  T RND\n");
  fprintf(stream, "%3d %4d %3d %5d %6.3f %6.3f %6.3f %2.0f %5d %3d %6d %3d\n",
					 n, bsize, K, xMax, lambda, alpha, theta, C, numruns, persist, T, RAND);
  time(&finish);
  fprintf(stream, "\n    Started : %s", ctime(&firstbegin));
  fprintf(stream, "    Finished: %s\n", ctime(&finish));
  fclose(stream);
  return(0);
}

void Init(void) {
  time(&firstbegin);
  stream = fopen("evol.fil", "w+");
  fprintf(stream, "  run slope dev part mon mtrd usmx Mgd Csur   Psur  Esur Nsh  NS  BS  RMS_0  RMS_1  dyr  myr\n");
  fclose(stream);
  r=0; d[0]=0; s[0]=0; Pinv[0]=0.0;
  for(i=0;i<=n;i++) {
    numprod[i]=0; numcons[i]=0;
    for(j=0;j<=(n-1)*bsize;j++) {
      produces[i][j]=0;
      consumes[i][j]=0;
    }
  }
  for(i=1;i<=n;i++) for(j=1;j<=n;j++) for(k=1;k<=groupsize(i,j);k++) {
    r++; s[r]=i; d[r]=j; q[r]=(s[r]>d[r]);
    numprod[i]++; produces[i][numprod[i]]=r;
    numcons[j]++; consumes[j][numcons[j]]=r;
  }
  Fmon=bsize*(n-2)*(n-1);
  for(h=0;h<=1;h++) {vol[h][0]=0; avp[h][0]=0;}
}

int groupsize(int v, int w) {
  if(v==w) return(0);
  else return(bsize);
}

void InitRun(void) {
  printf("Number  Using  Using  Using  Using  Using  Using                \n");
  printf("Active  Money  good1  good2  good3  good4  good5   Year   NS   \n");
  endcount=0; fulldev=0; devyear=-1; monyear=-1; devcount=0;
  SEED_RANDOM(RAND+run-1);
  clock_begin = clock();
  t=0;	NS=0;
  for(r=0;r<=m;r++) {
    for(h=0;h<=1;h++) sh[h][r]=0;
    familyshop[r]=0;
  }
  for(k=0;k<=K;k++) {
    active[k]=0; owner[k]=0;
    for(h=0;h<=1;h++) {
      P[h][k]=0; y[h][k]=0; tr[h][k]=0; g[h][k]=0;
    }
  }
  lineup();
}

void lineup(void) {
  int i,j,k;
  int start[m+1];
  line[0]=0;
  start[0]=0;
  for(i=1;i<=m;i++) start[i]=i;
  for(j=1;j<=m;j++) {
    k=RANDOM_N(m-j+1)+1;
    line[j]=start[k];
    for(i=k;i<=m-j;i++) start[i]=start[i+1];
  }
}

int Research(void) {     // Assume that we have tested U(r) first. Don't need to use m1
  for(h=0;h<=1;h++) targ[h]=RANDOM_N(xMax)+1;   // or m0 after that until we get to search.
  double U;
  double P0 = F(targ[0],targ[1],f(d[r]));  // price of s[r]
  double P1 = F(targ[1],targ[0],f(s[r]));  // price of d[r]
  U=0;                              // test a comrade and a friend
  fr = comrade();
  Ucomp=Usample();
  if(d[fr]==d[r]) U = P0;
  else if(g[m1][sh[1][fr]]==d[r]) U=P0*P[m1][sh[1][fr]];
  if(U<Ucomp) {
    U=0;
    fr = soulmate();
    Ucomp=Usample();
    if(s[fr]==s[r]) U = P0;
    else if(g[m0][sh[0][fr]]==s[r]) U=P[1-m0][sh[0][fr]]*P0;
    if(U<Ucomp) return(0);
  }
  U=0;                      // test a stranger who likes s[r] and one who produces d[r]
  int k=1+RANDOM_N(numcons[s[r]]-1);
  fr= consumes[s[r]][k];
  Ucomp=Usample();
  if(s[fr]==d[r]) U=P1;
  else if(g[m0][sh[0][fr]]==d[r]) U=P[1-m0][sh[0][fr]]*P1;
  if(U<Ucomp) {
    U=0;
    k=1+RANDOM_N(numprod[d[r]]-1);
    fr= produces[d[r]][k];
    Ucomp=Usample();
    if(d[fr]==s[r]) U=P1;
    else if( g[m1][sh[1][fr]]==s[r]) U=P1*P[m1][sh[1][fr]];
    if(U<Ucomp) return(0);
  }
  return(1);
}

int comrade(void) {
  int k=1+RANDOM_N(numprod[s[r]]-1);
  int fr= produces[s[r]][k];
  if(fr>=r) fr=produces[s[r]][k+1];
  return(fr);
}

int soulmate(void) {
  int k=1+RANDOM_N(numcons[d[r]]-1);
  int fr= consumes[d[r]][k];
  if(fr>=r) fr=consumes[d[r]][k+1];
  return(fr);
}

double Usample(void) {
  m0=(g[0][sh[0][fr]]==s[fr]); m1=(g[0][sh[1][fr]]==d[fr]);
  double X=0;
  if(sh[0][fr]>0) {
    if(g[m0][sh[0][fr]]==d[fr]) X=P[1-m0][sh[0][fr]];
    else if(g[m0][sh[0][fr]]==g[m1][sh[1][fr]])
      X=P[1-m0][sh[0][fr]]*P[m1][sh[1][fr]];
  }
  return(X);
}

double Utility(void) {
  m0=(g[0][sh[0][r]]==s[r]); m1=(g[0][sh[1][r]]==d[r]);
  double X=0;
  if(sh[0][r]>0) {
    if(g[m0][sh[0][r]]==d[r])	X=P[1-m0][sh[0][r]];
    else if(g[m0][sh[0][r]]==g[m1][sh[1][r]])
      X=P[1-m0][sh[0][r]]*P[m1][sh[1][r]];
  }
  return(X);
}

void addshop(int ss) {
  if(active[ss])
  if(g[0][ss]==s[r] || g[1][ss]==s[r] || g[0][ss]==d[r] || g[1][ss]==d[r]){
    int b=0; for(a=0;a<=1+extra;a++) b+=(c[a]==ss);
    if(b==0) {extra++; c[1+extra]=ss;}
  }
}

void tryBarter(void) {
  bestbarter=0; Ubarter=0;
  for(a=2;a<=1+extra;a++)
  if(g[q[r]][c[a]]==s[r] && g[1-q[r]][c[a]]==d[r]) {
    if( P[q[r]][c[a]]> MAX(Ucomp,Ubarter)) {
      bestbarter=c[a]; Ubarter=P[q[r]][c[a]];
    }
    dropshop();
  }
}

void dropshop(void) {
  for(int b=a;b<1+extra;b++) c[b]=c[b+1];
  a--; extra--;
}

void tryOne(void) {
  for(a=2;a<=1+extra;a++)
  if(g[0][c[a]]==s[r] || g[1][c[a]]==s[r]) {
    ma=(g[0][c[a]]==s[r]);
    if(g[ma][c[a]]==g[m1][c[1]] || P[1-m0][c[0]]==0) {
      if( P[1-m0][c[0]] < P[1-ma][c[a]] ) {
        c[0]=c[a]; m0=ma;
        Ucomp = (g[ma][c[a]]==g[m1][c[1]] ? P[1-ma][c[a]]*P[m1][c[1]] : 0);
      }
      dropshop();
    }
  }
  else {
    ma=(g[0][c[a]]==d[r]);
    if(g[ma][c[a]]==g[m0][c[0]] || P[m1][c[1]]==0) {
      if( P[m1][c[1]] < P[ma][c[a]] ) {
        c[1]=c[a]; m1=ma;
        Ucomp = (g[ma][c[a]]==g[m0][c[0]] ? P[1-m0][c[0]]*P[ma][c[a]] : 0);
      }
      dropshop();
    }
  }
}

void tryTwo(void) {
  for(a=2;a<=1+extra;a++)
  if(g[0][c[a]]==s[r]||g[1][c[a]]==s[r]) {
    for(b=2;b<=1+extra;b++)
    if(a!=b && (g[0][c[b]]==d[r]||g[1][c[b]]==d[r])) {
      ma=(g[0][c[a]]==s[r]); mb=(g[0][c[b]]==d[r]);
      if(g[ma][c[a]]==g[mb][c[b]] && Ucomp < P[1-ma][c[a]]*P[mb][c[b]]) {
        Ucomp = P[1-ma][c[a]]*P[mb][c[b]];
        c[0]=c[a]; c[1]=c[b];
        m0=ma; m1=mb;
      }
    }
  }
}

int Calc1(void) {
  part=0; moneytraders=0; for(j=1;j<=n;j++) usingmoney[j]=0;
  for(r=1;r<=m;r++) if(sh[0][r]>0) {
    a=sh[0][r]; b=sh[1][r]; ma=(g[0][a]==s[r]);mb=(g[0][b]==d[r]);
    part+=( (g[ma][a]==d[r]||g[ma][a]==g[mb][b]));
    if(g[ma][a]==g[mb][b]) {
      moneytraders++;
      usingmoney[g[mb][b]]++;
    }
  }
  if(fulldev==0) {
    if(part >= .99*m) {
      if(devcount==0) devyear=t/50;
      devcount++;
    }
    else devcount=0;
    if(devcount>=persist) fulldev=1;
  }
  usingmax=0; moneygood=0;
  for(i=1;i<=n;i++) if(usingmoney[i]>usingmax) {
    usingmax=usingmoney[i];
    moneygood=i;
  }
  if(usingmoney[moneygood] >= .99*Fmon) {
    if(endcount==0) monyear=t/50;
    endcount++;
  }
  else endcount=0;
  if(endcount<persist) return(0);
  else return (1);
}

void Calc2(void) {
  if(monetary==0) monyear=-1;
  if(fulldev==0) devyear=-1;
  BS=0; for(k=1;k<=K;k++)
    if(active[k] && g[0][k]!=moneygood && g[1][k]!=moneygood) BS++;
  W=1.0-(f(moneygood)+C)/bsize;
  if(W>0.0)
    for(i=1;i<=n;i++)
      Pinv[i]=((m/n)-(f(i)+C))/((m/n)-(n-2)*(f(moneygood)+C));
  SurpSME=m-n*f1-(slope/2)*n*(n-1)-(n-2)*f(moneygood);
  RMSE();
  Csurp=0;
  for(r=1;r<=m; r++) Csurp += Utility();
  Psurp=0; Nshop=0;
  for(k=1;k<=K;k++) if(active[k]) {
    for(h=0;h<=1;h++) Psurp += (y[h][k]-f(g[h][k])-P[1-h][k]*y[1-h][k]);
    if(y[q[owner[k]]][k]>1||y[1-q[owner[k]]][k]>0) Nshop++;
  }
}

void RMSE(void) {
  if(W>0.0) {
    for(h=0;h<=1;h++) for(i=1;i<=n;i++){vol[h][i]=0; avp[h][i]=0;}
      for(k=1;k<=K;k++) if(g[0][k]==moneygood || g[1][k]==moneygood) {
      ma=(g[1][k]==moneygood); i=g[1-ma][k];
      vol[0][i]+=y[1-ma][k];
      if(vol[0][i]) avp[0][i]+=(y[1-ma][k]/vol[0][i])*(P[1-ma][k]-avp[0][i]);
      vol[1][i]+=y[ma][k];
      if(vol[1][i]) avp[1][i]+=(y[ma][k]/vol[1][i])*(P[ma][k]-avp[1][i]);
    }
    for(i=1;i<=n;i++) {
      avp[0][i]/=W;
      avp[1][i]/=Pinv[i];
    }
    for(h=0;h<=1;h++) {
      R[h]=0;
      for(i=1;i<=n;i++) if(i!=moneygood) R[h]+=(1-avp[h][i])*(1-avp[h][i]);
      R[h]/=(n-1);
      R[h]=sqrt(R[h]);
    }
  }
  else {
    for(h=0;h<=1;h++) R[h]=-1.0;
  }
}

uint32_t pcg32_next(pcg32_t *rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1ULL);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

void pcg32_seed(pcg32_t *rng, uint64_t initstate, uint64_t initseq) {
    rng->state = 0U;
    rng->inc   = (initseq << 1u) | 1u;
    pcg32_next(rng);
    rng->state += initstate;
    pcg32_next(rng);
}
