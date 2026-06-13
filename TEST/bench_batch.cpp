// ============================================================================
// TEST/bench_batch.cpp
//
// BATCH evaluation: one polynomial structure, many (default 7200) coefficient
// sets -- the kinematic-scan regime the tropical pipeline is built for.  The
// fan + flattened sectors are computed ONCE (gen_kin.wl); here we ask how best
// to push 7200 coefficient sets through them.
//
// Spec: P = c0 + c1 x1 + c2 x2 + c3 x3 + c4 x4,  A=0,  B=-6,  n=4.
// Exact per coeff set:  I = 1 / (120 * c0^2 * c1 c2 c3 c4).
//
// Strategies compared (all single-threaded for fair timing):
//   1 MC_perkp        plain MC, independent per kp           (the shipped model)
//   2 QMC_shared      Sobol points drawn once, shared across all kp
//   3 QMC_shared_basis  + compute the coefficient-INDEPENDENT monomial basis
//                       once per sample, reuse for all kp     (the amortization)
//   4 Vegas_perkp     CUBA Vegas, one call per (kp,sector)    (naive CUBA batch)
//   5 Vegas_ncomp     CUBA Vegas, ONE call/sector with ncomp=2*Nkp (shared samples)
//   6 Cuhre_ncomp     CUBA Cuhre, ncomp=2*Nkp
//
// Reports wall time, throughput (kp/s), and accuracy (median/max rel err vs
// the closed form) at a fixed per-sector budget.
//
// Build: g++ -O3 -std=c++17 -I/opt/homebrew/include -DSECTOR_HEADER='"..."' \
//        bench_batch.cpp -L/opt/homebrew/lib -lcuba -lm -o bench_batch
// ============================================================================
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <complex>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <boost/random/sobol.hpp>
extern "C" {
#include <cuba.h>
}
#include SECTOR_HEADER   // integrand_table[], integrand_dim[], N_INTEGRANDS,
                         // N_PARAMS, MAX_DIM, KS_* shared-basis tables

using clk = std::chrono::steady_clock;
static double secs(clk::time_point t0){ return std::chrono::duration<double>(clk::now()-t0).count(); }
static inline double u01(uint64_t u){ return (double)(u>>11)*(1.0/9007199254740992.0); }

static int     NKP = 7200;
static int     MS  = 2000;          // samples / maxeval per sector
static std::vector<double> PAR;     // NKP * N_PARAMS
static std::vector<double> EXACT;   // NKP

static void make_coeffs(){
  PAR.assign((size_t)NKP*N_PARAMS,0.0);
  EXACT.assign(NKP,0.0);
  std::mt19937_64 g(123);
  std::uniform_real_distribution<double> ud(-0.7,0.7);   // c in [~0.5, ~2]
  double fact=1.0; for(int k=2;k<=CASE_DIM+1;++k) fact*=k;   // (n+1)!
  for(int kp=0;kp<NKP;++kp){
    double* c=&PAR[(size_t)kp*N_PARAMS];
    c[0]=1.0;                                   // c0
    for(int p=1;p<N_PARAMS;++p) c[p]=std::exp(ud(g));
    // exact = 1 / ( (n+1)! * c0^2 * prod_{i=1}^n c_i )
    double prod=c[0]*c[0]; for(int p=1;p<N_PARAMS;++p) prod*=c[p];
    EXACT[kp]=1.0/(fact*prod);
  }
}

static void accuracy(const char* tag,const std::vector<double>& est,double wall){
  std::vector<double> rel(NKP);
  double maxr=0,meanr=0;
  for(int kp=0;kp<NKP;++kp){ rel[kp]=std::fabs(est[kp]-EXACT[kp])/std::fabs(EXACT[kp]);
    maxr=std::max(maxr,rel[kp]); meanr+=rel[kp]; }
  meanr/=NKP;
  std::vector<double> s=rel; std::sort(s.begin(),s.end());
  double med=s[NKP/2];
  printf("%-18s  wall=%8.3f s   thru=%9.0f kp/s   relerr med=%.2e mean=%.2e max=%.2e\n",
         tag,wall,NKP/wall,med,meanr,maxr);
}

// ---------------- 1. MC per kp (the shipped model) -----------------------
static void run_mc_perkp(){
  auto t0=clk::now();
  std::vector<double> est(NKP,0.0);
  std::mt19937_64 rng(42);
  std::uniform_real_distribution<double> U(0,1);
  double y[MAX_DIM];
  for(int kp=0;kp<NKP;++kp){
    const double* pr=&PAR[(size_t)kp*N_PARAMS];
    double tot=0;
    for(int s=0;s<N_INTEGRANDS;++s){
      int d=integrand_dim[s]; double mean=0;
      for(int k=0;k<MS;++k){
        for(int i=0;i<d;++i) y[i]=U(rng);
        mean+=(integrand_table[s](y,pr).real()-mean)/(k+1);
      }
      tot+=mean;
    }
    est[kp]=tot;
  }
  accuracy("1.MC_perkp",est,secs(t0));
}

// ---------------- 2. QMC, Sobol points shared across kp ------------------
static void run_qmc_shared(){
  auto t0=clk::now();
  std::vector<double> est(NKP,0.0);
  std::mt19937_64 sg(7); std::uniform_real_distribution<double> U(0,1);
  double y[MAX_DIM], shift[MAX_DIM];
  for(int s=0;s<N_INTEGRANDS;++s){
    int d=integrand_dim[s];
    for(int i=0;i<d;++i) shift[i]=U(sg);
    boost::random::sobol eng(d);
    std::vector<double> acc(NKP,0.0);
    for(int m=0;m<MS;++m){
      for(int i=0;i<d;++i){ double v=u01(eng())+shift[i]; v-=std::floor(v); y[i]=v; }
      for(int kp=0;kp<NKP;++kp)
        acc[kp]+=integrand_table[s](y,&PAR[(size_t)kp*N_PARAMS]).real();
    }
    for(int kp=0;kp<NKP;++kp) est[kp]+=acc[kp]/MS;
  }
  accuracy("2.QMC_shared",est,secs(t0));
}

// -------- 3. QMC shared samples + shared (coefficient-independent) basis ----
static void run_qmc_shared_basis(){
  auto t0=clk::now();
  std::vector<double> est(NKP,0.0);
  std::mt19937_64 sg(7); std::uniform_real_distribution<double> U(0,1);
  double logy[MAX_DIM], shift[MAX_DIM];
  for(int s=0;s<KS_NSEC;++s){
    int d=KS_DIM, m0=KS_MOFF[s], m1=KS_MOFF[s+1], nm=m1-m0;
    double pref=KS_PREFAC[s];
    for(int i=0;i<d;++i) shift[i]=U(sg);
    boost::random::sobol eng(d);
    std::vector<double> basis(nm);
    std::vector<double> acc(NKP,0.0);
    for(int m=0;m<MS;++m){
      for(int i=0;i<d;++i){ double v=u01(eng())+shift[i]; v-=std::floor(v);
        logy[i]=(v>1e-300)?std::log(v):-700.0; }
      // coefficient-INDEPENDENT monomial basis: computed ONCE for all kp
      for(int j=0;j<nm;++j){
        const double* a=&KS_ALPHA[(size_t)(m0+j)*KS_DIM];
        double e=0; for(int i=0;i<d;++i) e+=a[i]*logy[i];
        basis[j]=std::exp(e);
      }
      // per kp: cheap dot-product + one pow
      for(int kp=0;kp<NKP;++kp){
        const double* pr=&PAR[(size_t)kp*N_PARAMS];
        double P=0;
        for(int j=0;j<nm;++j){
          const double* cr=&KS_CROW[(size_t)(m0+j)*KS_NPARAMS];
          double c=KS_CCONST[m0+j];
          for(int p=0;p<KS_NPARAMS;++p) c+=cr[p]*pr[p];
          P+=c*basis[j];
        }
        acc[kp]+=pref*std::pow(P,KS_B);
      }
    }
    for(int kp=0;kp<NKP;++kp) est[kp]+=acc[kp]/MS;
  }
  accuracy("3.QMC_shared_basis",est,secs(t0));
}

// ---------------- CUBA plumbing ----------------
static const double* g_par=nullptr;  static int g_sector=0, g_dim=0;
static int g_kp0=0, g_chunk=0;          // ncomp batching processes a kp-chunk
static const int CUBA_MAXCOMP=512;      // safe ncomp chunk (the usable limit
                                        // shrinks with dim: 1024 crashes Vegas
                                        // at dim 8; 512 is safe through dim 8)
// per-kp integrand (Vegas_perkp): ncomp=1 (integrand is real here)
static int cb_one(const int*,const cubareal x[],const int*,cubareal f[],void*){
  double y[MAX_DIM]; for(int i=0;i<g_dim;++i) y[i]=(double)x[i];
  f[0]=integrand_table[g_sector](y,g_par).real(); return 0;
}
// ncomp-batched integrand: ncomp = chunk size, shared sample point across the
// chunk's kp (component c -> kinematic point g_kp0+c)
static int cb_ncomp(const int*,const cubareal x[],const int*,cubareal f[],void*){
  double y[MAX_DIM]; for(int i=0;i<g_dim;++i) y[i]=(double)x[i];
  for(int c=0;c<g_chunk;++c)
    f[c]=integrand_table[g_sector](y,&g_par[(size_t)(g_kp0+c)*N_PARAMS]).real();
  return 0;
}

static void run_vegas_perkp(){
  auto t0=clk::now();
  std::vector<double> est(NKP,0.0);
  cubareal integ[1],err[1],prob[1]; int neval,fail;
  for(int kp=0;kp<NKP;++kp){
    g_par=&PAR[(size_t)kp*N_PARAMS];
    double tot=0;
    for(int s=0;s<N_INTEGRANDS;++s){
      g_sector=s; g_dim=integrand_dim[s];
      Vegas(g_dim,1,cb_one,nullptr,1, 1e-4,1e-12, 0,0, 0,MS, 1000,500,1000,
            0,nullptr,nullptr,&neval,&fail,integ,err,prob);
      tot+=integ[0];
    }
    est[kp]=tot;
  }
  accuracy("4.Vegas_perkp",est,secs(t0));
}

static void run_cuba_ncomp(bool cuhre){
  auto t0=clk::now();
  std::vector<double> est(NKP,0.0);
  std::vector<cubareal> integ(CUBA_MAXCOMP),err(CUBA_MAXCOMP),prob(CUBA_MAXCOMP);
  int neval,fail,nregions;
  // fair per-kp budget: MS shared sample points/sector (each serves all kp in
  // the chunk), so every kp gets an MS-sample estimate, same as MC/QMC.
  int maxev=MS;
  g_par=PAR.data();
  // CUBA's ncomp is capped (MAXCOMP=1024), so batch kp in chunks of <=1024.
  for(int s=0;s<N_INTEGRANDS;++s){
    g_sector=s; g_dim=integrand_dim[s];
    for(int k0=0;k0<NKP;k0+=CUBA_MAXCOMP){
      int cs=std::min(CUBA_MAXCOMP,NKP-k0);
      g_kp0=k0; g_chunk=cs;
      if(cuhre)
        Cuhre(g_dim,cs,cb_ncomp,nullptr,1, 1e-4,1e-12, 0, 0,maxev, 0,
              nullptr,nullptr,&nregions,&neval,&fail,integ.data(),err.data(),prob.data());
      else
        Vegas(g_dim,cs,cb_ncomp,nullptr,1, 1e-4,1e-12, 0,0, 0,maxev,
              1000,500,1000, 0,nullptr,nullptr,&neval,&fail,integ.data(),err.data(),prob.data());
      for(int c=0;c<cs;++c) est[k0+c]+=integ[c];
    }
  }
  accuracy(cuhre?"6.Cuhre_ncomp":"5.Vegas_ncomp",est,secs(t0));
}

int main(int argc,char**argv){
  setvbuf(stdout,nullptr,_IONBF,0);   // unbuffered: don't lose lines on crash
  if(argc>1) NKP=atoi(argv[1]);
  if(argc>2) MS =atoi(argv[2]);
  const int zero=0; cubacores(&zero,&zero);
  make_coeffs();
  printf("BATCH test: Nkp=%d coefficient sets, dim=%d, %d sectors, "
         "%d samples(or maxeval)/sector\n", NKP, CASE_DIM, N_INTEGRANDS, MS);
  printf("Spec P=c0+sum c_i x_i, B=-6; exact=1/(120 c0^2 prod c_i). "
         "Single-threaded.\n");
  printf("%s\n", std::string(96,'-').c_str());
  run_mc_perkp();
  run_qmc_shared();
  run_qmc_shared_basis();
  run_vegas_perkp();
  run_cuba_ncomp(false);
  run_cuba_ncomp(true);
  return 0;
}
