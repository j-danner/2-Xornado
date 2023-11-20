// Copyright (c) 2022-2023 Julian Danner <julian.danner@uni-passau.de>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
// the Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

//std
#include <atomic>
#include <stdint.h>
#include <cassert>
#include <chrono>
#include <vector>
#include <iostream>
#include <iomanip>
#include <unordered_map>
//other
#include <omp.h>
#include "robin_hood-3.11.5/robin_hood.h"


//verbosity output
#ifdef VERBOSITY
  #define VERB(lvl, msg) if(opt.verb >= lvl) { std::cout << msg << std::endl; }
#else
  #define VERB(lvl, msg)
#endif

#define SIGMA(i) (var_t)(((var_t) i) ^ ((var_t) 1))


//type for variable numbering (16bit should suffice)
typedef uint16_t var_t;
//typedef uint_fast16_t var_t;


//select vector impl to use
template<class T>
using vec = std::vector<T>;


enum class dec_heu { fv, mp, mr, mbn, lex};
enum class fls_alg { no, trivial, trivial_cc, full};
enum class upd_alg { ts, hf, par, hfd};
enum class sc {active, inactive};
enum class constr { simple, extended};
enum class preproc { no, scc, fls_scc, fls_scc_ee };

/**
 * @brief class that handles reordering according to guessing path
 */
class reordering {
  private:
    //TODO use faster hashmap
  #ifdef NDEBUG
    robin_hood::unordered_flat_map<var_t,var_t> P;
  #else
    std::unordered_map<var_t,var_t> P;
  #endif

  public:
    reordering() {};
    reordering(const reordering& o) : P(o.P) {};
    reordering(reordering&& o) : P(std::move(o.P)) {};

    std::size_t size() const noexcept { return P.size(); };

    void insert(const var_t& ind, const var_t& pos) {
      if(at(pos)==ind) return;
      const auto P_ind = at(ind);
      const auto P_pos = at(pos);
      P[pos] = P_ind;
      P[ind] = P_pos;
    };
    const var_t& at(const var_t& ind) const noexcept { return P.contains(ind) ? P.at(ind) : ind; };
};


/**
 * @brief struct that holds options for the various heuristic choices
 * 
 */
struct options {
    var_t num_vars = 0;
    var_t num_cls = 0;

    dec_heu dh = dec_heu::mp;
    fls_alg fls = fls_alg::no;
    int fls_s = 1;
    upd_alg upd = upd_alg::ts;
    sc score = sc::inactive;
    constr ext = constr::extended;
    preproc pp = preproc::no;
    
    int jobs = omp_get_num_threads();
    
    int verb = 0;

    int timeout = 0;

    reordering P;

    //TODO selected heuristics!

    //default settings
    options() : num_vars(0), num_cls(0) {};
    options(var_t n_vars) : num_vars(n_vars), num_cls(0) {};
    options(var_t n_vars, var_t n_cls) : num_vars(n_vars), num_cls(n_cls) {};
    options(var_t n_vars, var_t n_cls, dec_heu dh_, fls_alg fls_, upd_alg upd_, int jobs_, int verb_, int timeout_) : num_vars(n_vars), num_cls(n_cls), dh(dh_), fls(fls_), upd(upd_), jobs(jobs_), verb(verb_), timeout(timeout_) {};
    options(var_t n_vars, var_t n_cls, dec_heu dh_, fls_alg fls_, int fls_s_, upd_alg upd_, sc score_, constr ext_, preproc pp_, int jobs_, int verb_, int timeout_) : num_vars(n_vars), num_cls(n_cls), dh(dh_), fls(fls_), fls_s(fls_s_), upd(upd_), score(score_), ext(ext_), pp(pp_), jobs(jobs_), verb(verb_), timeout(timeout_) {};
    options(var_t n_vars, var_t n_cls, dec_heu dh_, fls_alg fls_, int fls_s_, upd_alg upd_, sc score_, constr ext_, preproc pp_, int jobs_, int verb_, int timeout_, reordering P_) : num_vars(n_vars), num_cls(n_cls), dh(dh_), fls(fls_), fls_s(fls_s_), upd(upd_), score(score_), ext(ext_), pp(pp_), jobs(jobs_), verb(verb_), timeout(timeout_), P(P_) {};
};


/**
 * @brief struct returned by solver, contains bool sat telling if intance is satisfiable; if it is, also contains a solution
 * 
 */
class stats {
  public:
    bool finished = false;
    bool sat = false;
    vec<bool> sol;
    std::atomic<bool> cancelled = false;

    unsigned long no_dec = 0;
    unsigned long no_confl = 0;
    unsigned long no_vert_upd = 0;
    unsigned long no_restarts = 0;
    unsigned long no_graph_upd = 0;
    unsigned long no_crGCP = 0;
    unsigned long total_upd_no_v = 0;
    unsigned long total_upd_xsys_size = 0;
    //newly llongnt pure-xors via scc, fls, upd
    unsigned long new_px_scc = 0;
    unsigned long new_px_fls = 0;
    unsigned long new_px_upd = 0;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    void print_stats() const {
      std::cout << "c v_upd     : " << no_vert_upd << std::endl;
      std::cout << "c crGCP     : " << no_crGCP << std::endl;
      std::cout << "c restarts  : " << no_restarts << std::endl;
      std::cout << "c decisions : " << no_dec << std::endl;
      std::cout << "c conflicts : " << no_confl << std::endl;
    };

    void print_final() const {
      float total_time = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0f;
      std::cout << std::fixed << std::setprecision(3);

      std::cout << "c dec/sec    : "  << no_dec/total_time << std::endl;
      std::cout << "c v_upd/sec  : " << no_vert_upd/total_time << std::endl;
      std::cout << "c " << std::endl;
      std::cout << "c v_upd/dec  : " << ((double) no_vert_upd)/((double) no_dec) << std::endl;
      std::cout << "c " << std::endl;
      std::cout << "c avg graph size : " << ((double) total_upd_no_v)/((double) no_graph_upd) << std::endl;
      std::cout << "c avg LinEqs size  : " << ((double) total_upd_xsys_size)/((double) no_graph_upd) << std::endl;
      std::cout << "c " << std::endl;

      std::cout << "c lins from upd  : " << new_px_upd << std::endl;
      std::cout << "c lins from SCC  : " << new_px_scc << std::endl;
      std::cout << "c lins from FLS  : " << new_px_fls << std::endl;
      std::cout << "c " << std::endl;

      std::cout << "c vertex upd : " << no_vert_upd << std::endl;
      std::cout << "c graph upd  : " << no_graph_upd << std::endl;
      std::cout << "c crGCP      : " << no_crGCP << std::endl;
      //std::cout << "c restarts   : " << no_restarts << std::endl;
      std::cout << "c decisions  : " << no_dec << std::endl;
      std::cout << "c conflicts  : " << no_confl << std::endl;
      std::cout << "c Total time : " << total_time << " [s]" << std::endl;

      print_sol();
    }
    
    void reorder_sol(const reordering& P) {
      if(sol.size()==0) return;
      vec<bool> Psol(sol);
      for(var_t i=1; i <= sol.size(); ++i) {
        Psol[i-1] = sol[P.at(i)-1];
      }
      sol = std::move(Psol);
    }

    void print_sol() const {
      if(finished) {
          if(sat) {
              std::cout << "s SATISFIABLE" << std::endl;
              std::cout << "v ";
              for (var_t i = 1; i <= sol.size(); i++) {
                  std::cout << (sol[i-1] ? "" : "-") << std::to_string( i ) << " ";
              }
              std::cout << "0" << std::endl;
          } else {
              std::cout << "s UNSATISFIABLE" << std::endl;
          }
      } else {
              std::cout << "c timeout reached or interupted!" << std::endl;
              std::cout << "s INDEFINITE" << std::endl;
      }
    };
    
    stats() {};
    ~stats() { /*std::cout << "destroying stats!" << std::endl;*/ };
    stats(stats& o) noexcept : finished(o.finished), sat(o.sat), sol(o.sol), no_dec(o.no_dec), no_confl(o.no_confl), no_vert_upd(o.no_vert_upd), no_restarts(o.no_restarts), new_px_scc(o.new_px_scc), new_px_fls(o.new_px_fls), new_px_upd(o.new_px_upd), begin(o.begin), end(o.end) {
      cancelled.store( o.cancelled.load() );
    }
    stats(stats&& o) noexcept : finished(std::move(o.finished)), sat(std::move(o.sat)), sol(std::move(o.sol)), no_dec(std::move(o.no_dec)), no_confl(std::move(o.no_confl)), no_vert_upd(std::move(o.no_vert_upd)), no_restarts(std::move(o.no_restarts)), new_px_scc(std::move(o.new_px_scc)), new_px_fls(std::move(o.new_px_fls)), new_px_upd(std::move(o.new_px_upd)), begin(std::move(o.begin)), end(std::move(o.end))  {
      cancelled.store( o.cancelled.load() );
    }
    stats(unsigned int no_dec_, unsigned int no_confl_, const vec<bool>& sol_) : sat(true), sol(sol_), no_dec(no_dec_), no_confl(no_confl_) {};
    stats(unsigned int no_dec_, unsigned int no_confl_) : sat(false), no_dec(no_dec_), no_confl(no_confl_) {};
};

