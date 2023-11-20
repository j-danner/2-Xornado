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

#define FULL_REDUCTION

#include <stack>
#include <math.h>
#include <map>
#include <queue>
#include <list>
#include <memory>

#include "solve.hpp"
#include "misc.hpp"
#include "graph/graph.hpp"
#include "LA/lineral.hpp"
#include "LA/lineqs.hpp"

#define Lsys xsys_stack.back()
#define linsys xsys_stack.back().back()

#include "vl/vl.hpp"

/**
 * @brief class for implication graph structures offering SCC-analysis and FL-search combined with quick O(1) backtracking
 * 
 */
class impl_graph : public graph
{
  private:
    vert_label vl;
    /**
     * @brief stack of maps from vertices to linerals
     */
    std::stack< vert_label_repr > vl_stack;

    /**
     * @brief stack of graph_repr for backtracking
     */
    std::stack< graph_repr > graph_stack;

    /**
     * @brief stack of lists of xsyses for backtracking
     */
    std::list< std::list<LinEqs> > xsys_stack;

#ifndef FULL_REDUCTION
    /**
     * @brief current assignments
     */
    vec<lineral> assignments;
#endif

    /**
     * @brief options for heuristics of dpll-solver (and more)
     */
    options opt;

    /**
     * @brief solving stats
     */
    stats s;

    /**
     * @brief 'activity' of each variable; used for decision heuristic
     * @note entries must be strictly positive! (otherwise max_path/max_reach might fail!)
     */
    vec<unsigned int> activity_score;
    unsigned int bump = 1;
    float decay = 0.9;

    //SCC-helper funcs
    void scc_dfs_util(const var_t rt, vec<lineral>& linerals, var_t v, vec<bool>& visited, std::list< std::pair<var_t,var_t> >& merge_list) const;
    void scc_fillOrder(const var_t v, vec<bool>& visited, std::stack<var_t> &Stack) const;

    typedef LinEqs (impl_graph::*upd_t)(stats& s, const LinEqs&);
    typedef LinEqs (impl_graph::*fls_t)() const;
    typedef std::pair<LinEqs,LinEqs> (impl_graph::*dec_heu_t)() const;

    void crGCP(stats& s, const upd_t upd, const fls_t fls, const bool scheduled_fls);
    void crGCP(stats& s, const upd_t upd = &impl_graph::update_graph, const fls_t fls = &impl_graph::fls_no ) { crGCP(s,upd,fls, true); };
    void crGCP_no_schedule(stats& s, const upd_t upd = &impl_graph::update_graph, const fls_t fls = &impl_graph::fls_no ) { crGCP(s,upd,fls, false); };

    void bump_score(const LinEqs& new_xsys);
    void decay_score();

    /**
     * @brief computes all linerals implied by lit (calls crGCP with fls_no! modifies graph but backtracks afterwards!)
     * 
     * @param lit lineral that is assumed to be true
     * @return LinEqs system of implied linerals
     */
    LinEqs implied_xlits(lineral& lit) {
      //(1) save state
      auto g_state = std::move(get_state());
      auto vl_state = std::move(vl.get_state());
      xsys_stack.emplace_back( std::list<LinEqs>() );
      add_new_xsys( lit );

      //(2) call crGCP
      stats s;
      crGCP_no_schedule(s);
      //sum over all list els in xsys_stack.top
      LinEqs implied_lits;
      for(const auto& sys : xsys_stack.back()) implied_lits += sys;

      //(3) backtrack state
      vl.backtrack( std::move(vl_state), vl_stack.size() );
      //revert assignments
      xsys_stack.pop_back();
      backtrack( std::move(g_state) );
      assert( assert_data_structs() );

      return implied_lits;
    };

    //memory-friendly sum of two linerals
    inline lineral Vxlit_sum(const var_t v1, const var_t v2) const {
      bool v1_contained = vl.contains(v1);
      bool v2_contained = vl.contains(v2);

      if(v1_contained) {
        if(v2_contained) {
          return std::move( vl.sum(v1,v2) );
        } else {
          lineral tmp = std::move( vl.sum(v1,SIGMA(v2)) );
          tmp.add_one();
          return tmp;
        }
      } else {
        if(v2_contained) {
          lineral tmp = std::move( vl.sum(SIGMA(v1),v2) );
          tmp.add_one();
          return tmp;
        } else {
          return std::move( vl.sum(SIGMA(v1),SIGMA(v2)) );
        }
      }
    }

  public:
    /**
     * @brief Construct a new impl graph where each vector in clss represents a xor-clause; they must be of length at most two!
     * 
     * @param clss vector of lineral-vectors that represent the clauses
     * @param opt_ options for heuristics, also includes number of vars
     */
    impl_graph(const vec< vec<lineral> >& clss, const options& opt_);

    /**
     * @brief Construct a new impl graph
     * 
     * @param parsed_xnf pair of options and clauses, as returned by parse_file
     */
    impl_graph(parsed_xnf& p_xnf) : impl_graph(p_xnf.cls, options(p_xnf.num_vars, p_xnf.num_cls)) {};

    //copy ctor
    impl_graph(const impl_graph& ig) : graph(ig), vl(ig.vl), vl_stack(ig.vl_stack), graph_stack(ig.graph_stack), xsys_stack(ig.xsys_stack), opt(ig.opt), activity_score(ig.activity_score) {};

    ~impl_graph() = default;

    inline void add_new_xsys(const LinEqs& L) {
    #ifndef FULL_REDUCTION
      //update assignments
      for(const auto [lt,idx] : L.get_pivot_poly_idx()) {
          assert( assignments[lt].is_zero() );
          assignments[lt] = L.get_linerals(idx);
      }
    #endif
      xsys_stack.back().emplace_back(L);
    };

    //various implementations to choose from for the core update functions!
    LinEqs update_graph(stats& s, const LinEqs& L);

    LinEqs update_graph_par(stats& s, const LinEqs& L);

    LinEqs update_graph_hash_fight(stats& s, const LinEqs& L);

    LinEqs update_graph_hash_fight_dev(stats& s, const LinEqs& L);
    
    //wrappers for update_funcs when stats are irrelevant
    inline LinEqs update_graph(const LinEqs& L) { return std::move( update_graph(s, L) ); };
    inline LinEqs update_graph_par(const LinEqs& L) { return std::move( update_graph_par(s, L) ); };
    inline LinEqs update_graph_hash_fight(const LinEqs& L) { return std::move( update_graph_hash_fight(s, L) ); };
    inline LinEqs update_graph_hash_fight_dev(const LinEqs& L) { return std::move( update_graph_hash_fight_dev(s, L) ); };

    //in-processing
    LinEqs scc_analysis();

    /**
     * @brief compute roots of graph
     * 
     * @return std::list<var_t> contains all verts with indegree 0
     */
    std::list<var_t> get_roots() const;

    
    /**
     * @brief finds all connected components, and assigns all verts in one comp the same label; in each comp there is one vert with label[IL[v]]=v
     * 
     * @return vec<var_t> labels; vert v has label label[IL[v]]
     * @note assumes graph is DAG
     */
    vec<var_t> label_components() const;

    /**
     * @brief returns the number of connected components
     * 
     * @return var_t number of connected components
     */
    var_t get_number_connected_components() const;

    /**
     * @brief compute topological ordering of vertices
     * 
     * @return vec<var_t> TO of vertices; empty iff graph has a cycle (or no_v==0)
     */
    vec<var_t> get_TO() const;
    
    /**
     * @brief checks if graph is DAG
     * 
     * @return true iff graph is DAG
     */
    bool is_DAG() const { return (no_v==0) || !get_TO().empty(); };

    /**
     * @brief checks whether there is a path v->w via DFS 
     * 
     * @param v src
     * @param w dst
     * @return true iff there is a path v->w
     */
    bool is_descendant(const var_t v, const var_t w);

    std::string graph_stats() const { 
      if(opt.verb < 120) {
        return "c graph stats: #V "+std::to_string(no_v)+" #E "+std::to_string(no_e)+", #roots "+std::to_string(get_roots().size())+", #CC "+std::to_string(get_number_connected_components())+", "+(is_DAG() ? "DAG" : "no DAG")+", "+(linsys.is_consistent() ? "consistent" : "inconsistent");
      } else {
        return "c graph stats: #V "+std::to_string(no_v)+" #E "+std::to_string(no_e)+", #roots "+std::to_string(get_roots().size())+", #CC "+std::to_string(get_number_connected_components())+", "+(is_DAG() ? "DAG" : "no DAG")+", "+(linsys.is_consistent() ? "consistent" : "inconsistent")
        + "\n" + to_str();
      }
    };


    LinEqs fls_no() const;
    LinEqs fls_trivial() const; 
    LinEqs fls_trivial_cc() const;
    LinEqs fls_full() const;
    //currently unused, as computationally expensive
    LinEqs fls_full_implied();

    //preprocess
    inline void preprocess();

    //decision heuristics

    /**
     * @brief branch on first vertex (i.e. vert at first position in L)
     */
    std::pair< LinEqs, LinEqs > first_vert() const;

    /**
     * @brief branch on largest tree, i.e., guess the whole to be correct
     */
    std::pair< LinEqs, LinEqs > max_reach() const;

    /**
     * @brief branch on making the largest bottleneck
     */
    std::pair< LinEqs, LinEqs > max_bottleneck() const;

    /**
     * @brief branch on lexicographically next un-assigned idx
     */
    std::pair< LinEqs, LinEqs > lex() const;

    /**
     * @brief branch on making the longest path a cycle
     * @note if FULL_REDUCTION is not defined, and guess was already previously made, it relies on max_reach heuristic!
     */
    std::pair< LinEqs, LinEqs > max_path() const;

    /**
     * @brief branch on making the path of highest score a cycle; note: longer paths are preferred!
     */
    std::pair< LinEqs, LinEqs > max_score_path() const;

    //solve-main
    stats dpll_solve() { return dpll_solve(s); };
    stats dpll_solve(stats& s);

    inline var_t get_dl() const { return graph_stack.size()-1; };
    
    vec< vec<lineral> > to_xcls() const;

    std::string to_xnf_string() const;
    std::string to_str() const noexcept override;

    std::string to_str_() const;

    impl_graph& operator=(const impl_graph& ig) = delete;

    bool assert_data_structs() const noexcept;

    options* get_opts() { return &opt; };
    const options* get_const_opts() const { return &opt; };
    
  #ifndef FULL_REDUCTION
    void print_assignments(std::string lead = "") const noexcept {
      VERB(80, lead);
      VERB(80, lead+" assignments");
      VERB(80, lead+" lt ass");
      for(var_t i = 0; i<assignments.size(); ++i) {
          VERB(80, lead+" " + std::to_string(i) + " " + assignments[i].to_str());
      }
    };
  #endif
};