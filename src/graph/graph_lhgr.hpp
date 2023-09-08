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

#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <ranges>

#include "../misc.hpp"

/**
 *  @brief implementation of a HGR-based graph representation for fast backtracking
 *         'A Hybrid Graph Representation for Exact Graph Algorithms'
 */

// struct that contains all information required for backtracking the graph
class graph_lhgr_repr {
  public:
    //number of active vertices
    var_t no_v;
    //number of active edges
    var_t no_e;

    //degree vector -- VD[v] is degree of vertex v
    vec<var_t> VD_out;

    //color vector -- VC[v] is color of vertex v
    vec<var_t> VC;

    //ctor for graph_lhgr_repr
    graph_lhgr_repr(const var_t _no_v, const var_t _no_e, const vec<var_t>& _VD_out, const vec<var_t>& _VC) noexcept : no_v(_no_v), no_e(_no_e), VD_out(_VD_out), VC(_VC) {};
    graph_lhgr_repr(const graph_lhgr_repr& o) noexcept : no_v(o.no_v), no_e(o.no_e), VD_out(o.VD_out), VC(o.VC) {};
    graph_lhgr_repr(graph_lhgr_repr&& o) noexcept : no_v(std::move(o.no_v)), no_e(std::move(o.no_e)), VD_out(std::move(o.VD_out)), VC(std::move(o.VC)) {};
    ~graph_lhgr_repr() = default;
};

/**
 * @brief class for skew-symmetric graphs supporting quick O(d(v)) vertex merging and O(d(v)) vertex removal algorithms, and O(1) undo operations
 * 
 */
class graph_lhgr {
  protected:
    /* data */
    /* 
     *  number of (active) vertices
     */
    var_t no_v;

    /* 
     *  number of edges
     */
    var_t no_e;
    
    /*
     *  L list of vertices
     *  O( no_v )
     */
    vec<var_t> L;

    /*
     *  IL - index of vertex in L
     *  O( no_v )
     */ 
    vec<var_t> IL;
    
    /*
     *  AL - adjacency list
     *  O( no_e )
     */ 
    vec< vec<var_t> > AL_out;
    
    /*
     *  IAL - 'inverse' adjacency list
     *  O( no_e )
     *  satisfies:
     *     sigma( AL_out[ sigma(AL_out[v,i]), ILA_in[v,i] ] ) = v
     */ 
    vec< vec<var_t> > IAL_in;
    
    
    /*
     *  CAL - singly linked list of vertices by color
     *  O( no_v )
     */ 
    vec< std::list<var_t> > CAL;

    /*
     *  VC - list of vertex colors
     *  O( no_v )
     */ 
    vec<var_t> VC;

    /*
     *  CD - array of color out-degrees
     *  O( no_v )
     */ 
    vec<var_t> CD_out;

    /*
     *  VD - list of vertex out-degrees
     *  O( no_v )
     */ 
    vec<var_t> VD_out;

    /**
     * @brief skey-symmetry w.r.t sigma, we have (v,w) in E iff (SIGMA(w)(v)) is in E
     * 
     * @note this allows us to only store out-going edges for each node, as AL_in[v] = sigma( AL_out[ sigma(v) ] )
     */
    //vec<var_t> sigma;

    /**
     * @brief removes edge src->AL_out[src][idx] from the graph; does not touch its symmetrical edge!
     * 
     * @param src source of edge
     * @param idx index of target in AL_out[src]
     */
    void remove_edge_(const var_t src, const var_t idx) noexcept;

  public:
    /**
     *  construct (lean) hybrid graph representation
     *  O( no_v + no_e )
     * 
     *  @param no_v_ number of vertices
     *  @param no_e_ number of edges
     *  @param E     vector of vector of idx of edges, i.e., we have edge (u,v) iff v is in E[u],
     *               note that the vertices have to be in {0,..,no_v_-1} and there should be EXACTLY no_e_ edges in total.
     */ 
    graph_lhgr(const vec< std::pair<var_t,var_t> >& E, const var_t no_v) noexcept;
    
    graph_lhgr(const graph_lhgr& g) noexcept : no_v(g.no_v), no_e(g.no_e), L(g.L), IL(g.L), AL_out(g.AL_out), IAL_in(g.IAL_in), CAL(g.CAL), VC(g.VC), CD_out(g.CD_out), VD_out(g.VD_out) {};

    graph_lhgr() noexcept = default;
    
    ~graph_lhgr() = default;

    /**
     *  construct (lean) hybrid graph representation
     *  O( no_v + no_e )
     * 
     *  @param no_v_ number of vertices
     *  @param E     vector of vector of idx of edges, i.e., we have edge (u,v) iff v is in E[u],
     *               note that the vertices have to be in {0,..,no_v_-1} and there should be EXACTLY no_e_ edges in total.
     */ 
    void init(const vec< std::pair<var_t,var_t> >& E, const var_t no_v_) noexcept;

    /**
     * @brief adds new edge from src to dst; and its skew-symmetric edge
     * 
     * @param src source vert
     * @param dst dest vert
     * @note breaks graph_lhgr_repr produced with get_state(); i.e., after adding an edge backtrack cannot be used with older states!
     */
    void add_edge(const var_t src, const var_t dst) noexcept;
    
    /**
     * @brief Get number of (active) verts of graph
     * 
     * @return var_t number of verts
     */
    inline var_t get_no_v() const noexcept { return no_v; };
    
    /**
     * @brief Get number of (active) edges of graph
     * 
     * @return var_t number of edges
     */
    inline var_t get_no_e() const noexcept { return no_e; };

    /**
     * @brief Get range over all active verts
     * 
     * @return auto range iterating over all active verts
     */
    inline const auto get_v_range() const noexcept {
      return L | std::ranges::views::take(no_v);
    };
    
    /**
     * @brief Get vector of active vertices
     * 
     * @return vector of out-neighbours
     */
    vec<var_t> get_v_vector() const noexcept {
      vec<var_t> vs(no_v);
      std::copy(L.begin(), L.begin()+no_v, vs.begin());
      return vs;
    }

    /**
     * @brief get the a representation of the graph that allows O(1) backtracking
     * 
     * @return graph_lhgr_repr 
     */
    inline graph_lhgr_repr get_state() const noexcept { return graph_lhgr_repr(no_v, no_e, VD_out, VC); };

    /**
     *  @brief backtrack to graph represented by graph_orig
     * 
     *  @param graph_orig representation of graph to backtrack to
     */
    void backtrack(graph_lhgr_repr&& graph_orig) noexcept;

    /**
     *  @brief removes the i-th outgoing edge of src and its symmetrical counterparts
     * 
     *  @param src source vertex
     *  @param idx index of out-edge to be removed (w.r.t. order in AL_out_logical)
     *  @warning use with care! function is not compatible with vertex-removal+backtracking (as order in AL might change)!
     *  @note implementation is in O(1)
     */
    void remove_edge(const var_t src, const var_t idx) noexcept;
    
    /**
     *  @brief removes all out-going and incoming edges of v (faster than removing the edges one-by-one) and their symmetrical counterparts
     * 
     *  @param v v vertex
     *  @note implementation is in O( CD_out[v]+CD_out[sigma(v)] )
     */
    void remove_all_edges(const var_t v) noexcept;

    /**
     *  @brief removes vertex v and its symmetric counterpart from the graph
     * 
     *  @param v vertex (and its color) to be removed
     *  @note also removes all outgoing edges
     *  @note implementation is in O( CD_out[v]+CD_out[sigma(v)] )
     */
    void remove_vert(const var_t v) noexcept;

    /**
     * @brief merges vertices v1 and v2, color of v2 will be set to color of v1 (also merges the symmetric counterparts!)
     * 
     * @param v1 vertex to merge with
     * @param v2 vertex to be merged with
     * @note implementation is in O( CD_out[v1]+CD_out[sigma(v1)]+CD_out[v2]+CD_out[sigma(v2)] );
     */
    void merge_verts(const var_t v1, const var_t v2) noexcept;

    /**
     * @brief get out-degree of vertex (i.e. degree of its color!)
     * 
     * @param v vertex (color)
     * @return out-degree of vert (color)
     */
    inline var_t get_out_degree(const var_t v) const noexcept { return CD_out[ VC[v] ]; };

    /**
     * @brief Get range over out-neighbours of v 
     * 
     * @param v vertex v
     * @return range over out-neighbors
     */
    inline auto get_out_neighbour_range(const var_t v) const noexcept {
      return CAL[ VC[v] ] | std::views::transform(
                             [this](const auto& w)
                             { return AL_out[w] | std::ranges::views::take(VD_out[w]); }
                            )
                          | std::views::join
                          | std::views::transform(
                             [this](const auto& w)
                             { return VC[w]; }
                          );
    };

    /**
     * @brief Get vector of out-neighbours of vertex v
     * 
     * @param v vertex v
     * @return vector of out-neighbours
     */
    vec<var_t> get_out_neighbour_vector(const var_t v) const noexcept {
      auto r = get_out_neighbour_range(v);
      vec<var_t> out_n;
      std::ranges::copy(r, std::back_inserter(out_n));
      return out_n;
    }
    
    /**
     * @brief get in-degree of vertex (i.e. degree of its color!)
     * 
     * @param v vertex (color)
     * @return in-degree of vert (color)
     */
    inline var_t get_in_degree(const var_t v) const noexcept { return CD_out[ VC[ SIGMA(v) ] ]; };
    
    /**
     * @brief Get range over in-neighbours of v 
     * 
     * @param v vertex v
     * @return range over in-neighbors
     */
    inline auto get_in_neighbour_range(const var_t v) const noexcept {
      return CAL[ VC[SIGMA(v)] ] | std::views::transform(
                                    [this](const auto& w)
                                    { return AL_out[w] | std::ranges::views::take(VD_out[w]); }
                                   )
                                 | std::views::join
                                 | std::views::transform(
                                    [this](const auto& w)
                                    { return VC[SIGMA(w)]; }
                                 );
    };

    /**
     * @brief Get vector of in-neighbours of vertex v
     * 
     * @param v vertex v
     * @return vector of in-neighbours
     */
    vec<var_t> get_in_neighbour_vector(const var_t v) const noexcept {
      auto r = get_in_neighbour_range(v);
      vec<var_t> in_n;
      std::ranges::copy(r, std::back_inserter(in_n));
      return in_n;
    }

    /**
     * @brief assert that internal repr of graph is consistent
     * 
     * @return true iff no checks fail; otherwise throws an assertion exception
     */
    virtual bool assert_data_structs() const noexcept;


    //get unique string repr of graph! (list of lexicographically ordered edges!)
    virtual std::string to_str() const noexcept;

    inline bool operator ==(const graph_lhgr &other) const noexcept { return to_str() == other.to_str(); };

    graph_lhgr& operator=(graph_lhgr& g) {
      no_v = g.no_v;
      no_e = g.no_e;
      L = g.L;
      IL = g.IL;
      AL_out = g.AL_out;
      IAL_in = g.IAL_in;
      CAL = g.CAL;
      VC = g.VC;
      CD_out = g.CD_out;
      VD_out = g.VD_out;

      return *this;
    }
};