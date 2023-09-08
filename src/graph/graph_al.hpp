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
#include <set>

#include "../misc.hpp"

//#include <unordered_set>

//type of 'adjacency_list'
typedef std::set<var_t> adj_l;
//typedef std::unordered_set<var_t> adj_l; //SIGNIFICANTLY SLOWER!

/**
 *  @brief implementation of skew-symmetric graph representation
 */

// struct that contains all information required for backtracking the graph
class graph_al_repr {
  public:
    //number of active vertices
    var_t no_v;
    //number of active edges
    var_t no_e;

    //adjacency list
    vec<adj_l> AL_out;

    //ctor for graph_al_repr
    graph_al_repr(const var_t _no_v, const var_t _no_e, const vec<adj_l>& _AL_out) noexcept : no_v(_no_v), no_e(_no_e), AL_out(_AL_out) {};
    graph_al_repr(const graph_al_repr& o) noexcept : no_v(o.no_v), no_e(o.no_e), AL_out(o.AL_out) {};
    graph_al_repr(graph_al_repr&& o) noexcept : no_v(std::move(o.no_v)), no_e(std::move(o.no_e)), AL_out(std::move(o.AL_out)) {};
    ~graph_al_repr() = default;
};

/**
 * @brief class for skew-symmetric graphs supporting quick O(d(v)) vertex merging and O(d(v)) vertex removal algorithms, and O(1) undo operations
 * 
 */
class graph_al {
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
    vec<adj_l> AL_out;

    /**
     * @brief skey-symmetry w.r.t sigma, we have (v,w) in E iff (SIGMA(w)(v)) is in E
     * 
     * @note this allows us to only store out-going edges for each node, as AL_in[v] = sigma( AL_out[ sigma(v) ] )
     */
    //vec<var_t> sigma;

  public:
    /**
     *  construct graph
     *  O( no_v + no_e )
     * 
     *  @param no_v_ number of vertices
     *  @param no_e_ number of edges
     *  @param E     vector of vector of idx of edges, i.e., we have edge (u,v) iff v is in E[u],
     *               note that the vertices have to be in {0,..,no_v_-1} and there should be EXACTLY no_e_ edges in total.
     */ 
    graph_al(const vec< std::pair<var_t,var_t> >& E, const var_t no_v) noexcept;

    graph_al(const graph_al& g) noexcept : no_v(g.no_v), no_e(g.no_e), L(g.L), IL(g.L), AL_out(g.AL_out) {};

    graph_al() = default;

    ~graph_al() = default;

    /**
     *  construct graph
     *  O( no_v + no_e )
     * 
     *  @param no_v_ number of vertices
     *  @param no_e_ number of edges
     *  @param E     vector of vector of idx of edges, i.e., we have edge (u,v) iff v is in E[u],
     *               note that the vertices have to be in {0,..,no_v_-1} and there should be EXACTLY no_e_ edges in total.
     */ 
    void init(const vec< std::pair<var_t,var_t> >& E, const var_t no_v_) noexcept;


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
     * @return graph_al_repr 
     */
    inline graph_al_repr get_state() const noexcept { return graph_al_repr(no_v, no_e, AL_out); }; //TODO!!

    /**
     *  @brief backtrack to graph represented by graph_orig
     * 
     *  @param graph_orig representation of graph to backtrack to
     */
    void backtrack(graph_al_repr&& graph_orig) noexcept;

    /**
     *  @brief removes the i-th outgoing edge of src and its symmetrical counterparts
     * 
     *  @param src source vertex
     *  @param dst dst vertex
     */
    void remove_edge(const var_t src, const var_t dst) noexcept;

    /**
     *  @brief removes all out-going and incoming edges of v (faster than removing the edges one-by-one) and their symmetrical counterparts
     * 
     *  @param v v vertex
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
     */
    void merge_verts(const var_t v1, const var_t v2) noexcept;

    /**
     * @brief get out-degree of vertex (i.e. degree of its color!)
     * 
     * @param v vertex (color)
     * @return out-degree of vert (color)
     */
    inline var_t get_out_degree(const var_t v) const noexcept { return AL_out[ v ].size(); };

    /**
     * @brief Get range over out-neighbours of v 
     * 
     * @param v vertex v
     * @return range over out-neighbors
     */
    inline auto get_out_neighbour_range(const var_t v) const noexcept {
      return AL_out[v];
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
    inline var_t get_in_degree(const var_t v) const noexcept { return AL_out[ SIGMA(v) ].size(); };

    /**
     * @brief Get range over in-neighbours of v 
     * 
     * @param v vertex v
     * @return range over in-neighbors
     */
    inline auto get_in_neighbour_range(const var_t v) const noexcept {
      return AL_out[SIGMA(v)] | std::views::transform(
                                 [this](const auto& w)
                                 { return SIGMA(w); }
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

    inline bool operator ==(const graph_al &other) const noexcept { return to_str() == other.to_str(); };

    graph_al& operator=(graph_al& g) noexcept {
      no_v = g.no_v;
      no_e = g.no_e;
      L = g.L;
      IL = g.IL;
      AL_out = g.AL_out;

      return *this;
    }
};
