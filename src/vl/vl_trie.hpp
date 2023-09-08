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
//trie-implementation that maps linerals to vertices and vice-versa

#include <vector>
#include <stack>
#include <list>
#include <map>
#include <iterator>

#include "../misc.hpp"
#include "../LA/lineral.hpp"
#include "../LA/lineqs.hpp"

#include "../robin_hood-3.11.5/robin_hood.h"
//#include "../parallel-hashmap-1.34/phmap.h"

#define ROOT 0

#ifdef NDEBUG
  template<class K, class V>
  //using child_map = std::unordered_map<K, V>;
  //using child_map = std::map<K, V>;
  using child_map = robin_hood::unordered_flat_map<K,V>;
  //using child_map = robin_hood::unordered_map<K,V>;
  //using child_map = phmap::parallel_flat_hash_map<K, V, phmap::priv::hash_default_hash<K>, phmap::priv::hash_default_eq<K>, phmap::priv::Allocator<std::pair<const K, V>>, 4, phmap::NullMutex>;
#else 
  template<class K, class V>
  using child_map = std::map<K,V>;
#endif

/**
 * @brief type of node-counter in trie
 * @note make sure it is large enough to fit all nodes! ...NO OVERFLOW CHECK!
 */
typedef unsigned int n_t;

struct trie_insert_return_type {
    const bool inserted;
    const bool found_plus_one;
    const var_t vert;

    trie_insert_return_type(const bool& _inserted, const bool& _found_plus_one, const var_t _vert) : inserted(_inserted), found_plus_one(_found_plus_one), vert(_vert) {};
    ~trie_insert_return_type() {};
};

struct trie_repr {
  const child_map<var_t,n_t> v_node;
  const n_t num_vs = 0;

  trie_repr() noexcept {};
  trie_repr(const trie_repr& o) noexcept : v_node(o.v_node), num_vs(o.num_vs) {};
  trie_repr(trie_repr&& o) noexcept : v_node(std::move(o.v_node)), num_vs(std::move(o.num_vs)) {};
  trie_repr(const child_map<var_t,n_t>& _v_node, const n_t _num_vs) noexcept : v_node(_v_node), num_vs(_num_vs) {};
  ~trie_repr() {};
};


struct node {
  /**
   * @brief node_idx of parent
   */
  n_t parent;

  /**
   * @brief label of node
   */
  var_t label;

  /**
   * @brief children of node; maps label onto childrens' node_idx
   */
  child_map<var_t,n_t> children; 

  node(const n_t _parent, const var_t _label) noexcept : parent(_parent), label(_label) {};
  node(const node& o) noexcept : parent(o.parent), label(o.label), children(o.children) {};
  node(node&& o) noexcept : parent(o.parent), label(o.label), children(std::move(o.children)) {};
  node& operator =(const node& o) noexcept { children=o.children; parent=o.parent; label=o.label; return *this; };
  node& operator =(node&& o) noexcept { children=std::move(o.children); parent=std::move(o.parent); label=std::move(o.label); return *this; };
};


class vl_trie
{
  private:
    /**
     * @brief vector of nodes
     */
    vec< node > nodes; //TODO change to stack and change nodes to have pointers to parents and children?

    /**
     * @brief v_node[v] is node_idx for vertex v
     */
    child_map<var_t,n_t> v_node;

    /**
     * @brief map assigning node_idx their vertex - if they have one!
     */
    child_map<n_t,var_t> assigned_vert;

    /**
     * @brief number of variables of the linerals
     */
    var_t num_vars;

    /**
     * @brief number of vertices stored in trie
     */
    n_t num_vs = 0;

    /**
     * @brief unused node_idxs
     * @note filled by remove_node
     */
    std::stack<n_t> unused_node_idxs;

    /**
     * @brief nodes_dl[dl] stores a list containing all node_idxs of nodes that were added in decision level dl
     */
    std::stack< std::list<n_t> > nodes_in_dl;

    void register_node(const n_t node_idx, const var_t dl) {
      while(dl >= nodes_in_dl.size()) nodes_in_dl.emplace( std::list<n_t>() );
      nodes_in_dl.top().push_back(node_idx);
    }
    
    /**
     * @brief adds a new node to the trie
     * 
     * @param parent_idx parent of new node
     * @param label label of new node
     * @param dl current decision level
     * @return n_t node_idx of new node
     */
    inline n_t add_node(const n_t parent_idx, const var_t label, const var_t dl) noexcept {
      n_t node_idx;
      if (unused_node_idxs.empty()) {
        node_idx = nodes.size();
        nodes.emplace_back( parent_idx, label );
        register_node( node_idx, dl );
      } else {
        //get unused node_idx
        node_idx = unused_node_idxs.top();
        unused_node_idxs.pop();
        assert( node_idx < nodes.size() );
        //change its members correspondingly!
        nodes[node_idx].parent = parent_idx;
        nodes[node_idx].label = label;
        nodes[node_idx].children.clear();
        register_node( node_idx, dl );
      }
      [[maybe_unused]] const auto insert = nodes[parent_idx].children.emplace( label, node_idx );
      assert(insert.second);
      return node_idx;
    }
    
    /**
     * @brief removes a node from the trie; may disconnect parts of the trie!
     * 
     * @note may also change nothing at all!
     * @param node_idx parent of new node
     */
    inline void remove_node(const n_t node_idx) noexcept {
      //rm link from parent
      nodes[nodes[node_idx].parent].children.erase( nodes[node_idx].label );
      //rm link to all children! //should be done before using this node again!
      //nodes[node_idx].children.clear();
      //re-init parent and node_label! 
      //indicate node_idx as unused!
      unused_node_idxs.push( node_idx );
    }

    inline void assign_vert(const n_t n, const var_t v) { assigned_vert[n] = v; v_node[v] = n; num_vs++; };

    inline bool has_assigned_vert(const n_t n) const { return assigned_vert.contains(n); };
    
    /**
     * @brief prunes trie to given decision-level dl, i.e., removes all nodes added later than dl
     * 
     * @param dl decision-level for which all nodes with strictly higher lvl are removed
     */
    void prune(const var_t dl) noexcept;

    inline var_t get_vert(const n_t n) const noexcept { return assigned_vert.at(n); };

  public:
    vl_trie() noexcept : vl_trie(1) {};

    vl_trie(const var_t _num_vars) noexcept : num_vars(_num_vars) {
      unused_node_idxs = std::stack<n_t>();
      nodes_in_dl = std::stack< std::list<n_t> >();
      //add root, i.e., add new els to children, parents
      nodes.push_back( node(ROOT,ROOT) );
      register_node(ROOT, 0);
    };

    vl_trie([[maybe_unused]] const var_t num_verts, const var_t _num_vars) noexcept : vl_trie(_num_vars)  { nodes.reserve(num_verts); };

    vl_trie(const vl_trie& tr) noexcept : nodes(tr.nodes), v_node(tr.v_node), assigned_vert(tr.assigned_vert), num_vars(tr.num_vars), num_vs(tr.num_vs), unused_node_idxs(tr.unused_node_idxs) {};
    
    vl_trie(vl_trie&& tr) noexcept : nodes(std::move(tr.nodes)), v_node(std::move(tr.v_node)), assigned_vert(std::move(tr.assigned_vert)), num_vars(std::move(tr.num_vars)), num_vs(std::move(tr.num_vs)), unused_node_idxs(std::move(tr.unused_node_idxs)) {};

    ~vl_trie() {};

    trie_repr get_state() const { return std::move(trie_repr(v_node, num_vs)); };

    void backtrack(trie_repr&& r, const var_t dl) noexcept;

    inline var_t size() const noexcept { return num_vs; };

    inline var_t get_num_nodes() const noexcept { return nodes.size()-unused_node_idxs.size(); };

    /**
     * @brief inserts lit if not yet present
     * 
     * @param v vertex index
     * @param lit lineral to be inserted
     * @param dl current decision level - defaults to 0
     * @return trie_insert_return_type field inserted is true iff lit could be assigned to v; field node_idx points to the node representing lit, to get corr vert use get_vert
     */
    const trie_insert_return_type insert(const var_t v, const lineral& lit, const var_t dl);

    /**
     * @brief erase vertex v from trie (along with its label)
     * 
     * @param v vertex to be removed
     * @return true iff vert could be erased
     */
    bool erase(const var_t v);

    /**
     * @brief updates lit assigned to v to l
     * 
     * @param v vertex to update
     * @param l label to change to
     * @param dl current decision level - if not provided defaults to 0
     * @return var_t vert where lit is stored (or lit+1); bool true iff vert points to lit+1
     */
    std::pair<var_t,bool> update(const var_t v, const lineral& l, const var_t dl);

    /**
     * @brief retrieves lineral of vertex
     * 
     * @param v vertex to search literal of
     * @note complexity is linear in size of literal to be found (if there is none, constant)
     * @return lineral literal of vertex v
     */
    lineral operator[](const var_t v) const;
    
    /**
     * @brief retrieves lineral of vertex
     * 
     * @param v vertex to search literal of
     * @note complexity is linear in size of literal to be found (if there is none, constant)
     * @return lineral literal of vertex v
     */
    lineral at(const var_t v) const;

    /**
     * @brief retrieves vertex of lit
     * 
     * @param lit literal to search vertex of
     * @note complexity is amortized linear in lit.get_idxs().size()
     * @return var_t vertex of lit; 0 if there is none!
     */
    var_t operator[](const lineral& lit) const;

    /**
     * @brief retrieves vertex of lit
     * 
     * @param lit literal to search vertex of
     * @note complexity is amortized linear in lit.get_idxs().size()
     * @return var_t vertex of lit;
     */
    var_t at(const lineral& lit) const;
    
    /**
     * @brief retrieves vertex of lit -- checks
     * 
     * @param lit literal to search vertex of
     * @note complexity is amortized linear in lit.get_idxs().size()
     * @return var_t vertex of lit (or lit+1); bool true iff lit was found
     */
    std::pair<var_t,bool> at_(const lineral& lit) const;

    /**
     * @brief compute string repr of object, listing all stored pairs (v, lit) in sorted order.
     * 
     * @return std::string representation of object
     */
    std::string to_str() const;

    /**
     * @brief check if given literal has an assigned vert in trie
     * 
     * @param lit literal to check containment
     * @return true iff literal has vert in trie
     */
    bool contains(const lineral& lit) const;

    
    inline var_t V(const lineral &l) const {
      const auto [v,b] = at_(l);
      //b is true if l was found at v; otherwise l+1 was found at v.
      return b ? v : SIGMA(v); 
    }

    //check if V contains l or l+1
    inline bool V_contains(const lineral &l) const {
      return contains(l) || contains(l.plus_one());
    }

    inline bool Vxlit_contains(const var_t &v) const {
      return contains(v);
    }

    inline lineral Vxlit(const var_t &v) const {
      assert( contains(v) || contains(SIGMA(v)) );
      if(contains(v)) {
        return std::move( at(v) );
      } else {
        lineral lit = at( SIGMA(v) );
        lit.add_one();
        return lit;
      }
    }

    inline var_t Vxlit_LT(const var_t &v) const {
      auto it = begin( contains(v) ? v : SIGMA(v) );
      return *it!=0 ? *it : *(++it);
    }

    
    /**
     * @brief returns a tuple allowing to find the vertex representing zero - if it exists
     * 
     * @return std::tuple<bool,bool,var_t> 1st bool true if zero exists; 2nd bool true iff 1 was found; 3rd var_t vertex found
     */
    std::tuple<bool,bool,var_t> if_exists_get_zero_v() const {
      if(has_assigned_vert(ROOT)) {
        return {true, false, get_vert(ROOT)};
      } else if(nodes[ROOT].children.contains(0) && has_assigned_vert(nodes[ROOT].children.at(0))) {
        //there is one!
        return {true, true, get_vert( nodes[ROOT].children.at(0) )};
      } else {
        //we have neither 0 nor 1...
        return {false, false, -1};
      }
    }

    /**
     * @brief check if given vert has a corr literal in trie
     * 
     * @param v vert to check containment
     * @return true iff vert has literal in trie
     */
    inline bool contains(const var_t v) const { return v_node.contains(v); };

    vl_trie& operator =(const vl_trie& o) noexcept {
      v_node = child_map<var_t,n_t>(o.v_node);
      nodes = vec< node >(o.nodes);
      assigned_vert = child_map<n_t,var_t>(o.assigned_vert);
      num_vars = o.num_vars;
      num_vs = o.num_vs;
      return *this;
    };

    vl_trie& operator =(const vl_trie&& o) noexcept {
      v_node = child_map<var_t,n_t>(std::move(o.v_node));
      nodes = vec< node >(std::move(o.nodes));
      assigned_vert = child_map<n_t,var_t>(std::move(o.assigned_vert));
      num_vars = std::move(o.num_vars);
      num_vs = std::move(o.num_vs);
      return *this;
    };


    /**
     * @brief const forward iterator, can be used instead of operator[](var_t) to avoid unnecessary copy and creation of lineral (!)
     */
    class const_iterator {
      private:
        const vl_trie* t;
        n_t curr_node;
      
      public:
        typedef std::forward_iterator_tag iterator_category;  //usually std::forward_iterator_tag or similar
        typedef var_t value_type; //almost always T
        typedef std::ptrdiff_t difference_type; //almost always ptrdiff_t
        typedef const var_t& reference; //almost always T& or const T&
        typedef const var_t* pointer; //almost always T* or const T*

        const_iterator(const vl_trie* const _t, const n_t _start_node) : t(_t), curr_node(_start_node) {};
        const_iterator(const const_iterator& o) : t(o.t), curr_node(o.curr_node) {};
        
        const node& get_node() const { return t->nodes[curr_node]; };
        var_t get_node_idx() const { return curr_node; };

        const_iterator& operator=(const const_iterator& o) { t=o.t; curr_node=o.curr_node; return *this; };
        const_iterator& operator++() { curr_node = t->nodes[curr_node].parent; return *this; }; //prefix increment
        reference operator*() const { return t->nodes[curr_node].label; };
        friend void swap(const_iterator& lhs, const_iterator& rhs) { std::swap(lhs.t, rhs.t); std::swap(lhs.curr_node, rhs.curr_node); };

        const_iterator operator++(int) { const_iterator tmp = *this; ++(*this); return tmp; }; //postfix increment
        pointer operator->() const { return &(t->nodes[curr_node].label); };
        friend bool operator==(const const_iterator& lhs, const const_iterator& rhs) { return lhs.curr_node==rhs.curr_node && lhs.t==rhs.t; };
        friend bool operator!=(const const_iterator& lhs, const const_iterator& rhs) { return lhs.curr_node!=rhs.curr_node || lhs.t!=rhs.t; };

    };

    /**
     * @brief const_iterator starting at node representing vertex v
     * 
     * @param v vertex
     * @return const_iterator starting at node repr vertex v
     */
    const_iterator begin(const var_t v) const { return const_iterator(this, v_node.at(v)); }

    /**
     * @brief const_iterator repr root, i.e., end of every path
     * 
     * @return const_iterator repr root
     */
    const_iterator end() const { return const_iterator(this, ROOT); }

    lineral sum(const var_t lhs, const var_t rhs) const;
};