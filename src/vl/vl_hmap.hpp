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

#include <memory>
#include <stack>

#include "../LA/lineral.hpp"

#include "../parallel-hashmap/phmap.h"
#include "../robin_hood-3.11.5/robin_hood.h"

//choose hmap implementation
template<class K, class V>
//using hmap = std::unordered_map<K, V>;
//using hmap = robin_hood::unordered_flat_map<K, V>;
using hmap = phmap::parallel_flat_hash_map<K, V,
                                          //phmap::priv::hash_default_hash<K>,
                                          robin_hood::hash<K>,
                                          phmap::priv::hash_default_eq<K>,
                                          phmap::priv::Allocator<std::pair<const K, V>>,
                                          4,
                                          std::mutex>;

struct vl_hmap_repr {
  var_t lvl;

  vl_hmap_repr(const var_t _lvl) noexcept : lvl(_lvl) {};
};

struct vl_hmap_insert_return_type {
    const bool inserted;
    const var_t vert;

    vl_hmap_insert_return_type(const bool& _inserted, const var_t _vert) : inserted(_inserted), vert(_vert) {};
    ~vl_hmap_insert_return_type() {};
};

class vl_hmap
{
  private:
    /**
     * @brief map from vertices to linerals
     */
    std::stack< hmap<var_t, lineral> > v_to_xl_stack;

    /**
     * @brief map from linerals to vertices
     */
    std::stack< hmap<lineral, var_t> > xl_to_v_stack;

  public:
    vl_hmap() noexcept : vl_hmap(1) {};

    vl_hmap(const var_t _num_verts) noexcept : v_to_xl_stack(std::stack< hmap<var_t,lineral> >()), xl_to_v_stack(std::stack< hmap<lineral,var_t> >()) {
        v_to_xl_stack.emplace( hmap<var_t,lineral>(_num_verts) );
        xl_to_v_stack.emplace( hmap<lineral,var_t>(_num_verts) );
    };
    
    vl_hmap(const var_t _num_verts, [[maybe_unused]] const var_t _num_vars) noexcept : vl_hmap(_num_verts)  {};

    vl_hmap(const vl_hmap& vl) noexcept : v_to_xl_stack(vl.v_to_xl_stack), xl_to_v_stack(vl.xl_to_v_stack) {};
    
    vl_hmap(vl_hmap&& vl) noexcept : v_to_xl_stack(std::move(vl.v_to_xl_stack)), xl_to_v_stack(std::move(vl.xl_to_v_stack)) {};
    
    ~vl_hmap() {};

    inline void put_Vxlit(hmap<var_t,lineral>&& _Vxlit) noexcept { v_to_xl_stack.top() = std::move(_Vxlit); };
    inline void put_V(hmap<lineral,var_t>&& _V) noexcept { xl_to_v_stack.top() = std::move(_V); };

    inline vl_hmap_repr get_state() {
      v_to_xl_stack.push( v_to_xl_stack.top() );
      xl_to_v_stack.push( xl_to_v_stack.top() );
      return vl_hmap_repr(v_to_xl_stack.size());
    };
    
    inline void backtrack(vl_hmap_repr&& r, [[maybe_unused]] const var_t dl) noexcept {
      while(xl_to_v_stack.size()>=r.lvl) {
        xl_to_v_stack.pop();
        v_to_xl_stack.pop();
      }
    };

    inline var_t size() const noexcept { return xl_to_v_stack.top().size(); };

    /**
     * @brief inserts lit if not yet present
     * 
     * @param v vertex index
     * @param lit lineral to be inserted
     * @param dl current decision level
     * @return vl_hmap_insert_return_type field inserted is true iff lit could be assigned to v; field vert points to the vertex representing lit
     */
    const vl_hmap_insert_return_type insert(const var_t v, lineral&& lit, const var_t dl);

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
    std::pair<var_t,bool> update(const var_t v, lineral&& l, const var_t dl);

    /**
     * @brief retrieves lineral of vertex
     * 
     * @param v vertex to search literal of
     * @note complexity is linear in size of literal to be found (if there is none, constant)
     * @return lineral literal of vertex v
     */
    inline lineral operator[](const var_t v) const noexcept { return v_to_xl_stack.top().at(v); };

    /**
     * @brief retrieves vertex of lit
     * 
     * @param lit literal to search vertex of
     * @note complexity is amortized linear in lit.get_idxs().size()
     * @return var_t vertex of lit; 0 if there is none!
     */
    inline var_t operator[](const lineral& lit) const noexcept { return xl_to_v_stack.top().at(lit); };

    /**
     * @brief compute string repr of object, listing all stored pairs (v, lit) in sorted order.
     * 
     * @return std::string representation of object
     */
    std::string to_str() const;

    /**
     * @brief check if given literal has an assigned vert
     * 
     * @param lit literal to check containment
     * @return true iff literal has vert in trie
     */
    inline bool contains(const lineral& lit) const noexcept { return xl_to_v_stack.top().contains(lit); };
    
    /**
     * @brief check if given vert has a corr literal
     * 
     * @param v vert to check containment
     * @return true iff vert has literal
     */
    inline bool contains(const var_t v) const noexcept { return v_to_xl_stack.top().contains(v); };


    //#define V(l) V_stack.top().at(l)
    inline var_t V(const lineral &l) const {
      if(!l.has_constant()) {
        auto search = xl_to_v_stack.top().find( l );
        assert(search != xl_to_v_stack.top().end());
        return search->second;
      } else {
        return SIGMA( xl_to_v_stack.top().find( l.plus_one() )->second );
      }
    }

    //check if V contains l
    inline bool V_contains(const lineral &l) const {
      if(!l.has_constant()) {
        return xl_to_v_stack.top().contains( l );
      } else {
        return xl_to_v_stack.top().contains( l.plus_one() );
      }
    }

    //#define v_to_xl(v) v_to_xl_stack.top().at(v)
    inline lineral Vxlit(const var_t &v) const {
      auto search = v_to_xl_stack.top().find( v );
      if(search == v_to_xl_stack.top().end()) {
        auto sigma_search = v_to_xl_stack.top().find( SIGMA(v) );
        assert(sigma_search != v_to_xl_stack.top().end());
        return (sigma_search->second).plus_one(); //TODO avoid copy!
      } else {
        return search->second;
      }
    }

    inline bool Vxlit_contains(const var_t &v) const {
      auto search = v_to_xl_stack.top().find( v );
      if(search == v_to_xl_stack.top().end()) {
        return v_to_xl_stack.top().find( SIGMA(v) ) != v_to_xl_stack.top().end(); //TODO avoid copy!
      } else {
        return true;
      }
    }

    inline var_t Vxlit_LT(const var_t &v) const {
      return v_to_xl_stack.top().at( contains(v) ? v : SIGMA(v) ).LT();
    }

    vl_hmap& operator =(vl_hmap& o) noexcept {
      v_to_xl_stack = o.v_to_xl_stack;
      xl_to_v_stack = o.xl_to_v_stack;
      return *this;
    };

    vl_hmap& operator =(const vl_hmap&& o) noexcept {
      v_to_xl_stack = std::move(o.v_to_xl_stack);
      xl_to_v_stack = std::move(o.xl_to_v_stack);
      return *this;
    };

    inline lineral sum(const var_t lhs, const var_t rhs) const noexcept { return v_to_xl_stack.top().at(lhs)+v_to_xl_stack.top().at(rhs); };
};
