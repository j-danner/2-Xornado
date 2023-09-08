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

#include <string>
#include <map>

#include "lineral.hpp"

#include "../robin_hood-3.11.5/robin_hood.h"

#include <m4ri/m4ri.h>

#ifdef NDEBUG
  template<class K, class V>
  using pivot_map = robin_hood::unordered_flat_map<K,V>;
#else

  template<class K, class V>
  //using pivot_map = std::unordered_map<K,V>;
  using pivot_map = std::map<K,V>;
#endif

class LinEqs
{
  private:
    vec< lineral > linerals;

    pivot_map<var_t, var_t> pivot_poly_idx;

    void rref();
  public:
    LinEqs() noexcept { linerals = vec<lineral>(0); };
    LinEqs(const lineral& lit) noexcept : linerals(vec<lineral>({lit})) { rref(); };
    LinEqs(lineral&& lit) noexcept : linerals(vec<lineral>({std::move(lit)})) { rref(); };
    LinEqs(const vec<lineral>& xlits_) noexcept : linerals(xlits_) { rref(); };
    LinEqs(vec<lineral>&& xlits_) noexcept : linerals(std::move(xlits_)) { rref(); };
    LinEqs(const LinEqs& o) noexcept : linerals(o.linerals), pivot_poly_idx(o.pivot_poly_idx) {};
    LinEqs(LinEqs&& o) noexcept : linerals(std::move(o.linerals)), pivot_poly_idx(std::move(o.pivot_poly_idx)) {};
    ~LinEqs() = default;

    /**
     * @brief reduces a given lineral by the linsys
     * 
     * @param l given lineral
     * @return lineral LT-reduced with the linerals in the linsys
     */
    lineral reduce(const lineral& l) const;

    /**
     * @brief updates xsyses LTs modulo l
     * 
     * @param l lit to reduce with
     */
    void lt_update(const lineral& l);
    
    /**
     * @brief updates xsyses LTs modulo l
     * 
     * @param assignments to reduce with
     */
    void lt_update(const vec<lineral>& assignments);
    
    /**
     * @brief updates xsyses LTs modulo l
     * 
     * @param assignments to reduce with
     * @param assignments_dl dl of each assignment
     * @param dl dl up to which assignments are considered
     */
    void lt_update(const vec<lineral>& assignments, const vec<var_t>& assignments_dl, const var_t dl);
    
    void update(const vec<lineral>& assignments, const vec<var_t>& assignments_dl, const var_t dl);

    inline lineral get_non_zero_el() const { assert(!pivot_poly_idx.empty()); return linerals[pivot_poly_idx.begin()->second]; };

    bool is_consistent() const { return pivot_poly_idx.find(0) == pivot_poly_idx.end(); };

    /**
     * @brief evaluates LinEqs with tuple sol
     * 
     * @param sol 
     * @return true if sol is a solution of the linsys
     * @return false if sol is not a solution of the linsys
     */
    bool eval(const vec<bool>& sol) const;

    /**
     * @brief returns a solution of the linsys, 'extends' sol_ it to a solution of this linsys, i.e., the LTs of this LinEqs are determined based on the values of sol_
     * 
     * @param sol_  (partial) solution to be extended
     * @return vec<bool> solution of linsys, if there exists one
     * @note throws an assertion-failure if linsys is inconsistent
     */
    void solve(vec<bool>& sol_) const;

    std::string to_str() const;

    inline int dim() const { return pivot_poly_idx.size(); };
    
    inline int size() const { return linerals.size(); };
    
    inline const vec<lineral>& get_linerals() const { return linerals; };
    inline const lineral& get_linerals(var_t i) const { return linerals[i]; };
    inline const pivot_map<var_t,var_t>& get_pivot_poly_idx() const { return pivot_poly_idx; };

    inline bool operator ==(const LinEqs& other) const { return to_str()==other.to_str(); };
    LinEqs& operator =(const LinEqs& other) { linerals = other.linerals; pivot_poly_idx = other.pivot_poly_idx; return *this; };
    LinEqs& operator =(LinEqs&& other) { linerals = std::move(other.linerals); pivot_poly_idx = std::move(other.pivot_poly_idx); return *this; };

    bool contains_lt(const var_t lt) const { return pivot_poly_idx.contains(lt); };
    
	  LinEqs operator+(const LinEqs &other) const;
    //in-place operation (!)
    LinEqs& operator +=(const LinEqs& other);	

    void clear() { linerals.clear(); pivot_poly_idx.clear(); };
};


vec<lineral> intersect(const LinEqs& U, const LinEqs& W);

std::pair<bool, lineral> intersectaffineVS(const LinEqs& U, const LinEqs& W);

vec<lineral> extend_basis(const vec<lineral>& B, const LinEqs& L);
