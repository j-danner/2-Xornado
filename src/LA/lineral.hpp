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

//from std
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <execution>
#include <memory>

#include "../misc.hpp"
//#include "LinEqs.hpp"
//forward declaration of class LinEqs
class LinEqs;

enum class cnst { zero, one };

//sparse implementation of a xor-literal
class lineral
{
    private:
        bool p1;
        //sparse repr of literal
        vec< var_t > idxs; /**<  List of sorted indices of the terms. */

    public:
        lineral() noexcept : p1(false), idxs(vec<var_t>({})) {};
        explicit lineral(const cnst zero_one) noexcept : p1(zero_one == cnst::one), idxs(vec<var_t>({})) {};
        lineral(lineral&& l) noexcept : p1(std::move(l.p1)), idxs(std::move(l.idxs)) {};
        lineral(const lineral& l) noexcept : p1(l.p1), idxs(l.idxs) {}; // no init required, as l.idxs is already sorted (i.e. initialized!)
        //b can be set to true if idxs_ is already sorted...
        lineral(const vec< var_t >& idxs_, const bool b = false) noexcept : p1(false), idxs(std::move(idxs_)) {
          if(!b){ init(); }
          else if( idxs.size()>0 && idxs[0]==0 ) { idxs.erase(idxs.begin()); p1^=true; }
        };
        lineral(vec< var_t >&& idxs_, const bool b = false) noexcept : p1(false), idxs(std::move(idxs_)) {
          if(!b){ init(); }
          else if( idxs.size()>0 && idxs[0]==0 ) { idxs.erase(idxs.begin()); p1^=true; }
        };
        lineral(const vec< var_t >& idxs_, const bool p1_, const bool b) noexcept : p1(p1_), idxs(idxs_) { if(!b){ init(); } };
        lineral(vec< var_t >&& idxs_, const bool p1_, const bool b) noexcept : p1(p1_), idxs(std::move(idxs_)) { if(!b){ init(); } };

        ~lineral() = default;

        inline void init() noexcept { 
            //sort
            std::sort(std::execution::par, idxs.begin(), idxs.end()); 
            //remove duplicates //TODO do we need this?
            idxs.erase( std::unique( idxs.begin(), idxs.end() ), idxs.end() );
            if( idxs.size()>0 && idxs[0]==0 ) { idxs.erase(idxs.begin()); p1^=true; }
            assert( idxs.empty() || idxs[0]!=0);
        }

        inline void reset() { p1=false; idxs.clear(); assert(is_zero()); };

        inline bool is_one() const { return p1 && (idxs.empty()); };
        inline bool is_zero() const { return !p1 && idxs.empty(); };

        inline bool has_constant() const { return p1; };

        inline var_t LT() const { return idxs.empty() ? 0 : idxs[0]; };

        size_t hash() const;

        inline lineral plus_one() const { return lineral( idxs, !p1, true ); };

        inline lineral add_one() { p1 ^= true; return *this; };

        bool reduce(const LinEqs& sys);
        bool reduce(const vec<lineral>& assignments, const vec<var_t>& assignments_dl, const var_t lvl);
        
        bool reduce(const vec<lineral>& assignments);
        vec<var_t> reducers(const vec<lineral>& assignments) const;
        
        /**
         * @brief lt-reduce lineral with assignments, i.e., reduce with assignments as long as LT can be reduced
         * 
         * @param assignments assignments to reduce with
         * @return true iff lineral was changed
         */
        bool lt_reduce(const vec<lineral>& assignments);

        inline vec<var_t> get_idxs() const { vec<var_t> r = idxs; if(p1){ r.insert(r.begin(), 0); } return r; };
        inline const vec<var_t>& get_idxs_() const { return idxs; };

        typedef vec<var_t>::const_iterator iterator;
        iterator begin() const { return idxs.begin(); };
        iterator end() const { return idxs.end(); };

        inline int size() const { return idxs.size(); };

        std::string to_str() const;
        std::string to_xnf_str() const;
        std::string to_full_str(var_t num_vars) const;

        //overloaded operators
	      lineral operator+(const lineral &other) const;
        //in-place operation (!)
        lineral& operator +=(const lineral& other);	
        inline lineral& operator =(const lineral& other) noexcept { idxs = other.idxs; p1 = other.p1; return *this; };
        inline lineral& operator =(const lineral&& other) noexcept { idxs = std::move(other.idxs); p1 = std::move(other.p1); return *this; };
        //lineral& operator =(lineral&& other) : idxs(std::move(other.idxs)) { return *this; }; //NOTE fails to compile...

        void swap(lineral& other) { std::swap(idxs, other.idxs); };

        inline bool operator ==(const lineral& other) const { return (p1==other.p1) && (idxs==other.idxs); };
        bool operator <(const lineral& other) const;
        inline bool operator[](const var_t idx) const { return idx==0 ? p1 : std::binary_search(idxs.begin(), idxs.end(), idx); };
        std::ostream& operator<<(std::ostream& os) const;

        bool eval(const vec<bool> &sol) const { bool out = !p1; for(const auto &i : idxs) out ^= sol[i-1]; return out; };
        void solve(vec<bool>& sol_) const { if(LT()>0) { sol_[LT()-1] = eval(sol_) ? sol_[LT()-1] : !sol_[LT()-1]; } };
};

namespace std {
  template <>
  struct hash<lineral> {
    std::size_t operator()(const lineral& k) const {
      return k.hash();
    }
  };
}

