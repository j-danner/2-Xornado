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

#include <iostream>
#include <algorithm>
#include <iterator>
#include <parallel/algorithm>
#include <functional>

#include <list>

#include "../misc.hpp"

#include "lineral.hpp"
#include "lineqs.hpp"
#include "omp.h"

//implementation inspired by the one of 3BA by Jan Horacek

#define DIFF diff_[omp_get_thread_num()]

// this suppress creating the new objects again and again
// (each thread has their own diff-vec)
vec< vec<var_t> > diff_( omp_get_max_threads() );


size_t lineral::hash() const {
    size_t h = idxs.size() + (p1 ? 1 : 0);
    h = p1 ? h : h^~0;
    for (auto &&i : idxs) {  
        h = (h << i) ^ ~i;
    }
    return h;
};

//gcc-dependent integer log2 func
#define LOG2(X) ((int) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

bool lineral::reduce(const LinEqs& sys) {
    bool changed = false;
    if( size() > LOG2(size())*sys.size() ) {
        //complexity to find correct update linerals: O( log( this.size() ) * sys.size() )
        for (const auto &lt_row_idx : sys.get_pivot_poly_idx()) {
            const var_t lt      = lt_row_idx.first;
            const var_t row_idx = lt_row_idx.second;
            if( (*this)[lt] ) {
                *this += sys.get_linerals( row_idx );
                changed = true;
            }
        }
    } else {
        //complexity to find correct update linerals: amortized O( this.size() )
        auto upd_idxs = std::list<var_t>();
        const auto& pivot_poly_idx = sys.get_pivot_poly_idx();
        for(const auto& l : idxs) {
            auto search = pivot_poly_idx.find(l);
            if( search != pivot_poly_idx.end() ) upd_idxs.push_back( search->second );
        }
        for(const auto& row_idx: upd_idxs) *this += sys.get_linerals( row_idx );
        changed = !upd_idxs.empty();
    }
    return changed;
};

bool lineral::reduce(const vec<lineral>& assignments) {
    bool ret = false;
    var_t offset = 0;
    while(offset<idxs.size()) {
        if( assignments[ idxs[offset] ].LT()>0 ) {
            ret = true;
            *this += assignments[ idxs[offset] ];
        } else {
            ++offset;
        }
    }
    return ret;
};

bool lineral::reduce(const vec<lineral>& assignments, const vec<var_t>& assignments_dl, const var_t lvl) {
    bool ret = false;
    var_t offset = 0;
    while(offset<idxs.size()) {
        if( assignments[ idxs[offset] ].LT()>0 && assignments_dl[ idxs[offset] ] <= lvl ) {
            ret = true;
            *this += assignments[ idxs[offset] ];
        } else {
            ++offset;
        }
    }
    return ret;
};


vec<var_t> lineral::reducers(const vec<lineral>& assignments) const {
    vec<var_t> ret;
    lineral l(*this);
    for(var_t offset = 0; offset<l.idxs.size(); offset++) {
        if( assignments[ l.idxs[offset] ].LT()>0 ) {
            ret.emplace_back( l.idxs[offset] );
            l += assignments[ l.idxs[offset] ];
            offset--;
        }
    }
    return ret;
};

bool lineral::lt_reduce(const vec<lineral>& assignments) {
    const bool ret = !assignments[LT()].is_zero();
    while( !assignments[LT()].is_zero() ) (*this)+=assignments[LT()];
    return ret;
};

std::string lineral::to_str() const {
    //if empty
    if(idxs.size() == 0 && !has_constant()) return "0";
    //else construct string
    std::string str;
    for (var_t i = 0; i < idxs.size(); i++)
    {
        str.append("x"+std::to_string( idxs[i] )+"+");
    }
    if(has_constant()) {
        str.append("1");
    } else {
        str.pop_back();
    }
    return str;
};

std::string lineral::to_xnf_str() const {
    //if empty
    if(idxs.size() == 0 && !has_constant()) return "";
    //else construct string
    std::string str;
    if(!has_constant()) {
        str.append("-");
    }
    for (var_t i = 0; i < idxs.size(); i++) {
        str.append( std::to_string( idxs[i] )+"+" );
    }
    if(idxs.size()>0) str.pop_back();
    return str;
};

std::string lineral::to_full_str(var_t num_vars) const{ 
    std::string str(num_vars, '0');
    for (auto &&i : idxs) {
        str[i] = '1';
    }
    if(has_constant()) str[0]='1';
    std::rotate(str.begin(), str.begin()+1, str.end());

    return str;
};

//overloaded operators
lineral lineral::operator+(const lineral &other) const {
    /* \warning we assume that both linerals have same num_vars (!) */
    DIFF.clear(); // DIFF is declared global and static, this saves creating new DIFFs for each calling
    std::set_symmetric_difference(std::execution::par, idxs.begin(), idxs.end(), other.idxs.begin(), other.idxs.end(), std::back_inserter(DIFF));
    //NOTE back_insterter might lead to repeated reallocations!
    //idxs = DIFF;

    return lineral(DIFF, p1^other.p1, true); //call ctor that does NOT sort DIFF
};

//in-place operation (!)
lineral& lineral::operator +=(const lineral& other) {
    if(other.size()==0) { p1^=other.p1; return *this; }

    DIFF.clear(); // DIFF is declared global and static, this saves creating new DIFFs for each calling
    std::set_symmetric_difference(std::execution::par, idxs.begin(), idxs.end(), other.idxs.begin(), other.idxs.end(), std::back_inserter(DIFF));
    std::swap(idxs, DIFF);

    p1 ^= other.p1;

    return *this;
};


bool lineral::operator <(const lineral& other) const {
    //get min of sizes
    var_t m = idxs.size() > other.idxs.size() ? other.idxs.size() : idxs.size();
    for (var_t idx = 0; idx < m; ++idx) {
        if(idxs[idx] > other.idxs[idx]) return false;
    }
    return true;
};

std::ostream& lineral::operator<<(std::ostream& os) const {
    os << to_str();
    return os;
};
