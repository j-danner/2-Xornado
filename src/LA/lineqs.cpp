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

#include "lineqs.hpp"

#include <unordered_map>
#include <map>
#include <sstream>



void LinEqs::rref() {
    pivot_poly_idx.clear();
    for (var_t i = 0; i < linerals.size(); i++) {
        //reduce new row
        for (const auto &[lt,row_idx] : pivot_poly_idx) {
            if(linerals[i][lt]) {
                linerals[i] += linerals[ row_idx ];
            }
        }
        if(!linerals[i].is_zero() ) {
            //if non-zero, add to LT_to_row_idx-map
            const var_t new_lt = linerals[i].LT();
            //add new LT to map
            pivot_poly_idx[ new_lt ] = i;

            if (new_lt == 0) continue; //no full-reduction for constant!
            //full-reduction of previous pivot-rows, i.e., reduce all previously found rows:
            for (const auto &lt_row_idx : pivot_poly_idx)
            {
                const int r_idx = lt_row_idx.second;
                if( r_idx != i && (linerals[r_idx])[new_lt] ) linerals[r_idx] += linerals[i];
            }
        } else {
            //if zero, remove row from linerals + adjust running var i
            linerals.erase( linerals.begin()+i );
            i--;
        }
    }
};

lineral LinEqs::reduce(const lineral& l) const {
    //TODO optimize by reducing given row -- without need to create copy of row!
    lineral l_(l);
    for (const auto &lt_row_idx : pivot_poly_idx) {
        const var_t lt = lt_row_idx.first;
        if(l_[lt]) {
            const int row_idx = lt_row_idx.second;
            l_ += linerals[ row_idx ];
        }
    }
    return l_;
}


void LinEqs::lt_update(const lineral& l) {
    const auto search = pivot_poly_idx.find( l.LT() );
    if(search != pivot_poly_idx.end()) {
        const auto i = search->second;
        //LT found -- start reduction!
        linerals[i] += l;
        pivot_poly_idx.erase( search );
        //rm from pivot_poly_idx, then reduce with other eqs
        for (const auto &lt_row_idx : pivot_poly_idx) {
            const var_t lt = lt_row_idx.first;
            if(linerals[i][lt]) {
                const int row_idx = lt_row_idx.second;
                linerals[i] += linerals[ row_idx ];
            }
        }
        //if non-zero, add back to pivot_poly_idx
        if(!linerals[i].is_zero()) {
            pivot_poly_idx[linerals[i].LT()] = i;
        }
    }
};

void LinEqs::lt_update(const vec<lineral>& assignments) {
    for(auto& l : linerals) l.reduce(assignments);

    pivot_poly_idx.clear();
    for (var_t i = 0; i < linerals.size(); i++) {
        //reduce new row
        for (const auto &[lt,row_idx] : pivot_poly_idx) {
            //if(linerals[i][lt]) {
            if(linerals[i].LT() == lt) {
                linerals[i] += linerals[ row_idx ];
            }
        }
        if(!linerals[i].is_zero() ) {
            //if non-zero, add to LT_to_row_idx-map
            const var_t new_lt = linerals[i].LT();
            //add new LT to map
            pivot_poly_idx[ new_lt ] = i;
        }
    }
    return;

    for (var_t i = 0; i < linerals.size(); i++) {
        //reduce new row
        if( linerals[i].reduce(assignments) ) {
            for (const auto &[lt,row_idx] : pivot_poly_idx) {
                if(linerals[i][lt]) {
                    linerals[i] += linerals[ row_idx ];
                }
            }
            if(!linerals[i].is_zero() ) {
                //if non-zero, add to LT_to_row_idx-map
                const var_t new_lt = linerals[i].LT();
                //add new LT to map
                pivot_poly_idx[ new_lt ] = i;

                if (new_lt == 0) continue; //no full-reduction for constant!
                //full-reduction of previous pivot-rows, i.e., reduce all previously found rows:
                while( pivot_poly_idx.contains( linerals[i].LT() ) ) linerals[i] += linerals[ pivot_poly_idx[linerals[i].LT()] ];
            } else {
                //if zero, remove row from linerals + adjust running var i
                //linerals.erase( linerals.begin()+i );
                //i--;
            }
        }
    }

};


void LinEqs::lt_update(const vec<lineral>& assignments, const vec<var_t>& assignments_dl, const var_t dl) {
    for(auto& l : linerals) l.reduce(assignments, assignments_dl, dl);
    rref();
    return;

    pivot_map<var_t,var_t> ppi_cpy = pivot_poly_idx;
    for(const auto& [lt,r_idx] : ppi_cpy) {
        if(assignments[lt].is_zero() || assignments_dl[lt]>dl) continue;
        //rm from pivot_poly_idx, then reduce with assignments
        pivot_poly_idx.erase( lt );
        //reduce with assignments as long as possible!
        while( !assignments[ linerals[r_idx].LT() ].is_zero() && assignments_dl[linerals[r_idx].LT()] <= dl) {
            // (1) reduce with assignments
            while(!assignments[ linerals[r_idx].LT() ].is_zero() && assignments_dl[linerals[r_idx].LT()] <= dl) {
                linerals[r_idx] += assignments[linerals[r_idx].LT()];
            }
            // (2) reduce with other eqs in linerals -- if they have same LT
            const auto search = pivot_poly_idx.find( linerals[r_idx].LT() );
            if(search != pivot_poly_idx.end()){ linerals[r_idx] += linerals[search->second]; }
        }
        //if non-zero, add back to pivot_poly_idx
        if(!linerals[r_idx].is_zero()) {
            assert( !pivot_poly_idx.contains(linerals[r_idx].LT()) );
            pivot_poly_idx[linerals[r_idx].LT()] = r_idx;
        }
    }
};

void LinEqs::update(const vec<lineral>& assignments, const vec<var_t>& assignments_dl, const var_t dl) {
    for(auto& l : linerals) {
        l.reduce(assignments, assignments_dl, dl);
    }
    pivot_poly_idx.clear();
    rref();
};


LinEqs LinEqs::operator+(const LinEqs &other) const {
    LinEqs cpy(*this);
    cpy += other;
    return cpy;
};

LinEqs& LinEqs::operator +=(const LinEqs& other) {
    const auto orig_xlits_size = linerals.size();
    linerals.reserve( linerals.size() + other.linerals.size() );
    std::copy(other.linerals.begin(), other.linerals.end(), std::back_inserter(linerals));

    for (var_t i = orig_xlits_size; i < linerals.size(); i++) {
        //reduce new row
        for (const auto &[lt,row_idx] : pivot_poly_idx)
        {
            if(linerals[i][lt]) {
                linerals[i] += linerals[ row_idx ];
            }
        }
        if(!linerals[i].is_zero() ) {
            //if non-zero, add to LT_to_row_idx-map
            const var_t new_lt = linerals[i].LT();
            //add new LT to map
            pivot_poly_idx[ new_lt ] = i;

            //full-reduction of previous pivot-rows, i.e., reduce all previously found rows:
            for (const auto &lt_row_idx : pivot_poly_idx)
            {
                const int r_idx = lt_row_idx.second;
                if( r_idx != i && (linerals[r_idx])[new_lt] ) linerals[r_idx] += linerals[i];
            }
        } else {
            //if zero, remove row from linerals + adjust running var i
            linerals.erase( linerals.begin()+i );
            i--;
        }
    }
    return *this;
};


bool LinEqs::eval(const vec<bool>& sol) const {
    return std::all_of(linerals.begin(), linerals.end(), [&sol](lineral l) { return l.eval(sol); } );
};

void LinEqs::solve(vec<bool>& sol_) const {
    if(linerals.size()==0) return;
    //TODO can be done in parallel!!
    for (const auto &lt_row_idx : pivot_poly_idx) {
        const var_t lt = lt_row_idx.first;
        const int row_idx = lt_row_idx.second;
        //update sol_[lt]: if sol_ is a zero of the lineral do nothing, otherwise flip it
        sol_[lt-1] = linerals[row_idx].eval(sol_) ? sol_[lt-1] : !sol_[lt-1];
        //sol_[lt] = !sol_[lt] ^ linerals[row_idx].eval(sol_);
    };
};

std::string LinEqs::to_str() const {
    vec< std::string > str_xlits( linerals.size() );
    auto to_str = [](const lineral l) -> std::string {return l.to_str();};
    std::transform(std::execution::par, linerals.begin(), linerals.end(), str_xlits.begin(), to_str);
    std::sort(std::execution::par, str_xlits.begin(), str_xlits.end());
    //rotate if 1 is first element, i.e., if LinEqs is inconsistent!
    if(!is_consistent()) std::rotate(str_xlits.begin(), str_xlits.begin()+1, str_xlits.end());
    std::stringstream ss;
    std::copy(str_xlits.begin(), str_xlits.end(), std::ostream_iterator<std::string>(ss, " "));
    std::string result = ss.str();
    if (!result.empty()) {
        result.resize(result.length() - 1); // trim trailing space
        return result;
    } else {
        return "0";
    }
};


//Intersect xsyses

vec<lineral> intersect(const LinEqs& U, const LinEqs& W) {
    //Zassenhaus Alg: Put U, W in Matrix [ U U \\ W 0 ] and compute rref [ A 0 \\ 0 B ]. Then B corr to basis of U \cap W.
    //NOTE assumes that n_vars is less than half of max val (of var_t type)

    //if U contains 1 return W and vice versa
    if(!U.is_consistent()) return W.get_linerals();
    if(!W.is_consistent()) return U.get_linerals();
    
    //rewrite linerals s.t. they have a continous range of idxs
    vec<var_t> supp = vec<var_t>({0});
    vec<var_t> tmp;
    for(const auto& l : U.get_linerals()) {
        std::set_union( supp.begin(), supp.end(), l.begin(), l.end(), std::back_insert_iterator(tmp) );
        std::swap(tmp, supp);
        tmp.clear();
    }
    for(const auto& l : W.get_linerals()) {
        std::set_union( supp.begin(), supp.end(), l.begin(), l.end(), std::back_insert_iterator(tmp) );
        std::swap(tmp, supp);
        tmp.clear();
    }
    //now supp contains all lits that occur in U and W
    std::unordered_map<var_t, var_t> Isupp;
    for(var_t i=0; i<supp.size(); ++i) Isupp[ supp[i] ] = i;
    const var_t n_vars = supp.size();

    rci_t nrows = U.dim() + W.dim();
    rci_t ncols = 2*(n_vars+1);

    mzd_t* M = mzd_init(nrows, ncols);
    assert( mzd_is_zero(M) );

    //fill with U
    rci_t r = 0;
    for(const auto& l : U.get_linerals()) {
        if(l.is_zero()) continue;
        if(l.has_constant()) {
            mzd_write_bit(M, r, 0, 1);
            mzd_write_bit(M, r, n_vars+1, 1);
        }
        for(const auto& i : l.get_idxs_()) {
            assert(i>0);
            assert(Isupp[i]+n_vars+1<ncols);
            mzd_write_bit(M, r, Isupp[i], 1);
            mzd_write_bit(M, r, Isupp[i]+n_vars+1, 1);
        }
        ++r;
    }
    //fill with W
    for(const auto& l : W.get_linerals()) {
        if(l.is_zero()) continue;
        if(l.has_constant()) mzd_write_bit(M, r, 0, 1);
        for(const auto& i : l.get_idxs_()) {
            assert(i>0);
            mzd_write_bit(M, r, Isupp[i], 1);
        }
        ++r;
    }
    assert(r == nrows);

    //compute rref
    const rci_t rank = mzd_echelonize(M, true);
    //read results
    vec<lineral> int_lits;
    r = rank-1;
    while(r>0) {
        mzd_t* r_ = mzd_init_window(M, r, 0, r+1, n_vars+1);
        if(!mzd_is_zero(r_)) { mzd_free_window(r_); break;}
        mzd_free_window(r_);

        vec<var_t> idxs;
        for(rci_t c=(n_vars+1)+1; c<ncols; ++c) {
            if( mzd_read_bit(M, r, c) ) idxs.push_back(supp[c-n_vars-1]);
        }
        int_lits.emplace_back( std::move(idxs), (bool) mzd_read_bit(M, r, n_vars+1), true );
        --r;
    }

    mzd_free(M);
    
    return int_lits;
}

std::pair<bool, lineral> intersectaffineVS(const LinEqs& U, const LinEqs& W) {
    //replace U and W with matrix representations; find l s.t. l in U and l+1 in W, i.e.,
    // if there is [a b] s.t. [U^T W^T] * [a \\ b] = 1
    //we compute M = [U^T W^T] and solve with righthandside 1

    if(!W.is_consistent()) return {true, lineral()};
    if(!U.is_consistent()) return {true, lineral()}; //TODO not entirely shure about this...
    
    //rewrite linerals s.t. they have a continous range of idxs
    vec<var_t> supp = vec<var_t>({0});
    vec<var_t> tmp;
    for(const auto& l : U.get_linerals()) {
        std::set_union( supp.begin(), supp.end(), l.begin(), l.end(), std::back_insert_iterator(tmp) );
        std::swap(tmp, supp);
        tmp.clear();
    }
    for(const auto& l : W.get_linerals()) {
        std::set_union( supp.begin(), supp.end(), l.begin(), l.end(), std::back_insert_iterator(tmp) );
        std::swap(tmp, supp);
        tmp.clear();
    }
    //now supp contains all lits that occur in U and W
    std::unordered_map<var_t, var_t> Isupp;
    for(var_t i=0; i<supp.size(); ++i) Isupp[ supp[i] ] = i;
    const var_t n_vars = supp.size();

    rci_t ncols = U.size() + W.size();
    rci_t nrows = n_vars;

    mzd_t* M = mzd_init(nrows, ncols);
    assert( mzd_is_zero(M) );

    //fill with U^T
    rci_t r = 0;
    for(const auto& l : U.get_linerals()) {
        if(l.has_constant()) {
            mzd_write_bit(M, 0, r, 1);
        }
        for(const auto& i : l.get_idxs_()) {
            assert(i>0);
            mzd_write_bit(M, Isupp[i], r, 1);
        }
        ++r;
    }
    //fill with W^T
    for(const auto& l : W.get_linerals()) {
        if(l.has_constant()) mzd_write_bit(M, 0, r, 1);
        for(const auto& i : l.get_idxs_()) {
            assert(i>0);
            mzd_write_bit(M, Isupp[i], r, 1);
        }
        ++r;
    }
    assert(r == ncols);

    //solve for M x = b for b =1
    mzd_t* b = mzd_init(std::max(ncols,nrows), 1);
    mzd_write_bit(b, 0, 0, 1); //uses that supp[0]==0

    //find solution
    //mzd_print(M);
    //std::cout << std::endl;
    //mzd_print(b);
    //std::cout << std::endl;
    const auto ret = mzd_solve_left(M, b, 0, true);
    if(ret==-1) return {false, lineral()};
    assert(ret == 0);
    //mzd_print(b);

    //take first half of one solution

    //comp lineral l in U with l+1 in W
    lineral out = lineral();
    for(rci_t r = 0; r<U.size(); ++r) {
        if(mzd_read_bit(b, r, 0)) out += U.get_linerals(r);
    }

    mzd_free(b);
    mzd_free(M);

    assert( U.reduce(out).is_zero() );
    assert( W.reduce(out).is_one() || !W.is_consistent());

    return {true,out};
}

vec<lineral> extend_basis(const vec<lineral>& B, const LinEqs& L) {
    LinEqs b(B);
    vec<lineral> out;
    for(const lineral& l : L.get_linerals()) {
        if(b.dim() == L.dim()) break; //correct dimension reached!
        lineral l_ = b.reduce(l);
        if(!l_.is_zero()) {
            //should we add b.reduce(l) instead of l (?)
            //out.push_back(l);
            out.push_back(l_);
            b += LinEqs(l_);
        }
    }
    assert(b.dim() == L.dim());
    return out;
}

//xcls sres_opt(xcls& cl1, xcls& cl2) {
//    //rewrite cl1 and cl2 s.t. we can do s-resolution
//    //compute intersection of cl1.assVS and cl2.assVS
//
//    if(cl1.is_zero()) return cl2;
//    if(cl2.is_zero()) return cl1;
//
//    const auto& [b,a] = intersectaffineVS( cl1.get_ass_VS(), cl2.get_ass_VS() );
//    if(!b) {
//        //no non-trivial reduction can be made --> concat ass_VS!
//        LinEqs new_ass_xsys = cl1.get_ass_VS() + cl2.get_ass_VS();
//        return xcls( new_ass_xsys );
//    }
//    //rewrite clauses to increase s-resolution!
//    //compute intersection:
//    vec<lineral> intVS = intersectVS( cl1.get_ass_VS(), cl2.get_ass_VS() );
//
//    //compute partial bases and their extensions
//    vec<lineral> &pB = intVS; //partial basis of cl1
//    pB.push_back(a);
//    //std::copy(intVS.begin(), intVS.end(), std::back_inserter(pB));
//    //extend bases
//    vec<lineral> cl1_ext = extend_basis(pB, cl1.get_ass_VS());
//    pB.back().add_one();
//    vec<lineral> cl2_ext = extend_basis(pB, cl2.get_ass_VS());
//
//    //compute s-resulution: if intVS = {l1, ..., lk} and cl1_ext = {f1,...,fn} and cl2_ext = {g1,...,gm}, then
//    //                      sres(cl1,cl2) = sres({ a, a+l1,...,a+lk, f1,...,fn },{ a+1, a+l1+1,...,a+lk+1, g1,...,gm })
//    //                                    = { a+(a+l1+1)+1,..., a+l(k-1)+(a+lk+1)+1, f1,...,fn, g1,...,gm }
//    //                                    = { l1,...,lk, f1,...,fn, g1,...,gm }
//    // more details:
//    // V_cl1 = <a, a+l1,...,a+lk, f1,...,fn> thus       cl1 = {a+1, a+l1+1,...,a+lk+1, f1+1,...,fn+1}
//    // V_cl2 = <a+1, a+1+l1,...,a+1+lk, g1,...,gm> thus cl2 = {a, a+l1,...,a+lk, g1+1,...,gm+1}
//    // S = sres(cl1,cl2) = { (a+1)+(a+l1),...,(a+1)+(a+lk), f1+1,...,fn+1, g1+1,...,gm+1 }
//    //                   = { l1+1,...,lk+1, f1+1,...,fn+1, g1+1,...,gm+1 }
//    //   V_S = < l1,...,lk, f1,...,fn, g1,...,gm >
//
//    vec<lineral> lits;
//    std::move(cl1_ext.begin(), cl1_ext.end(), std::back_inserter(lits));
//    std::move(cl2_ext.begin(), cl2_ext.end(), std::back_inserter(lits));
//    std::move(intVS.begin(), intVS.end()-1, std::back_inserter(lits));
//
//    return xcls( std::move(LinEqs(std::move(lits))) );
//};
