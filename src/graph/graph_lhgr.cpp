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

#include "graph_lhgr.hpp"

#include <execution>

graph_lhgr::graph_lhgr(const vec< std::pair<var_t,var_t> >& E, const var_t no_v_) noexcept : no_v(no_v_), no_e(2*E.size()) {
    init(E,no_v_);
};

void graph_lhgr::init(const vec< std::pair<var_t,var_t> >& E, const var_t no_v_) noexcept {
    //init data structures!
    no_e = 0;
    no_v = no_v_;

    //active vertices
    L = vec<var_t>(no_v);
    IL = vec<var_t>(no_v);

    //adjacency list
    AL_out = vec< vec<var_t> >(no_v);
    IAL_in = vec< vec<var_t> >(no_v);

    //degree
    VD_out = vec<var_t>(no_v, 0);

    //color
    CAL = vec< std::list<var_t> >(no_v);
    
    VC = vec<var_t>(no_v);
    CD_out = vec<var_t>(no_v);

    //fill data structures:

    //fill AL_out
    for (const auto &e : E) {
        const var_t src  = e.first;
        const var_t dst = e.second;
        add_edge(src,dst);
    }

    //init L, IL, CAL
    for (var_t v = 0; v < no_v; ++v) {
        //fill L and IL
        L[v] = v;
        IL[v] = v;
        
        //fill CAL
        CAL[v] = std::list<var_t>(1, v);
        
        //init VC
        VC[v] = v; //every node has their own color

        //init VD_out, VD_in, CD_out, CD_in
        CD_out[v] = VD_out[v];
        no_e += VD_out[v];
    }

    assert( graph_lhgr::assert_data_structs() );
};


void graph_lhgr::add_edge(const var_t src, const var_t dst) noexcept {
    assert(src != dst);
    assert(src < no_v);
    assert(dst < no_v);
    //add edge
    AL_out[src].push_back(dst);
    ++VD_out[src];
    //set IAL_in
    const var_t i = AL_out[src].size()-1; //idx of dst in src

    //add symmetric edge sigma(dst)->sigma(src) --- if distinct!
    var_t j = i;
    if(SIGMA(dst) != src) {
        AL_out[ SIGMA(dst) ].push_back( SIGMA(src) );
        ++VD_out[ SIGMA(dst) ];
        //set IAL_in
        j = AL_out[SIGMA(dst)].size()-1; //idx of sigma(src) in AL_out[ sigma(dst) ]
    }
    
    //store in IAL_in -- resize vector if necessary!
    if(i >= IAL_in[src].size()) IAL_in[src].resize(i+1);
    IAL_in[src][i] = j;
    if(SIGMA(dst) != src) {
        //store in IAL_in -- resize vector if necessary!
        if(j >= IAL_in[SIGMA(dst)].size()) IAL_in[SIGMA(dst)].resize(j+1);
        IAL_in[SIGMA(dst)][j] = i;
    }
    assert( AL_out[ SIGMA(AL_out[src][i]) ][IAL_in[src][i]] == SIGMA(src) );
    assert( AL_out[ SIGMA(AL_out[SIGMA(dst)][j]) ][IAL_in[SIGMA(dst)][j]] == SIGMA(SIGMA(dst)) );
};

void graph_lhgr::backtrack(graph_lhgr_repr&& graph_orig) noexcept {
    //restore graph_lhgr from graph_lhgr_lhgr_repr  
    no_v = std::move(graph_orig.no_v);
    no_e = std::move(graph_orig.no_e);
    VD_out = std::move(graph_orig.VD_out);
    VC = std::move(graph_orig.VC);

    //rebuild CAL and CD_out
    for (var_t c_idx = 0; c_idx < no_v; ++c_idx) {
        //reset color lists in CAL
        CAL[L[c_idx]].clear();

        //reset CD_out
        CD_out[L[c_idx]] = 0;
    }

    //apply known colors to construct CAL and generate CD_out
    for (var_t v = 0; v < VC.size(); ++v) {
        CAL[ VC[v] ].push_back(v);
        CD_out[ VC[v] ] += VD_out[v];
    }
};

//remove single edge src->AL_out[src][idx] (does not touch its symmetric edge!)
void graph_lhgr::remove_edge_(const var_t src, const var_t idx) noexcept {
    //do nothing if edge is already out-of-scope (might be relevant when this func is called from remove_all_edges where symmetrical edges are processed)
    if(idx >= VD_out[src]) return;

    //adjust degrees:
    //decrease color-degree of colors of src and dst
    --CD_out[VC[src]];
    //decrease degree of src and dst
    --VD_out[src];
    
    const var_t dst =  AL_out[src][idx];
    //new position in AL_out_logical[src]
    const var_t idx_ = VD_out[src];
    const var_t dst_ = AL_out[src][idx_];
    
    //decrease number of edges
    --no_e;
    
    //move dst out-of-range in AL_out[src]
    std::swap( AL_out[src][idx], AL_out[src][idx_] );
    //adapt IAL_in
    std::swap( IAL_in[ SIGMA(dst_) ][ IAL_in[src][idx_] ], IAL_in[ SIGMA(dst) ][ IAL_in[src][idx] ] );
    //could do without local var dst_ -- but probably compiler optimizes that away anyways...
    //std::swap( IAL_in[ SIGMA(AL_out[src][idx]) ][ IAL_in[src][idx_] ], IAL_in[ SIGMA(dst) ][ IAL_in[src][idx] ] );
    std::swap( IAL_in[src][idx], IAL_in[src][idx_] );
}

//removes idx-th going out from src as listed in AL_out_logical(src,idx) (ignores same-colored verts!)
void graph_lhgr::remove_edge(const var_t src, const var_t idx) noexcept {
    //remove src -> AL_out[src,idx]
    remove_edge_(src, idx);
    //remove sigma(AL_out[src,idx]) -> IAL_out
    remove_edge_(SIGMA(AL_out[src][VD_out[src]]), IAL_in[src][VD_out[src]]);
};

//removes all out-going edges of src
void graph_lhgr::remove_all_edges(const var_t v) noexcept {
    //loop over verts of same color and remove all of their out-edges
    for (const auto &src : CAL[VC[v]]) {
        //for every out-going edge (src,dst) in E, remove (src,dst) and (sigma(dst)(src));
        for (var_t idx = 0; idx < VD_out[src]; ++idx) {
            //remove edge sigma(dst)->sigma(scr) -- if is NOT not self-symmetric (i.e. an out-going edge itself)
            if( SIGMA(AL_out[src][idx]) != src ) remove_edge_(SIGMA(AL_out[src][idx]), IAL_in[src][idx]);
        };
    
        //remove all out-edges of src setting the degree to 0. (symmetric edges untouched!)
        CD_out[ VC[src] ] -= VD_out[src];
        no_e -= VD_out[src];
        VD_out[src] = 0;
    };
};

//remove color/vert c and SIGMA(c)
void graph_lhgr::remove_vert(const var_t c) noexcept {
    //std::cout << "no_v = " << std::to_string(no_v) << "; removing " << std::to_string(c) << " (" << std::to_string(SIGMA(c)) << ")" << std::endl;
    //remove c and SIGMA(c)
    for (const auto &v : std::vector<var_t>({c, SIGMA(c)})) {
        --no_v;
        //swap vertices in L and fix IL, move v 'out-of-scope'
        std::swap( L[ IL[VC[v]] ], L[ no_v ] );
        //fix IL accordingly
        std::swap( IL[ L[IL[VC[v]]] ], IL[ L[no_v] ] );

        this->remove_all_edges(v);

        CAL[v].clear();
    };
};


vec<bool> exists_edge_to_c( 100, false ); //TODO instead of large array use hash_set?
//merges colors v1 and v2, color of v2 will be set to color of v1; and color of SIGMA(v2) set to color of SIGMA(v1)
void graph_lhgr::merge_verts(const var_t v1_, const var_t v2_) noexcept {
    //std::cout << "no_v = " << std::to_string(no_v) << "; merging " << std::to_string(VC[v1]) << " and " << std::to_string(VC[v2]) << " (" << std::to_string(VC[SIGMA(v1)]) << " and " << std::to_string(VC[SIGMA(v2)]) << ")" << std::endl;
    const var_t v1 = VC[v1_];
    const var_t v2 = VC[v2_];
    //abort if already merged!
    if(VC[v1]==VC[v2]) return;

    assert(VC[v1] == v1);
    assert(VC[v2] == v2);
    assert(v1 != v2);
    assert(IL[v1] < no_v );
    assert(IL[v2] < no_v );
    assert( graph_lhgr::assert_data_structs() );

    --no_v;
    //swap vertices in L to move v2 'out-of-scope'
    std::swap( L[ IL[ v2 ] ], L[ no_v ] );
    //fix IL accordingly
    std::swap( IL[ L[IL[v2]] ], IL[ L[no_v] ] );
    assert( L[IL[v2]] == v2 );
    assert( IL[L[no_v]] == no_v );
    //adjust degrees of enlarged color v1
    CD_out[ v1 ] += CD_out[ v2 ];
    CD_out[ v2 ] = 0;
    
    //adjust CAL s.t. v2 is part of color list of v1
    //update vert cols
    for (auto &&v : CAL[v2]) VC[v] = v1;
    //update color list
    CAL[v1].splice(CAL[v1].end(), CAL[v2]);
    
    //perform the same steps for the symmetric nodes sigma(v1) and sigma(v2) -- if they are distinct!
    if(VC[SIGMA(v2)] != v1) {
        --no_v;
        //move SIGMA(v2) out-of-scope
        std::swap( L[ IL[ VC[SIGMA(v2)] ] ], L[ no_v ] );
        //fix IL accordingly
        std::swap( IL[ L[IL[VC[SIGMA(v2)]]] ], IL[ L[no_v] ] );
        assert( L[IL[VC[SIGMA(v2)]]] == VC[SIGMA(v2)] );
        assert( IL[L[no_v]] == no_v );
    
        //adjust degrees of enlarged color SIGMA(v1)
        CD_out[ VC[SIGMA(v1)] ] += CD_out[ VC[SIGMA(v2)] ];
        CD_out[ VC[SIGMA(v2)] ] = 0;
    
        //adjust CAL s.t. sigma(v2) is part of list of sigma(v1)
        //update vert cols
        var_t VC_sigma_v2 = VC[SIGMA(v2)];
        for (auto &&v : CAL[VC_sigma_v2]) VC[v] = VC[SIGMA(v1)];
        //update color list
        CAL[VC[SIGMA(v1)]].splice(CAL[VC[SIGMA(v1)]].end(), CAL[VC_sigma_v2]);
    }


    //adjust out-going edges, i.e., remove all duplicate edges between same colours
    //loop over verts of same color and remove each one of them
    exists_edge_to_c.resize(no_v,false);
    assert(std::all_of(exists_edge_to_c.begin(), exists_edge_to_c.end(), [](const bool v){ return !v; }));
    std::list<var_t> needs_reset;
    //std::fill(exists_edge_to_c.begin(),exists_edge_to_c.end(), false);
    //note IL[VC[w]] == IL[VC[w']] iff VC[w] == VC[w'], and we have IL[VC[w]] in [0,...,no_v-1]; hence we access color c via IL[c]
    exists_edge_to_c[ IL[v1] ] = true; //avoid edges from color to itself
    needs_reset.emplace_back( IL[v1] );
    for (const auto &w : CAL[ v1 ]) {
        //loop through out-edges
        for (var_t idx = 0; idx < VD_out[w]; ++idx) //cannot be run in parallel, as remove_edge decreases VD_out !
        {
            if(!exists_edge_to_c[ IL[ VC[AL_out[w][idx]] ] ]) {
                exists_edge_to_c[ IL[ VC[AL_out[w][idx]] ] ] = true; //only one edge is kept!
                needs_reset.emplace_back( IL[VC[AL_out[w][idx]]] );
            } else {
                //there is already an edge; i.e., remove it -- except it is self-symmetric (i.e. v->w and sigma(w)->sigma(v) coincide), then remove only one of the edges!
                if( VC[w] != VC[ SIGMA(AL_out[w][idx]) ] ) {
                    remove_edge(w, idx);
                } else {
                    remove_edge_(w, idx);
                }
                --idx; //decrease idx, as we still need to check whether AL_out[w][idx] needs to be removed (!)
            }
        }
    }
    //reset exists_edge_to_c
    for(const auto idx : needs_reset) exists_edge_to_c[ idx ] = false;
    needs_reset.clear();
    assert(std::all_of(exists_edge_to_c.begin(), exists_edge_to_c.end(), [](const bool v){ return !v; }));

    if(VC[SIGMA(v2)] != v1) {
        //adjust in-going edges, by remove all duplicate out-going edges between same colours of SIGMA(v1)
        //loop over verts of same color and remove each one of them
        exists_edge_to_c.resize(no_v, false);
        //note IL[VC[w]] == IL[VC[w']] iff VC[w] == VC[w'], and we have IL[VC[w]] in [0,...,no_v-1]; hence we access color c via IL[c]
        exists_edge_to_c[ IL[VC[SIGMA(v1)]] ] = true; //avoid edges from color to itself
        needs_reset.emplace_back( IL[VC[SIGMA(v1)]] );
        for (const auto &w : CAL[ VC[SIGMA(v1)] ]) {
            //loop through out-edges
            for (var_t idx = 0; idx < VD_out[w]; ++idx) //cannot be run in parallel, as remove_edge decreases VD_out !
            {
                if(!exists_edge_to_c[ IL[ VC[AL_out[w][idx]] ] ]) {
                    exists_edge_to_c[ IL[ VC[AL_out[w][idx]] ] ] = true; //only one edge is kept!
                    needs_reset.emplace_back( IL[VC[AL_out[w][idx]]] );
                } else {
                    //there is already an edge; i.e., remove it -- except it is self-symmetric (i.e. v->w and sigma(w)->sigma(v) coincide), then remove only one of the edges!
                    if( VC[w] != VC[ SIGMA(AL_out[w][idx]) ] ) {
                        remove_edge(w, idx);
                    } else {
                        remove_edge_(w, idx);
                    }
                    --idx; //decrease idx, as we still need to check whether AL_out[w][idx] needs to be removed (!)
                }
            }
        }
        //reset exists_edge_to_c
        for(const auto idx : needs_reset) exists_edge_to_c[ idx ] = false;
        needs_reset.clear();
        assert(std::all_of(exists_edge_to_c.begin(), exists_edge_to_c.end(), [](const bool v){ return !v; }));
    }
};

bool graph_lhgr::assert_data_structs() const noexcept {
    assert(no_e <= no_v*no_v-no_v);
    var_t total_d_orig_out = 0;
    var_t total_d_out = 0;
    for (var_t u = 0; u < L.size(); ++u) {
        //check validity of AL and AL_in (and JAL)
        //AL_out[ AL_in[u,i], IAL_in[u,i] ] = u for all u (not only active ones!) and all i (up to orig degree of u!)
        for (var_t i = 0; i < AL_out[u].size(); i++) {
            //std::cout << AL_out[SIGMA(AL_out[u)[i]] ][IAL_in[u][i]] << " == " << SIGMA(u) << std::endl;
            assert( AL_out[ SIGMA(AL_out[u][i]) ][IAL_in[u][i]] == SIGMA(u) );
        }
        //std::cout << std::endl;
        //get orig d_out
        const var_t d_orig_out = AL_out[u].size();

        //check upper bound on degree vecs
        assert(d_orig_out >= VD_out[u]);
        total_d_orig_out += d_orig_out;
        total_d_out += VD_out[u];

        //check sigma
        assert(SIGMA(SIGMA(u)) == u);
    };
    //total out-degree must be no_e
    assert( total_d_out <= total_d_orig_out );
    
    //check correct bounds of CAL, and resp cardinalities and degrees
    for (var_t c_idx = 0; c_idx < no_v; ++c_idx) {
        //get color at c_idx
        const var_t c = L[c_idx];
        var_t cd_out = 0;

        //loop over verts of same color
        for (const auto &v : CAL[c]) cd_out = cd_out + VD_out[v];

        assert( cd_out == CD_out[c] );
        assert( VC[c] == c);
    }

    //check that L and IL match!
    for (var_t u = 0; u < L.size(); ++u) {
        assert( L[ IL[u] ] == u && IL[ L[u] ] == u );
    }

    return true;
};

std::string graph_lhgr::to_str() const noexcept {
    std::map<var_t, vec<var_t> > edges;
    for (var_t c_idx = 0; c_idx < no_v; ++c_idx) {
        //get color at c_idx
        const var_t c = VC[ L[c_idx] ];
        //add all out-neighbors to edges[c_idx]
        edges[c] = get_out_neighbour_vector( c );
        //iterate over all verts of color c
        //for (auto &&v : CAL[c]) {
        //    //add out-neighbors of v to edges[c_idx]
        //    for (var_t i = 0; i < VD_out[v]; ++i) {
        //        edges[ c ].push_back( VC[ AL_out[v][i] ] );
        //    }
        //}
        assert( edges[c].size() == CD_out[c] );
        //sort color c_idx
        std::sort( std::execution::par, edges[c].begin(), edges[c].end() );
    }

    //generate string of edges with lexicographic ordering!
    vec< std::string > str_edges(0);
    for( std::map<var_t, vec<var_t>>::iterator iter = edges.begin(); iter != edges.end(); ++iter ) {
        vec<std::string> out_edges_str( iter->second.size() );
        //construct strings!
        const var_t src = iter->first;
        auto to_str = [src](const var_t dst) -> std::string {return "("+std::to_string(src)+","+std::to_string(dst)+")";};
        std::transform(std::execution::par, iter->second.begin(), iter->second.end(), out_edges_str.begin(), to_str);

        std::stringstream ss;
        std::copy(out_edges_str.begin(), out_edges_str.end(), std::ostream_iterator<std::string>(ss, " "));
        std::string out_edges = ss.str();
        if (!out_edges.empty()) {
            out_edges.resize(out_edges.length() - 1); // trim trailing space
            str_edges.push_back(out_edges);
        }
    }

    std::stringstream ss;
    std::copy(str_edges.begin(), str_edges.end(), std::ostream_iterator<std::string>(ss, "; "));
    std::string result = ss.str();
    if (!result.empty()) {
        result.resize(result.length() - 2); // trim trailing line feed
        return result;
    } else {
        return "";
    }
};