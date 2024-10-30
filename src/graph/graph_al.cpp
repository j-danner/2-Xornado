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

#include "graph_al.hpp"


graph_al::graph_al(const vec< std::pair<var_t,var_t> >& E, const var_t no_v_) noexcept : no_v(no_v_), no_e(2*E.size()) {
    init(E, no_v_);
};

void graph_al::init(const vec< std::pair<var_t,var_t> >& E, const var_t no_v_) noexcept {
    //init data structures!
    no_e = 0;
    no_v = no_v_;

    //active vertices
    L = vec<var_t>(no_v);
    IL = vec<var_t>(no_v);

    //adjacency list
    AL_out = vec<adj_l>(no_v);

    //fill data structures:

    //fill AL_out
    for (const auto &e : E) {
        const var_t src  = e.first;
        const var_t dst = e.second;
        //add edge
        //AL_out[src].push_back(dst);
        const auto [_,inserted] = AL_out[src].insert(dst);
        if(!inserted) continue;
        ++no_e;

        //add symmetric edge sigma(dst)->sigma(src) --- if distinct!
        if(SIGMA(dst) != src) {
            //AL_out[ SIGMA(dst) ].push_back( SIGMA(src) );
            const auto ins = AL_out[ SIGMA(dst) ].insert( SIGMA(src) );
            if(ins.second) ++no_e;
        }
    }

    //init L, IL, CAL
    for (var_t v = 0; v < no_v; ++v) {
        //fill L and IL
        L[v] = v;
        IL[v] = v;
    }

    assert( graph_al::assert_data_structs() );
};

void graph_al::backtrack(graph_al_repr&& graph_orig) noexcept { //TODO use move-ctor!
    //restore graph from graph_al_repr  
    no_v = std::move(graph_orig.no_v);
    no_e = std::move(graph_orig.no_e);
    //TODO!!
    AL_out = std::move(graph_orig.AL_out);
};


//removes idx-th going out from src as listed in AL_out_logical(src,idx) (ignores same-colored verts!)
void graph_al::remove_edge(const var_t src, const var_t dst) noexcept {
    [[maybe_unused]] auto el_er = AL_out[src].erase( dst );
    assert( el_er == 1 );
    no_e--;
    //remove symmetric counterpart
    no_e -= AL_out[SIGMA(dst)].erase( SIGMA(src) );
    
    assert( graph_al::assert_data_structs() );
    //TODO can we remove vert src or dst if one becomes isolated?!
};

//removes all out-going edges of src
void graph_al::remove_all_edges(const var_t v) noexcept {
    //remove all symmetric edges
    for(const auto& dst : AL_out[v]) {
        if(SIGMA(dst)!=v) {
            AL_out[SIGMA(dst)].erase( SIGMA(v) );
            no_e--;
        }
    }
    //clear vert!
    no_e -= AL_out[v].size();
    AL_out[v].clear();
};

//remove color/vert c and SIGMA(c)
void graph_al::remove_vert(const var_t c) noexcept {
    //std::cout << "no_v = " << std::to_string(no_v) << "; removing " << std::to_string(c) << " (" << std::to_string(SIGMA(c)) << ")" << std::endl;
    //remove c and SIGMA(c)
    for (const auto &v : std::vector<var_t>({c,SIGMA(c)})) {
        --no_v;
        //swap vertices in L and fix IL, move v 'out-of-scope'
        std::swap( L[ IL[v] ], L[ no_v ] );
        //fix IL accordingly
        std::swap( IL[ L[IL[v]] ], IL[ L[no_v] ] );

        remove_all_edges(v);
    };
};

//merges colors v1 and v2, color of v2 will be set to color of v1; and color of SIGMA(v2) set to color of SIGMA(v1)
void graph_al::merge_verts(const var_t v1, const var_t v2) noexcept {
    //std::cout << "no_v = " << std::to_string(no_v) << "; merging " << std::to_string(VC[v1]) << " and " << std::to_string(VC[v2]) << " (" << std::to_string(VC[SIGMA(v1)]) << " and " << std::to_string(VC[SIGMA(v2)]) << ")" << std::endl;
    //abort if 
    if(v1==v2 || IL[v1] >= no_v || IL[v2] >= no_v) return;
    assert(v1 != v2);
    assert(IL[v1] < no_v );
    assert(IL[v2] < no_v );
    assert( graph_al::assert_data_structs() );

    --no_v;
    //swap vertices in L to move v2 'out-of-scope'
    std::swap( L[ IL[ v2 ] ], L[ no_v ] );
    //fix IL accordingly
    std::swap( IL[ L[IL[v2]] ], IL[ L[no_v] ] );
    assert( L[IL[v2]] == v2 );
    assert( IL[L[no_v]] == no_v );
    
    //perform the same steps for the symmetric nodes sigma(v1) and sigma(v2) -- if they are distinct!
    if(SIGMA(v2) != v1) {
        --no_v;
        //move SIGMA(v2) out-of-scope
        std::swap( L[ IL[ SIGMA(v2) ] ], L[ no_v ] );
        //fix IL accordingly
        std::swap( IL[ L[IL[SIGMA(v2)]] ], IL[ L[no_v] ] );
        assert( IL[L[no_v]] == no_v );
    }

    //'bend' all incoming edges of v2 to go to v1 instead!
    for(const var_t v : get_in_neighbour_vector(v2)) { //cannot use range, as lazy eval produce bugs as to AL_out[SIGMA(v2)] being updated!
        AL_out[v].erase(v2);
        const auto ins = AL_out[v].insert(v1);
        if(!ins.second) no_e--;
    }
    if(SIGMA(v2) != v1) {
        //'bend' all incoming edges of SIGMA(v2) to go to SIGMA(v1)!
        for(const var_t v : get_in_neighbour_vector(SIGMA(v2))) { //cannot use range, as lazy eval produce bugs as to AL_out[v2] being updated!
            AL_out[v].erase(SIGMA(v2));
            const auto ins = AL_out[v].insert(SIGMA(v1));
            if(!ins.second) no_e--;
        }
    }

    //add out-edges of v1 to v2!
    AL_out[v1].merge(AL_out[v2]);
    //remove possible self-edges
    no_e -= AL_out[v1].erase(v1);
    no_e -= AL_out[v2].size();
    AL_out[v2].clear(); //adjacency list as set already avoids double entries! i.e., we can just clear AL_out[v2]!

    //perform same steps for symmetric edges!
    if(SIGMA(v2) != v1) {
        //add out-edges of SIGMA(v1) to SIGMA(v2)!
        AL_out[SIGMA(v1)].merge(AL_out[SIGMA(v2)]);
        no_e -= AL_out[SIGMA(v1)].erase(SIGMA(v1)); //remove possible self-edges
        no_e -= AL_out[SIGMA(v2)].size();
        AL_out[SIGMA(v2)].clear(); //adjacency list as set already avoids double entries! i.e., we can just clear AL_out[v2]!
    }
    
    assert( graph_al::assert_data_structs() );
};

bool graph_al::assert_data_structs() const noexcept {
    assert(no_e <= no_v*no_v);
    var_t total_d_out = 0;
    for (var_t u = 0; u < no_v; ++u) {
        total_d_out += AL_out[L[u]].size();
        for([[maybe_unused]] const auto dst : AL_out[L[u]]) assert(IL[dst] < no_v);
    };
    assert(total_d_out == no_e);

    //check sigma
    for (var_t u = 0; u < L.size(); ++u) {
        assert(SIGMA(SIGMA(u)) == u);
    };

    //check that L and IL match!
    for (var_t u = 0; u < L.size(); ++u) {
        assert( L[ IL[u] ] == u && IL[ L[u] ] == u );
    }

    return true;
};

std::string graph_al::to_str() const noexcept {
    std::map<var_t, vec<var_t> > edges;
    for (var_t c_idx = 0; c_idx < no_v; ++c_idx) {
        //get color at c_idx
        const var_t c = L[c_idx];
        //add all out-neighbors to edges[c_idx]
        edges[c] = get_out_neighbour_vector( c );
        //iterate over all verts of color c
        //for (auto &&v : CAL[c]) {
        //    //add out-neighbors of v to edges[c_idx]
        //    for (var_t i = 0; i < VD_out[v]; ++i) {
        //        edges[ c ].push_back( VC[ AL_out[v][i] ] );
        //    }
        //}
        //sort color c_idx
        std::sort( edges[c].begin(), edges[c].end() );
    }

    //generate string of edges with lexicographic ordering!
    vec< std::string > str_edges(0);
    for( std::map<var_t, vec<var_t>>::iterator iter = edges.begin(); iter != edges.end(); ++iter ) {
        vec<std::string> out_edges_str( iter->second.size() );
        //construct strings!
        const var_t src = iter->first;
        auto to_str = [src](const var_t dst) -> std::string {return "("+std::to_string(src)+","+std::to_string(dst)+")";};
        std::transform(iter->second.begin(), iter->second.end(), out_edges_str.begin(), to_str);

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
