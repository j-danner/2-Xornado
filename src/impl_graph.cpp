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
#include <stack>
#include <deque>
#include <queue>
#include <omp.h>
#include <stdexcept>

#include "robin_hood-3.11.5/robin_hood.h"
#include "impl_graph.hpp"

#include <unordered_set>

impl_graph::impl_graph(const vec< vec<lineral> >& clss_, const options& opt_) : graph(), opt(opt_) {
    vec< vec<lineral> > clss = clss_;

    bool repeat = true;
    while(repeat) {
        repeat = false;

        //init stacks
        graph_stack = std::stack< graph_repr >();
        xsys_stack = std::list< std::list<LinEqs> >();
        //init maps
        vl_stack = std::stack< vert_label_repr >();
    #ifndef FULL_REDUCTION
        //init assignments
        assignments = vec<lineral>(opt.num_vars+1, lineral());
    #endif

        auto E = vec< std::pair<var_t,var_t> >();
        E.reserve( (opt.ext==constr::extended ? 2 : 6) * clss.size() );

        //current number of verts
        var_t no_v = 0;

        //init stacks
        vl = vert_label(clss.size()* (opt.ext==constr::extended ? 6 : 2), opt.num_vars);
        //vl_stack.push( vl_trie(clss.size()*6, opt.num_vars) ); //TODO better guess on number of verts?!
        vec<lineral> _L = vec<lineral>();

        //run through xor-clauses to find lineq and construct first entries in trie, also generate sigma and E
        for (auto &cls : clss) {
            //we can only solve 2-XNFs!
            if(cls.size() == 1) { _L.push_back( cls.front() ); continue; }
            if(cls.size() > 2) throw std::runtime_error("given clauses are not in 2-XNF!");
            assert(cls.size() <=2);

            //now cls must consists of exactly 2 elements, say cls = {f,g}
            const lineral f = std::move(cls.front());
            const lineral g = std::move(cls.back());
            const lineral fpg = f+g;

            // if f=g, then we have a xor-literal
            if(f == g) _L.push_back( f );
            // if f = g+1, then the clause should be ignored
            if( (fpg).is_one() || f==g ) continue;
            //now clause is of degree 2 AND non-trivial

            lineral fp1 = f.plus_one();
            lineral gp1 = g.plus_one();
            lineral fpgp1 = fpg.plus_one();

            //generate f, g, f+1, g+1, f+g, f+g+1; their respective idxs, and sigma
            auto v_lits = (opt.ext==constr::extended) ? vec<lineral>({f, g, fpg}) : vec<lineral>({f, g});
            for (auto &l : v_lits)
            {
                if(vl.V_contains(l)) continue;
                if(l.has_constant()) l.add_one();
                auto ins = vl.insert(no_v, std::move(l), 0);
                no_v += ins.inserted ? 2 : 0;
            }
            //note: f*g = (f+g+1)*f = (f+g+1)*g, i.e. add edges:
            //now add edges: fp1 -> g; fpg -> f; fpg -> g   (symmetric edges are then constructed from sigma!)
            E.push_back( std::pair<var_t,var_t>( vl.V(fp1), vl.V(g) ) );
            if(opt.ext==constr::extended) {
                E.push_back( std::pair<var_t,var_t>( vl.V(fpg), vl.V(f) ) );
                E.push_back( std::pair<var_t,var_t>( vl.V(fpg), vl.V(g) ) );
            }
        }

        std::sort(E.begin(), E.end());
        std::unique(E.begin(), E.end());

        assert(no_v == 2*vl.size());

        //now E -- and its intepretations with V and Vxlit are constructed, so we can init the graph!
        init(E,no_v);

        if(s.cancelled.load())
        {
            VERB(10, "c cancelled during preprocessing");
            break;
        }

        //init backtrack-stacks
        xsys_stack.push_back( std::list<LinEqs>({ LinEqs(_L) }) );
        vl_stack.push( vl.get_state() );
        graph_stack.push( get_state() );

        //init activity_score, based on number of occurances as LTs
        activity_score = vec<unsigned int>(opt.num_vars+1, 1);
        //TODO run in parallel?! (or at least integreate in above construction of Vxlit?!)
        for (const auto &v : get_v_range()) {
            activity_score[ vl.Vxlit_LT(v) ]++;
        }
        assert( assert_data_structs() );
        
        //init done!
        VERB(70, graph_stats());
        //pre-process 
        preprocess();
        
        if(opt.pp != preproc::fls_scc_ee) continue;
        //if edge-extension is wanted; first perform crGCP
        //update clss to gcp vers
        clss.clear();
        clss = to_xcls(); //TODO inefficient!
        //compute topological order of graph
        const auto TO = get_TO();
    
        vec<LinEqs> D(no_v);
        for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
            var_t v = *v_it;
            const lineral f = vl.Vxlit(v);
            D[IL[v]] += LinEqs(f);
            for (const auto &w : get_in_neighbour_range(v)) D[IL[w]] += D[IL[v]];
            //D[IL[v]].clear();
        }
        //add new edges
        const auto roots = get_roots();
        var_t c_new_edges = 0;
        for(const auto& r1 : roots) {
            for(const auto& r2 : roots) {
                if(r1==r2) continue;
                const LinEqs tmp = D[IL[r1]] + D[IL[r2]];
                if(tmp.is_consistent()) continue;
                //const auto& [b,a] = intersectaffineVS(D[IL[r1]], D[IL[r2]]);
                //if(!b) continue;
                if(is_descendant(r1, SIGMA(r2))) continue;
                repeat = true;
                clss.emplace_back( vec<lineral>({vl.Vxlit(SIGMA(r1)), vl.Vxlit(SIGMA(r2))}) );
            }
        }
        VERB(70, "c deduced "+std::to_string(c_new_edges)+" new edges!");
        VERB(110, to_str());
    } while(repeat);
};

#ifdef USE_TRIE
    LinEqs impl_graph::update_graph(stats& s, const LinEqs& L) {
        return update_graph_hash_fight_dev(s,L);
    };
    
    LinEqs impl_graph::update_graph_par(stats& s, const LinEqs& L) {
        return update_graph(s,L);
    };
    
    //one global var for lineral lit, reduces memory allocations
    lineral lit;
    LinEqs impl_graph::update_graph_hash_fight(stats& s, const LinEqs& L) {
        s.no_graph_upd++;
        s.total_upd_no_v += no_v;
        s.total_upd_xsys_size += L.size();

        //return std::move( update_graph(s, L) );
        assert( L.is_consistent() );
        auto new_L = vec<lineral>();
        if(L.size() == 0) return LinEqs(std::move(new_L));
    
        //update in two steps:
        // (1) in parallel: reduce all lits and create new map Vxlit (linerals -> verts)
        // (2) sequentially: for each vertex v check whether 'V_stack.top().at( Vxlit(v) ) == v'; if not merge correspondingly!
    
        std::list<std::pair<var_t,var_t> > merge_list;
        //update current trie!
        //  -- (1) --
        for (auto v : get_v_range()) {
        #ifdef FULL_REDUCTION
            if(!vl.contains(v)) continue;
            lit = std::move( vl.Vxlit(v) );
            //reduce with linsys
            const bool update_req = lit.reduce( L );
        #else
            if(!vl.contains(v) || assignments[vl.Vxlit_LT(v)].is_zero()) continue;
            lit = std::move( vl.Vxlit(v) );
            //reduce with linsys
            const bool update_req = lit.lt_reduce( assignments );
            assert(lit.LT() != vl.Vxlit_LT(v));
        #endif
            if(update_req) {
                s.no_vert_upd++;
                //insert reduced lit in new_trie
                auto [v_upd,b] = vl.update(v, lit, get_dl());
                if( v_upd != v) merge_list.push_back( std::pair<var_t,var_t>{v, b ? SIGMA(v_upd) : v_upd} );
            }
        }
    
        // -- (2) --
        for (const auto &[v, v_] : merge_list) merge_verts(v_, v);
    
        //check if we can deduce linerals: i.e. check for literal '0'
        const auto [has_zero,is_one,v_] = vl.if_exists_get_zero_v();
        if(has_zero) {
            //flip if v_zero is actually vertex for 1 (!)
            var_t v_zero = is_one ? SIGMA(v_) : v_;
            for (const auto &v : get_out_neighbour_range( v_zero )) {
                new_L.push_back( vl.Vxlit(v) );
            #ifndef FULL_REDUCTION
                new_L.back().reduce(assignments);
            #endif
            }
            //remove vertex!
            remove_vert( v_zero );
            //remove lit/v from trie 
            vl.erase(v_);
        }
    
        assert( assert_data_structs() );
    
        return LinEqs(std::move(new_L));
    };
    
    LinEqs impl_graph::update_graph_hash_fight_dev(stats& s, const LinEqs& L) {
        s.no_graph_upd++;
        s.total_upd_no_v += (long) no_v;
        s.total_upd_xsys_size += (long) L.size();
    
        //return std::move( update_graph(L) );
        assert( L.is_consistent() );
        auto new_L = vec<lineral>();
        if(L.size() == 0) return LinEqs(std::move(new_L));
    
        //if(!assignments[9].is_zero()) {
        //    lit.reset();
        //}

        //update in two steps:
        // (1) in parallel: reduce all lits and create new map Vxlit (linerals -> verts)
        // (2) sequentially: for each vertex v check whether 'V_stack.top().at( Vxlit(v) ) == v'; if not merge correspondingly!
    
        std::list<std::pair<var_t,var_t> > merge_list;
        //update current trie!
        //  -- (1) --
        for (auto v : get_v_range()) {
        #ifdef FULL_REDUCTION
            if(!vl.contains(v)) continue;
            lit = std::move( vl.Vxlit(v) );
            //reduce with linsys
            const bool update_req = lit.reduce( L );
            if(update_req) {
                s.no_vert_upd++;
                //insert reduced lit in new_trie
                auto [v_upd,b] = vl.update(v, lit, get_dl());
                if( v_upd != v) merge_list.push_back( std::pair<var_t,var_t>{v, b ? SIGMA(v_upd) : v_upd} );
            }
        #else
            if(!vl.contains(v) || assignments[vl.Vxlit_LT(v)].is_zero()) continue;
            lit = std::move( vl.Vxlit(v) );
            assert(lit.LT() == vl.Vxlit_LT(v));
            //reduce with linsys
            const bool update_req = lit.lt_reduce( assignments );
            assert(lit.LT() != vl.Vxlit_LT(v));
            assert(update_req);
            s.no_vert_upd++;
            //insert reduced lit in new_trie
            auto [v_upd,b] = vl.update(v, lit, get_dl());
            if( v_upd != v) {
                merge_list.emplace_back( v, b ? SIGMA(v_upd) : v_upd );
                assert( lit == vl.Vxlit(merge_list.back().second) );
            }
        #endif
        }
    
        // -- (2) --
        for (const auto &[v, v_] : merge_list) merge_verts(v_, v);
    
        //check if we can deduce linerals: i.e. check for literal '0'
        const auto [has_zero,is_one,v_] = vl.if_exists_get_zero_v();
        if(has_zero) {
            //flip if v_zero is actually vertex for 1 (!)
            var_t v_zero = is_one ? SIGMA(v_) : v_;
            //collect all linerals reachable by v_zero
            robin_hood::unordered_flat_set<var_t> marked;
            std::list<var_t> q;
            q.push_back( v_zero );
            marked.insert( v_zero );
            while(!q.empty()) {
                const auto v = q.front();
                q.pop_front();
                //get corr lineral -- if SIGMA(v) is marked this means that l AND l+1 are implied, i.e., add 1 to new_L (!)
                if(marked.contains(SIGMA(v))) {
                    new_L.push_back( lineral(cnst::one) );
                    break;
                }
                //now v is certainly still contained in vl
                assert( vl.contains(v) || vl.contains(SIGMA(v)));
                new_L.emplace_back( std::move(vl.Vxlit(v)) );
            #ifndef FULL_REDUCTION
                new_L.back().reduce(assignments);
            #endif
                //put children in queue - if not yet marked
                for (const auto &w : get_out_neighbour_range( v )) {
                    if(!marked.contains(w)) {
                        q.push_back(w);
                        marked.insert(w); //ensure that w can only appears once in q
                    }
                }
                //remove vertex!
                assert( IL[v] < no_v );
                remove_vert( v );
                //remove lit/v from trie 
                vl.erase( vl.contains(v) ? v : SIGMA(v) );

                assert( assert_data_structs() );
            }
        }
    
        assert( assert_data_structs() );
        
        return LinEqs(std::move(new_L));
        //return std::move( update_graph(L) );
    };
    
#else
    LinEqs impl_graph::update_graph(stats& s, const LinEqs& L) {
        assert( L.is_consistent() );
        auto new_L = vec<lineral>();
        if(L.size() == 0) return std::move( LinEqs(std::move(new_L)) );

        //update in two steps:
        // (1) in parallel: reduce all lits and create new map Vxlit (linerals -> verts)
        // (2) sequentially: for each vertex v check whether 'V_stack.top().at( Vxlit(v) ) == v'; if not merge correspondingly!

        vec<std::pair< var_t,var_t> > merge_vec(no_v);
        //compute new V and Vxlit
        hmap<lineral, var_t> new_V(no_v/10);
        hmap<var_t, lineral> new_Vxlit(no_v/10);
        //  -- (1) -- parallel
        //#pragma omp parallel for
        for (auto v : get_v_range()) { //TODO update v then update SIGMA(v) based on v --> saves half the time!
            if(!vl.contains(v)) continue;
            lineral lit = vl.Vxlit( v );
            //reduce with linsys
            #ifdef FULL_REDUCTION
                const bool update_req = lit.reduce( L );
            #else
                const bool update_req = lit.lt_reduce( assignments );
            #endif
            if(lit.has_constant()) {
                //flip lit and v
                v = SIGMA(v);
                lit.add_one();
            }
            //TODO only update if update_req is true!
            if(update_req) s.no_vert_upd++;
            //insert new lineral in new_V
            auto inserted = new_V.emplace( lit, v );
            if(inserted.second) {
                //if inserted, then update Vxlit correspondingly
                //auto inserted2 = new_Vxlit.insert( v, lit );
                [[maybe_unused]] auto inserted2 = new_Vxlit.emplace( v, std::move(lit) ); 
                assert(inserted2.second);
            }
            merge_vec[ IL[v] ] = std::pair<var_t,var_t>{v, inserted.second ? v : (inserted.first)->second};
        }

        // -- (2) -- sequential
        for (const auto &[v, v_] : merge_vec) merge_verts(v_, v);

        //check if we can deduce linerals: i.e. check for literal '0'
        const auto search = new_V.find( lineral() );
        if( search != new_V.end() ) { //0-literal found!
            for (const auto &v : get_out_neighbour_range( search->second )) {
                new_L.push_back( new_Vxlit.contains(v) ? new_Vxlit.at(v) : new_Vxlit.at(SIGMA(v)).plus_one() );
            #ifndef FULL_REDUCTION
                new_L.back().reduce(assignments);
            #endif
            }
            //remove vertex!
            remove_vert( search->second );
            //remove lit/v from new_Vxlit and new_V (!)
            new_Vxlit.erase( search->second );
            new_V.erase( search );
        }
        
        vl.put_Vxlit( std::move(new_Vxlit) );
        vl.put_V( std::move(new_V) );

        assert( assert_data_structs() );

        return std::move( LinEqs(std::move(new_L)) );
    };

    LinEqs impl_graph::update_graph_par(stats& s, const LinEqs& L) {
        return std::move( update_graph(s,L) );
    };

    LinEqs impl_graph::update_graph_hash_fight(stats& s, const LinEqs& L) {
        return std::move( update_graph(s,L) );
    };
    
    LinEqs impl_graph::update_graph_hash_fight_dev(stats& s, const LinEqs& L) {
        return std::move( update_graph(s,L) );
    };

#endif


//implication graph analysis!
// C++ Implementation of Kosaraju's algorithm to print all SCCs -- adapted from 'https://www.geeksforgeeks.org/strongly-connected-components/'
void impl_graph::scc_fillOrder(const var_t v, vec<bool>& visited, std::stack<var_t> &Stack) const {
    // Mark the current node as visited and print it
    visited[IL[v]] = true;
    // Recur for all the vertices adjacent to this vertex
    for (const auto &w : get_out_neighbour_range(v)) {
        if(!visited[IL[w]]) scc_fillOrder(w, visited, Stack);
    }
    // All vertices reachable from v are processed by now, push v
    Stack.push(v);
}

// A recursive function to print reverse DFS starting from v (with root rt), i.e., in transpose graph and reading visited-vector negated
void impl_graph::scc_dfs_util(const var_t rt, vec<lineral>& linerals, var_t v, vec<bool>& visited, std::list< std::pair<var_t,var_t> >& merge_list) const {
    // Mark the current node as visited and print it
    visited[IL[v]] = false;
    if(v!=rt) {
        linerals.emplace_back( std::move( Vxlit_sum(rt,v) ) );
        merge_list.emplace_back( rt, v );
    }
    // Recur for all the vertices in-coming to this vertex
    for (const auto &&w : get_in_neighbour_range(v)) {
        if (visited[IL[w]]) scc_dfs_util(rt, linerals, w, visited, merge_list);
    }
}

// The main function that finds and prints all strongly connected
// components
// TODO rewrite complete tarjans algorithm to also check for trivial FLS ?
vec<lineral> linerals;
vec<lineral> out;
LinEqs impl_graph::scc_analysis() {
    std::stack<var_t> Stack;

    // Mark all the vertices as not visited (For first DFS)
    vec<bool> visited( no_v, false );

    // Fill vertices in stack according to their finishing times
    for (const auto &w : get_v_range()) {
        if(visited[IL[w]] == false) scc_fillOrder(w, visited, Stack);
    }

    //new linerals:
    linerals.clear();
    std::list< std::pair<var_t,var_t> > merge_list;
    // Now process all vertices in order defined by Stack
    while (Stack.empty() == false) {
        // Pop a vertex from stack
        const var_t v = Stack.top();
        Stack.pop();
        // get SCC of v
        if (visited[IL[v]] == true) {
            scc_dfs_util(v, linerals, v, visited, merge_list);
            //minor optimization to finding only one of each symmetrical components!
            visited[IL[SIGMA(v)]] = false;
        }
    }

#ifndef FULL_REDUCTION
    //filter out all already known linerals
    std::for_each(linerals.begin(), linerals.end(), [&](lineral& l){ l.reduce(assignments); } );
    out.clear();
    std::copy_if(linerals.begin(), linerals.end(), std::back_inserter(out), [](const lineral l){ return !l.is_zero(); });
    std::swap(linerals,out);
#endif

    LinEqs scc = LinEqs(std::move(linerals));
    //merge nodes if scc is consistent!
    if(scc.is_consistent()) {
        //merge SCCs (and remove labels)
        for(const auto& [rt,v] : merge_list) {
            merge_verts(rt, v);
            //erase label from vl
            if(vl.contains(v)) vl.erase(v);
            else if(vl.contains(SIGMA(v))) vl.erase(SIGMA(v)); 
        }
    }

    return scc;
};


#define IS_MARKED(a,i) !(a[i]==L.size())

vec<var_t> impl_graph::label_components() const {
    assert(no_v <= L.size());
    //construct map that maps vert number to root vert of same comp
    vec<var_t> label(no_v, L.size());

    std::list<var_t> q;
    for(const auto rt : get_v_range()) {
        if(IS_MARKED(label,IL[rt])) continue;
        //mark all verts reachable from rt via bfs
        q.push_back(rt);
        while(!q.empty()) {
            const auto v = q.front();
            q.pop_front();
            if(IS_MARKED(label,IL[v])) continue;
            label[IL[v]] = rt;
            for(const auto n : get_in_neighbour_range(v))  if(!IS_MARKED(label,IL[n])) q.push_back(n);
            for(const auto n : get_out_neighbour_range(v)) if(!IS_MARKED(label,IL[n])) q.push_back(n);
        }
    }
    
    return label;
}

var_t impl_graph::get_number_connected_components() const {
    const auto label = label_components();
    //number of connected components is exactly number of v's with label[IL[v]] = v
    var_t n_cc = 0;
    for(const auto& v : get_v_range()) if(label[IL[v]] == v) n_cc++;

    return n_cc;
}

std::list<var_t> impl_graph::get_roots() const {
    std::list<var_t> q_roots;
    
    // initialize the queue with all the vertices with no inbound edges
    for (const auto &v : get_v_range()) {
        if( get_in_degree(v) == 0 ) q_roots.push_back(v);
    }

    return q_roots;
}

// computes a topological order of the graph; assertion failure if graph is no DAG (!)
vec<var_t> impl_graph::get_TO() const {
    std::queue<var_t> q_roots;
    vec<var_t> in_degree_tmp(no_v);
    
    // initialize the queue with all the vertices with no inbound edges
    for (const auto &v : get_v_range()) {
        in_degree_tmp[IL[v]] = get_in_degree(v);
        if( get_in_degree(v) == 0 ) q_roots.push(v);
    }
    
    vec<var_t> to;
    to.reserve( no_v );
    while (!q_roots.empty()) {
        const var_t v = q_roots.front();
        q_roots.pop();
        to.emplace_back(v);
        // 'remove' v from the graph, by decreasing in_degree_tmp corr
        for (const auto &w : get_out_neighbour_range(v)) {
            --in_degree_tmp[IL[w]];
            if(in_degree_tmp[IL[w]] == 0) {
                q_roots.push(w);
            }
        }
    }
    //if size of result is not no_e, then there was a cycle, i.e., we had not DAG!
    if(to.size() < no_v) to.clear();
    assert(to.size() == 0 || to.size() == no_v);
    return to;
} 


bool impl_graph::is_descendant(const var_t src, const var_t dst) {
    if(src==dst) return true;
    std::list<var_t> queue;
    queue.push_back(src);
    while(!queue.empty()) {
        const var_t v = queue.back();
        queue.pop_back();
        for(const auto& n : get_out_neighbour_range(v)) {
            if(n==dst) return true;
            queue.push_back(n);
        }
    }
    return false;
};

//decision heuristics
std::pair< LinEqs, LinEqs > impl_graph::first_vert() const {
    //assert(no_v>0);
    //return std::move( std::pair< LinEqs, LinEqs >( LinEqs( std::move( vl.Vxlit(L[0]) ) ), LinEqs( std::move( vl.Vxlit(L[0]).add_one() ) ) ) );
    //guess single ind
    var_t i = 0;
    var_t lt = 0;
    while(lt == 0) {
        assert(no_v>i);
        lt = vl.Vxlit_LT(L[i]); 
        ++i;
    }
    lineral lt_lit(vec<var_t>({lt}));
#ifndef FULL_REDUCTION
    assert(assignments[lt].is_zero());
#endif
    return std::pair< LinEqs, LinEqs >( LinEqs( lt_lit ), LinEqs( lt_lit.plus_one() ) );
}

vec<lineral> tree_xlits;
vec<bool> marked;
std::stack<var_t> queue;
std::pair< LinEqs, LinEqs > impl_graph::max_reach() const {
    //find max tree by traversing TO in reverse
    //TODO can we adapt code from TO to do it in one go? (otherwise two full traversals of the graph are necessary...)
    vec<int> tree_score(no_v, 1);
    //init tree_score
    if(opt.score == sc::active) {
        for (const auto &v : get_v_range()) {
            tree_score[IL[v]] = activity_score[ vl.Vxlit_LT(v) ];
        }
    }
    
    var_t v_max_tree = L[0];
    //compute topological ordering of graph
    const auto TO = get_TO();
    //traverse through it in reverse, and calculate tree_size
    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
        const var_t v = *v_it;
        for (const auto &w : get_out_neighbour_range(v)) {
            tree_score[IL[v]] += tree_score[IL[w]];
        }
        if(tree_score[IL[v]] > tree_score[IL[v_max_tree]]) v_max_tree = v;
    }

    //compute tree LinEqs (all out-neighbours):
    tree_xlits.clear();
    marked.clear();
    marked.resize(no_v, false);
    while(!queue.empty()) queue.pop();
    queue.push(v_max_tree);
    marked[IL[v_max_tree]] = true;
    while(!queue.empty()) {
        var_t v = queue.top();
        queue.pop();
        //fast exit if SIGMA(v) was already marked!
        if(marked[IL[SIGMA(v)]]) {
            tree_xlits.clear();
            tree_xlits.emplace_back( lineral( cnst::one ) );
            break;
        };
        tree_xlits.emplace_back( std::move( vl.Vxlit(v) ) );
    #ifndef FULL_REDUCTION
        tree_xlits.back().reduce(assignments);
    #endif
        for(const auto &w : get_out_neighbour_range(v)) {
            if(!marked[IL[w]]) {
                marked[IL[w]] = true;
                queue.push(w);
            }
        }
    }
    assert(tree_xlits.size() < no_v); //ensure that no linerals are stored multiple times!
    const LinEqs tree_xsys = LinEqs( std::move(tree_xlits) );
    tree_xlits.clear();

    //compute inv-tree LinEqs (all in-neighbours):
    while(!queue.empty()) queue.pop();
    queue.push(v_max_tree);
    marked.clear();
    marked.resize(no_v, false);
    marked[IL[v_max_tree]] = true;
    while(!queue.empty()) {
        var_t v = queue.top();
        queue.pop();
        //fast exit if SIGMA(v) was already marked!
        if(marked[IL[SIGMA(v)]]) {
            tree_xlits.clear();
            tree_xlits.emplace_back( lineral( cnst::one ) );
            break;
        };
        tree_xlits.emplace_back( std::move( vl.Vxlit(v).add_one() ) );
    #ifndef FULL_REDUCTION
        tree_xlits.back().reduce(assignments);
    #endif
        for(const auto &w : get_in_neighbour_range(v)) {
            if(!marked[IL[w]]) {
                marked[IL[w]] = true;
                queue.push(w);
            }
        }
    }
    const LinEqs inv_tree_xsys = LinEqs( std::move(tree_xlits) );

    return std::pair< LinEqs, LinEqs >( std::move( tree_xsys ), std::move( inv_tree_xsys ) );
}

std::pair< LinEqs, LinEqs > impl_graph::max_bottleneck() const {
    //find max tree by traversing TO in reverse
    //TODO can we adapt code from TO to do it in one go? (otherwise two full traversals of the graph are necessary...)
    vec<int> bn_in_score(no_v, 1);
    vec<int> bn_out_score(no_v, 1);
    //init tree_score
    for (const auto &v : get_v_range()) {
        bn_in_score[IL[v]] = activity_score[ vl.Vxlit_LT(v) ];
        bn_out_score[IL[v]] = activity_score[ vl.Vxlit_LT(v) ];
    }
    
    var_t v_max_bn = L[0];
    //compute topological ordering of graph
    const auto TO = get_TO();
    //traverse through it twice, one time in order and once in reverse; in the latter find maximum
    for (auto v_it = TO.begin(); v_it != TO.end(); ++v_it) {
        const var_t v = *v_it;
        //TODO bn_score is actually an upper bound on the score; as 1 -> 2 -> 3 and 1 -> 3, makes 3 seem to have 3 incoming nodes, '1' is counted twice!
        for (const auto &w : get_in_neighbour_range(v)) {
            bn_in_score[IL[v]] += bn_in_score[IL[w]];
        }
    }
    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
        const var_t v = *v_it;
        //TODO bn_score is actually an upper bound on the score; as 1 -> 2 -> 3 and 1 -> 3, makes 3 seem to have 3 incoming nodes, '1' is counted twice!
        for (const auto &w : get_in_neighbour_range(v)) {
            bn_out_score[IL[v]] += bn_out_score[IL[w]];
        }
        //write total score into bn_in_score
        bn_in_score[IL[v]] += bn_out_score[IL[v]];
        if(bn_in_score[IL[v]] > bn_in_score[IL[v_max_bn]]) v_max_bn = v;
    }
    
    //compute tree LinEqs (all out-neighbours):
    tree_xlits.clear();
    while(!queue.empty()) queue.pop();
    queue.push(v_max_bn);
    marked.clear();
    marked.resize(no_v);
    while(!queue.empty()) {
        var_t v = queue.top();
        if(marked[IL[SIGMA(v)]]) {
            tree_xlits.clear();
            tree_xlits.emplace_back( lineral( cnst::one ) );
            break;
        };
        queue.pop();
        tree_xlits.emplace_back( std::move( vl.Vxlit(v) ) );
    #ifndef FULL_REDUCTION
        tree_xlits.back().reduce(assignments);
    #endif
        for(const auto &w : get_out_neighbour_range(v)) {
            if(!marked[IL[w]]) queue.push(w);
        }
    }
    const LinEqs tree_xsys = LinEqs( std::move(tree_xlits) );
    tree_xlits.clear();

    //compute inv-tree LinEqs (all in-neighbours):
    while(!queue.empty()) queue.pop();
    queue.push(v_max_bn);
    marked.clear();
    marked.resize(no_v);
    while(!queue.empty()) {
        var_t v = queue.top();
        if(marked[IL[SIGMA(v)]]) {
            tree_xlits.clear();
            tree_xlits.emplace_back( lineral( cnst::one ) );
            break;
        };
        queue.pop();
        tree_xlits.emplace_back( std::move( vl.Vxlit(v).add_one() ) );
    #ifndef FULL_REDUCTION
        tree_xlits.back().reduce(assignments);
    #endif
        for(const auto &w : get_in_neighbour_range(v)) {
            if(!marked[IL[w]]) queue.push(w);
        }
    }
    const LinEqs inv_tree_xsys = LinEqs( std::move(tree_xlits) );

    
    return std::pair< LinEqs, LinEqs >( std::move( tree_xsys ), std::move( inv_tree_xsys ) );
}

vec<bool> assigned;
std::pair< LinEqs, LinEqs > impl_graph::lex() const {
    assigned.resize(no_v, false);
    for(const auto& l_lineqs : xsys_stack) {
        for(const auto& l : l_lineqs) {
            for(const auto &[lt,idx] : l.get_pivot_poly_idx()) {
                assigned[lt] = l.get_linerals(idx).size() == 1;
            }
        }
    }
    for(var_t i=1; i<opt.num_vars; ++i) {
        if(!assigned[i]) {
            //guess single ind
            var_t lt = i;
            lineral lt_lit(vec<var_t>({lt}));
            return std::pair< LinEqs, LinEqs >( LinEqs( lt_lit ), LinEqs( lt_lit.plus_one() ) );
        }
    }
    assert(false);
}

std::pair< LinEqs, LinEqs > impl_graph::max_path() const {
    //find max path by traversing TO in reverse
    if(no_e == 0) return first_vert(); //contains no edges, i.e., longest path is of length 1, i.e., we guess a single vertex!
    //TODO can we adapt code from TO to do it in one go? (otherwise two full traversals of the graph are necessary...)
    vec<var_t> path_length(no_v, 1);
    vec<var_t> path_next(no_v);
    var_t v_max_path_src = L[0];
    //compute topological ordering of graph
    const auto TO = get_TO();
    //traverse through it in reverse, and calculate tree_size
    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
        const var_t v = *v_it;
        if(get_out_degree(v) == 0) path_next[IL[v]] = v;
        assert( path_length[IL[v]] == 1 );
        for (const auto &w : get_out_neighbour_range(v)) {
            //if incoming path is longer
            if(path_length[IL[w]]+1 > path_length[IL[v]]) {
                path_length[IL[v]] = path_length[IL[w]]+1;
                path_next[IL[v]] = w;
            }
        }
        if(path_length[IL[v]] > path_length[IL[v_max_path_src]]) v_max_path_src = v;
    }

    assert(path_length[IL[v_max_path_src]] > 1);
    VERB(40, "c chosen path has length " << path_length[IL[v_max_path_src]])
    //construct cycle and no_cycle LinEqs:
    auto cycle_xlits = vec<lineral>();
    cycle_xlits.reserve(path_length[IL[v_max_path_src]]);
    var_t v = v_max_path_src;
    for (var_t i = 0; i < path_length[IL[v_max_path_src]]; i++) {
        //cycle_xlits.emplace_back( std::move( vl.Vxlit(v)+vl.Vxlit(path_next[IL[v]]) ) );
        cycle_xlits.emplace_back( std::move( Vxlit_sum(v,path_next[IL[v]]) ) );
    #ifndef FULL_REDUCTION
        cycle_xlits.back().reduce(assignments);
    #endif
        //go to next vert in max-path
        v = path_next[IL[v]];
    }
    //max path is from f=src_lit to g=prev_lit; guess g->f, giving lits in cycle_xlits, or < (g+1)*f+1 >=< g,f+1 >
#ifdef FULL_REDUCTION
    const LinEqs no_cycle( { std::move( vl.Vxlit(v_max_path_src).add_one() ), std::move( vl.Vxlit(v) ) } );
#else
    vec<lineral> no_cycle_lits = { std::move( vl.Vxlit(v_max_path_src).add_one() ), std::move( vl.Vxlit(v) ) };
    std::for_each(no_cycle_lits.begin(), no_cycle_lits.end(), [&](lineral& l){ l.reduce(assignments); } );
    const LinEqs no_cycle( no_cycle_lits );
#endif
    const LinEqs cycle( std::move(cycle_xlits) );
    //if( cycle.size() < no_cycle.size() ) std::swap( cycle, no_cycle ); //NOTE negative impact on performance with random instances!
#ifdef FULL_REDUCTION
    return std::pair<LinEqs,LinEqs>(std::move(cycle), std::move(no_cycle));
#else
    //we might have already made this exact guess (or one that 'contains' it), however it did not 'remove' the current path and make it an SCC (as it should have!); so we need to guess differently! --> use max_reach heuristic, which certainly leads to propagation and reductions!
    if(no_cycle.size()==0 || cycle.size()==0) return max_reach();
    else return std::pair<LinEqs,LinEqs>(std::move(cycle), std::move(no_cycle));
#endif
}

std::pair< LinEqs, LinEqs > impl_graph::max_score_path() const {
    //find max path by traversing TO in reverse
    if(no_e == 0) return first_vert(); //contains no edges, i.e., longest path is of length 1, i.e., we guess a single vertex!
    //TODO can we adapt code from TO to do it in one go? (otherwise two full traversals of the graph are necessary...)
    vec<var_t> path_score(no_v, 1);
    vec<var_t> path_length(no_v, 1);
    vec<var_t> path_next(no_v);
    var_t v_max_path_src = L[0];
    //compute topological ordering of graph
    const auto TO = get_TO();
    //traverse through it in reverse, and calculate tree_size
    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
        const var_t v = *v_it;
        const int vert_score = activity_score[ vl.Vxlit_LT(v) ];
        path_score[IL[v]] = vert_score;
        if(get_out_degree(v) == 0) { 
            path_next[IL[v]] = v;
        } else {
            int best_out_score = 0;
            //find best outgoing path
            for (const auto &w : get_out_neighbour_range(v)) {
                if(path_score[IL[w]] > best_out_score) {
                    best_out_score = path_score[IL[w]];
                    path_next[IL[v]] = w;
                    path_length[IL[v]] = path_length[IL[w]] + 1;
                }
            }
            path_score[IL[v]] += best_out_score;
            //path_score[IL[v]] *= best_out_score;
        }
        if(path_score[IL[v]] > path_score[IL[v_max_path_src]]) v_max_path_src = v;
    }

    VERB(40, "c chosen path has score " << path_score[IL[v_max_path_src]] << " and length " << path_length[IL[v_max_path_src]])
    assert(path_score[IL[v_max_path_src]] > 0);
    if(path_length[IL[v_max_path_src]] > 1) {
        //construct cycle and no_cycle LinEqs:
        auto cycle_xlits = vec<lineral>();
        cycle_xlits.reserve(path_length[IL[v_max_path_src]]);
        var_t v = v_max_path_src;
        for (var_t i = 0; i < path_length[IL[v_max_path_src]]; i++) {
            //cycle_xlits.emplace_back( std::move( vl.Vxlit(v)+vl.Vxlit(path_next[IL[v]]) ) );
            cycle_xlits.emplace_back( std::move( Vxlit_sum(v,path_next[IL[v]]) ) );
        #ifndef FULL_REDUCTION
            cycle_xlits.back().reduce(assignments);
        #endif
            //go to next vert in max-path
            v = path_next[IL[v]];
        }
        //max path is from f=src_lit to g=prev_lit; guess g->f, giving lits in cycle_xlits, or < (g+1)*f+1 >=< g,f+1 >
    #ifdef FULL_REDUCTION
        const LinEqs no_cycle( { std::move( vl.Vxlit(v_max_path_src).add_one() ), std::move( vl.Vxlit(v) ) } );
    #else
        vec<lineral> no_cycle_lits = { std::move( vl.Vxlit(v_max_path_src).add_one() ), std::move( vl.Vxlit(v) ) };
        std::for_each(no_cycle_lits.begin(), no_cycle_lits.end(), [&](lineral& l){ l.reduce(assignments); } );
        const LinEqs no_cycle( no_cycle_lits );
    #endif
        const LinEqs cycle( std::move(cycle_xlits) );
        //if( cycle.size() < no_cycle.size() ) std::swap( cycle, no_cycle ); //NOTE negative impact on performance with random instances!
    #ifdef FULL_REDUCTION
        return std::pair<LinEqs,LinEqs>(std::move(cycle), std::move(no_cycle));
    #else
        //we might have already made this exact guess (or one that 'contains' it), however it did not 'remove' the current path and make it an SCC (as it should have!); so we need to guess differently! --> use max_reach heuristic, which certainly leads to propagation and reductions!
        if(cycle.size()==0) return max_reach();
        else return std::pair<LinEqs,LinEqs>(std::move(cycle), std::move(no_cycle));
    #endif
    } else {
        //path has length 1 (!), i.e., we only guess a single vert!
        return std::pair< LinEqs, LinEqs >( LinEqs( std::move( vl.Vxlit( v_max_path_src).add_one() ) ), LinEqs( std::move( vl.Vxlit( v_max_path_src) ) ) );
    }
}

//in-processing
LinEqs impl_graph::fls_no() const {
    return LinEqs();
};

LinEqs impl_graph::fls_trivial_cc() const {
    //((1)) label connected components
    const auto label = label_components();

    //((2)) find FLS starting only from roots r where r and sigma[r] are in the same component
    //(1) compute roots
    //auto roots = get_roots();
    std::list<var_t> roots;
    // initialize the queue with all the vertices with no inbound edges
    for (const auto &v : get_v_range()) {
        if( (get_in_degree(v) == 0) && (label[IL[v]] == label[IL[SIGMA(v)]]) ) roots.push_back(v);
    }

    //(2) from every root r perform dfs; on discovery of v mark it with r and check if SIGMA(v) was already marked with r;
    //    in that case the root r is a failed lineral (or one of its descendents!)
    std::list<var_t> failing_v;
    vec<bool> marked(no_v, false);
    vec<var_t> root(no_v, L.size());
    std::stack<var_t> dfs_q;
    for(const auto& r : roots) {
        dfs_q.push(r);
        while(!dfs_q.empty()) {
            const auto v = dfs_q.top();
            dfs_q.pop();
            if(marked[IL[v]]) continue;
            //discover v
            if(marked[IL[SIGMA(v)]] && root[IL[SIGMA(v)]] == r) {
                failing_v.emplace_back( v );
            }
            marked[IL[v]] = true;
            root[IL[v]] = r;
            for(const auto& n : get_out_neighbour_range(v)) {
                if(!marked[IL[n]]) dfs_q.push(n);
            }
        }
    }
    // (3) for each pair (n,v) in flits we get at least one failed lineral; 
    vec<bool> marked_sigma(no_v,false);

    //for each r in failing_v run a ascending bfs from r and SIGMA(r); all elements in the intersection are failed linerals!
    vec<lineral> f_xlits;
    for(const auto& r : failing_v) {
        for(var_t i=0; i<marked.size(); ++i) {
            marked[i] = false;
            marked_sigma[i] = false;
        }
        //dfs from r
        dfs_q.push(r);
        while(!dfs_q.empty()) {
            const auto v = dfs_q.top();
            dfs_q.pop();
            if(marked[IL[v]]) continue;
            //discover v
            marked[IL[v]] = true;
            for(const auto& n : get_in_neighbour_range(v)) if(!marked[IL[n]]) dfs_q.push(n);
        }
        //dfs from SIGMA(r)
        dfs_q.push(SIGMA(r));
        while(!dfs_q.empty()) {
            const auto v = dfs_q.top();
            dfs_q.pop();
            if(marked_sigma[IL[v]]) continue;
            //discover v
            marked_sigma[IL[v]] = true;
            if(marked[IL[v]]) {
                //FAILED LINERAL FOUND!
                f_xlits.emplace_back( std::move( vl.Vxlit(v) ) );
                f_xlits.back().add_one();
            #ifndef FULL_REDUCTION
                f_xlits.back().reduce(assignments);
            #endif
            }
            for(const auto& n : get_in_neighbour_range(v)) if(!marked_sigma[IL[n]]) dfs_q.push(n);
        }
    }
    return LinEqs( std::move(f_xlits) );
}

LinEqs impl_graph::fls_trivial() const {
    //(1) compute roots
    auto roots = get_roots();
    //(2) from every root r perform dfs; on discovery of v mark it with r and check if SIGMA(v) was already marked with r;
    //    in that case the root r is a failed lineral (or one of its descendents!)
    std::list<var_t> failing_v;
    vec<var_t> mark_root(no_v);
    vec<bool>marked(no_v, false);
    std::stack<var_t> dfs_q;
    for(const auto& r : roots) {
        dfs_q.push(r);
        while(!dfs_q.empty()) {
            const auto v = dfs_q.top();
            dfs_q.pop();
            if(marked[IL[v]]) continue;
            //discover v
            if(marked[IL[SIGMA(v)]] && mark_root[IL[SIGMA(v)]] == r) {
                failing_v.emplace_back( v );
            }
            marked[IL[v]] = true;
            mark_root[IL[v]] = r;
            for(const auto& n : get_out_neighbour_range(v)) {
                if(!marked[IL[n]]) dfs_q.push(n);
            }
        }
    }
    // (3) for each pair (n,v) in flits we get at least one failed lineral; 
    vec<bool> marked_sigma(no_v,false);

    //for each r in failing_v run a ascending bfs from r and SIGMA(r); all elements in the intersection are failed linerals!
    vec<lineral> f_xlits;
    for(const auto& r : failing_v) {
        for(var_t i=0; i<marked.size(); ++i) {
            marked[i] = false;
            marked_sigma[i] = false;
        }
        //dfs from r
        dfs_q.push(r);
        while(!dfs_q.empty()) {
            const auto v = dfs_q.top();
            dfs_q.pop();
            if(marked[IL[v]]) continue;
            //discover v
            marked[IL[v]] = true;
            for(const auto& n : get_in_neighbour_range(v)) if(!marked[IL[n]]) dfs_q.push(n);
        }
        //dfs from SIGMA(r)
        dfs_q.push(SIGMA(r));
        while(!dfs_q.empty()) {
            const auto v = dfs_q.top();
            dfs_q.pop();
            if(marked_sigma[IL[v]]) continue;
            //discover v
            marked_sigma[IL[v]] = true;
            if(marked[IL[v]]) {
                //FAILED lineral FOUND!
                f_xlits.emplace_back( std::move( vl.Vxlit(v) ) );
                f_xlits.back().add_one();
            #ifndef FULL_REDUCTION
                f_xlits.back().reduce(assignments);
            #endif
            }
            for(const auto& n : get_in_neighbour_range(v)) if(!marked_sigma[IL[n]]) dfs_q.push(n);
        }
    }
    return LinEqs( std::move(f_xlits) );
    /*
    vec<lineral> new_xlits;
    vec< robin_hood::unordered_flat_set<var_t> > reaches(no_v);
    //vec< std::set<var_t> > reaches(no_v);
    vec<bool> failed(no_v, false);

    const auto TO = get_TO();
    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
        const var_t v = *v_it;
        if(failed[IL[v]]) continue;
        reaches[ IL[v] ].insert( v );
        for (const auto &w : get_out_neighbour_range(v)) {
            if(failed[IL[w]]) failed[IL[v]]=true;
            if(failed[IL[v]]) break;
            assert(!failed[IL[v]]);
            //combine reaches
            for(const auto &l : reaches[ IL[w] ]) {
                if(reaches[IL[v]].contains(SIGMA(l))) {
                    //v is a failed lineral!
                    failed[IL[v]] = true;
                } else {
                    reaches[IL[v]].insert(l);
                }
            }
        }
    }

    for(const auto& v : get_v_range()) {
        if(failed[IL[v]]) {
            new_xlits.emplace_back( std::move( vl.Vxlit(v) ) );
            new_xlits.back().add_one();
        }
    }

    return std::move( LinEqs(new_xlits) );
    */
};

LinEqs impl_graph::fls_full() const {
    vec<lineral> new_xlits;
    //compute topological ordering of graph
    const auto TO = get_TO();
    
    vec<LinEqs> D(no_v);
    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
        var_t v = *v_it;
        lineral f = vl.Vxlit(v);
        D[IL[v]] += LinEqs(f);
        for (const auto &w : get_in_neighbour_range(v)) D[IL[w]] += D[IL[v]];
        if(!D[ IL[v] ].is_consistent() ) new_xlits.emplace_back( std::move( f.add_one() ) );
    }

    //now deduce <D[f]> cap <D[f+1]> for all f
    vec<bool> marked(no_v, false);
    for(const auto &v : get_v_range()) {
        if(marked[IL[SIGMA(v)]]) continue;
        marked[IL[v]] = true;
        const auto intVS = intersect(D[IL[v]], D[IL[SIGMA(v)]]);
        if(intVS.size()<1) continue;
        VERB(80, "c GFLS derived "+std::to_string(intVS.size())+" new eqs");
        new_xlits.insert(new_xlits.end(), intVS.begin(), intVS.end());
    }

#ifndef FULL_REDUCTION
    std::for_each(new_xlits.begin(), new_xlits.end(), [&](lineral& l){ l.reduce(assignments); });
#endif

    return LinEqs( std::move(new_xlits) );
};

LinEqs impl_graph::fls_full_implied() {
    vec<lineral> new_xlits;
    //compute topological ordering of graph
    const auto TO = get_TO();
    
    vec<LinEqs> D(no_v);
    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
        var_t v = *v_it;
        lineral f = vl.Vxlit(v);
        D[IL[v]] = implied_xlits(f);
    }

    //now deduce <D[f]> cap <D[f+1]> for all f
    vec<bool> marked(no_v, false);
    for(const auto &v : get_v_range()) {
        if(marked[IL[SIGMA(v)]]) continue;
        marked[IL[v]] = true;
        const auto intVS = intersect(D[IL[v]], D[IL[SIGMA(v)]]);
        if(intVS.size()<1) continue;
        VERB(80, "c GFLS derived "+std::to_string(intVS.size())+" new eqs");
        new_xlits.insert(new_xlits.end(), intVS.begin(), intVS.end());
    }

#ifndef FULL_REDUCTION
    std::for_each(new_xlits.begin(), new_xlits.end(), [&](lineral& l){ l.reduce(assignments); });
#endif

    return LinEqs( std::move(new_xlits) );
};

//LinEqs impl_graph::fls_full() {
//    vec<lineral> new_xlits;
//    vec<LinEqs> prev_xsys(no_v);
//    //compute topological ordering of graph
//    const auto TO = get_TO();
//    
//    for (auto v_it = TO.rbegin(); v_it != TO.rend(); ++v_it) {
//        var_t v = *v_it;
//        if(get_out_degree(v) == 0) {
//            prev_xsys[ IL[v] ] = LinEqs( std::move( vl.Vxlit(v) ) );
//        } else {
//            prev_xsys[ IL[v] ] = LinEqs( std::move( vl.Vxlit(v) ) );
//            for (const auto &w : get_out_neighbour_range(v)) prev_xsys[IL[v]] += prev_xsys[IL[w]];
//        }
//        if(!prev_xsys[ IL[v] ].is_consistent() ) new_xlits.emplace_back( std::move( vl.Vxlit(v).add_one() ) );
//    }
//
//    return std::move( LinEqs( std::move(new_xlits) ) );
//};

void impl_graph::bump_score(const LinEqs& new_xsys) {
    for (const auto &[lt,_] : new_xsys.get_pivot_poly_idx()) {
        assert(lt >= 0 && lt < activity_score.size());
        activity_score[ lt ] += bump;
    }
};

void impl_graph::decay_score() {
    for (auto &s : activity_score) {
        s = ceil(s*decay); //TODO make more efficient?! (since result of mult first is float and then is casted back to unsigned int...)
    }
};


//pre-process

void impl_graph::preprocess() {
    //set update/fls/decH funcs
    upd_t upd;
    switch (opt.upd) {
        case upd_alg::ts:
            upd = &impl_graph::update_graph; break;
        case upd_alg::hf:
            upd = &impl_graph::update_graph_hash_fight; break;
        case upd_alg::hfd:
            upd = &impl_graph::update_graph_hash_fight_dev; break;
        case upd_alg::par:
            upd = &impl_graph::update_graph_par; break;
        default:
            upd = &impl_graph::update_graph_hash_fight; break;
    }

    fls_t fls;
    switch (opt.pp) {
        case preproc::no:
            VERB(40, "c preprocess 'no'")
            //instant return!
            return;
        case preproc::scc:
            VERB(40, "c preprocess 'scc'")
            fls = &impl_graph::fls_no; break;
        case preproc::fls_scc:
            VERB(40, "c preprocess 'fls_scc'")
            fls = &impl_graph::fls_full; break;
        case preproc::fls_scc_ee:
            VERB(40, "c preprocess 'fls_scc_ee'")
            fls = &impl_graph::fls_full; break;
    }
    //perform crGCP with selected fls
    crGCP_no_schedule(s, upd, fls);
};

//make global var xsyses s.t. the do not need to be destroyed and constructed on each call to crGCP
LinEqs new_, scc, fls, upd;

void impl_graph::crGCP(stats& s, const upd_t upd_graph, const fls_t fls_alg, const bool scheduled_fls ) {
    if(!linsys.is_consistent()) return;
    ++s.no_crGCP;

    bool repeat = true;

    while(repeat) {
        repeat = false;

        do {
            upd = (this->*upd_graph)( s, linsys );
            VERB(40, "c       |---> deduced " << std::to_string(upd.size()) << " new eqs (upd)")
            if(upd.size() > 0) {
                s.new_px_upd += upd.size();
                add_new_xsys( upd );
                repeat = true;
                //if(opt.score == sc::active) {  bump_score(upd); }
            }
        } while (upd.size() > 0 && upd.is_consistent());
        if(!upd.is_consistent()) return;
        
        //in-processing!
        scc = scc_analysis();
        if(scc.size() > 0) {
            s.new_px_scc += scc.size();
            add_new_xsys( scc );
            repeat = true;
            //if(opt.score == sc::active) {  bump_score(scc); }
        }
        VERB(40, "c       |---> deduced " << std::to_string(scc.size()) << " new eqs (scc)")
        if(!scc.is_consistent()) return;
        
        if(scc.size() == 0 && (!scheduled_fls || (s.no_crGCP % opt.fls_s == 0))) {
            fls = (this->*fls_alg)();
            if(fls.size() > 0) {
                s.new_px_fls += fls.size();
                add_new_xsys( fls );
                repeat = true;
                //if(opt.score == sc::active) {  bump_score(fls); }
            }
            VERB(40, "c       |---> deduced " << std::to_string(fls.size()) << " new eqs (fls)")
            if(!fls.is_consistent()) return;
        }
    }
    assert( !linsys.is_consistent() || is_DAG() );
};



//implementation of graph-based dpll-solver
stats impl_graph::dpll_solve(stats& s) {
    VERB(25, "c dpll-solving start")
    //return UNSAT if linsys has no solution
    if(!linsys.is_consistent()) {
        s.sat = false;
        s.finished = true;
        return s;
    }

    //current decision lvl
    var_t dl = 0;

    //set update/fls/decH funcs
    upd_t upd;
    switch (opt.upd) {
        case upd_alg::ts:
            upd = &impl_graph::update_graph; break;
        case upd_alg::hf:
            upd = &impl_graph::update_graph_hash_fight; break;
        case upd_alg::hfd:
            upd = &impl_graph::update_graph_hash_fight_dev; break;
        case upd_alg::par:
            upd = &impl_graph::update_graph_par; break;
        default:
            upd = &impl_graph::update_graph_hash_fight; break;
    }
    
    fls_t fls;
    switch (opt.fls) {
        case fls_alg::no:
            fls = &impl_graph::fls_no; break;
        case fls_alg::trivial:
            fls = &impl_graph::fls_trivial; break;
        case fls_alg::trivial_cc:
            fls = &impl_graph::fls_trivial_cc; break;
        case fls_alg::full:
            fls = &impl_graph::fls_full; break;
        default:
            fls = &impl_graph::fls_no; break;
    }

    dec_heu_t decH;
    switch (opt.dh) {
        case dec_heu::fv:
            decH = &impl_graph::first_vert; break;
        case dec_heu::mp:
            decH = (opt.score == sc::active) ? &impl_graph::max_score_path : &impl_graph::max_path;
            break;
        case dec_heu::mr:
            decH = &impl_graph::max_reach; break;
        case dec_heu::mbn:
            decH = &impl_graph::max_bottleneck; break;
        case dec_heu::lex:
            decH = &impl_graph::lex; break;
        default:
            decH = &impl_graph::max_path; break;
    }

    //stack for LinEqs that store alternative dec
    auto backtrack_xsys = std::stack<LinEqs>();
    LinEqs new_xsys = LinEqs();

    //update graph -- before making decisions!
    VERB(45, graph_stats());
    crGCP(s, upd, fls);
    VERB(45, graph_stats());

    while( no_e > 0 || !linsys.is_consistent() ) {
        if( s.cancelled.load() ) {
            VERB(10, "c cancelled");
            return s;
        }
        //make decision / backtrack
        if( !linsys.is_consistent() ) {
        backtrack:
            VERB(25, "c " << std::to_string( dl ) << " : " << "conflict --> backtrack!")
            ++s.no_confl;
            //conflict!
            if( dl == 0) {
                //return UNSAT
                s.finished = true;
                s.sat = false;
                return s;
            }
            //decay + bump scores
            if(opt.score == sc::active) { //only update scores if used!
                for(const auto& sys : xsys_stack.back() ){
                    bump_score(sys);
                }
                decay_score();
            }
            //backtrack
            --dl;
            vl.backtrack( std::move(vl_stack.top()), dl );
            vl_stack.pop();
            //revert assignments
        #ifndef FULL_REDUCTION
            for(const auto& L : Lsys) {
                for(const auto& [lt,_] : L.get_pivot_poly_idx()) {
                    assignments[lt].reset();
                }
            }
        #endif
            xsys_stack.pop_back();
            backtrack( std::move(graph_stack.top()) );
            assert( assert_data_structs() );
            graph_stack.pop();
            //add forced alt decision
            add_new_xsys( std::move(backtrack_xsys.top()) );
            //Lsys.emplace_back( backtrack_xsys.top() );
            backtrack_xsys.pop();
            //check whether backtracking worked!
        } else {
            //make new decision!
            //use decisions heuristic to find next decision!
            std::pair<LinEqs,LinEqs> dec = (this->*decH)();
            //as long as at least one of the two decisions is inconsistent, instead of taking decision propagate the consistent decision, if there is one, and take another 'decision'; otherwise immediately backtrack
            while (!dec.first.is_consistent() || !dec.second.is_consistent()) {
                if(dec.first.is_consistent()) add_new_xsys(dec.first); //Lsys.emplace_back( dec.first );
                else if(dec.second.is_consistent()) add_new_xsys(dec.second); //Lsys.emplace_back( dec.second );
                else /* both decisions are inconsistent */ goto backtrack;
                VERB(40, "c       |---> deduced " << std::to_string(linsys.size()) << " new eqs (decH)")
                //update graph
                crGCP(s, upd, fls);
                //if now inconsistent; we need to backtrack; for this take decision '1,1'; otherwise take a different decision
                if(!linsys.is_consistent()) {
                    //we need to backtrack last decision
                    goto backtrack;
                } else {
                    if(no_e==0) goto comp_sol;
                    //take another decision
                    dec = (this->*decH)();
                }
            }

            ++dl;
            ++s.no_dec;
            //save state
            graph_stack.push( std::move(get_state()) );
             //duplicate top of vl_stack
            vl_stack.push( std::move(vl.get_state()) );

            VERB(25, "c " << std::to_string( dl ) << " : " << "decision " << std::to_string(s.no_dec) << " : " << std::to_string(dec.first.size()) << " or " << std::to_string(dec.second.size()) << " eqs")
            VERB(50, "c " << std::to_string( dl ) << " : " << "decision " << std::to_string(s.no_dec) << " namely [" << dec.first.to_str() << "] or [" << dec.second.to_str() << "]")
            xsys_stack.emplace_back( std::list<LinEqs>() );
            add_new_xsys( std::move(dec.first) );
            backtrack_xsys.emplace( std::move(dec.second) );
        }
        
        //update graph
        crGCP(s, upd, fls);

        VERB(45, graph_stats());
        assert((var_t) graph_stack.size() == dl+1);
        assert((var_t) vl_stack.size()  == dl+1);
        assert((var_t) xsys_stack.size()  == dl+1);
    }

    comp_sol:
    //solution can be deduced from xsys_stack!
    s.sol = vec<bool>(opt.num_vars);
    while(!xsys_stack.empty()) {
        const auto l_xsys = xsys_stack.back();
        //work through list in reverse order (to solve last linsys first!)
        for (auto sys = l_xsys.rbegin(); sys != l_xsys.rend(); ++sys) sys->solve( s.sol );
        xsys_stack.pop_back();
    }

    s.sat = true;
    s.finished = true;

    return s;
};


//note: should be re-written as it uses strings to ensure no clauses are stored multiple times...
vec< vec<lineral> > impl_graph::to_xcls() const {
    vec<vec<lineral>> xclss;
    auto xclss_str = robin_hood::unordered_flat_set<std::string>();
    //go through edges
    for(const auto& v : get_v_range()) {
        const lineral fp1 = vl.Vxlit(v).add_one();
        for(const auto& n : get_out_neighbour_range(v)) {
            const lineral g = vl.Vxlit(n);
            //append (f+1), g to xclss if it is not yet present!
            const auto fp1_s = fp1.to_str();
            const auto g_s = g.to_str();
            const auto fpgp1_s = (fp1+g).to_str();
            //if( xclss_str.contains(fp1_s+" "+g_s) || xclss_str.contains(g_s+" "+fp1_s) || xclss_str.contains(fpgp1_s+" "+fp1_s) || xclss_str.contains(fp1_s+" "+fpgp1_s) ) {
            if( xclss_str.contains(fp1_s+" "+g_s) || xclss_str.contains(g_s+" "+fp1_s)) {
                continue;
            }
            xclss.emplace_back( std::move( vec<lineral>({fp1,g}) ) );
            xclss_str.emplace( fp1_s+" "+g_s );
        }
    }
    //add linear polys
    for(const auto & sys : Lsys) {
        for(const auto& l : sys.get_linerals()) {
            xclss.emplace_back( std::move( vec<lineral>({l})) );
        }
    }

    return xclss;
};

std::string impl_graph::to_xnf_string() const {
    auto xclss_str = std::set<std::string>();
    var_t n_cls = 0;
    //go through edges
    for(const auto& v : get_v_range()) {
        const lineral fp1 = vl.Vxlit(v).add_one();
        for(const auto& n : get_out_neighbour_range(v)) {
            const lineral g = vl.Vxlit(n);
            //append (f+1), g to xclss if it is not yet present!
            const auto fp1_s = fp1.to_xnf_str();
            const auto g_s = g.to_xnf_str();
            const auto fpgp1_s = (fp1+g).to_xnf_str();
            //if( xclss_str.contains(fp1_s+" "+g_s) || xclss_str.contains(g_s+" "+fp1_s) || xclss_str.contains(fpgp1_s+" "+fp1_s) || xclss_str.contains(fp1_s+" "+fpgp1_s) ) {
            if( xclss_str.contains(fp1_s+" "+g_s) || xclss_str.contains(g_s+" "+fp1_s)) {
                continue;
            }
            xclss_str.emplace( fp1_s+" "+g_s );
            ++n_cls;
        }
    }
    //add linear polys
    for(const auto & sys : Lsys) {
        for(const auto& l : sys.get_linerals()) {
            xclss_str.emplace( l.to_xnf_str() );
        }
    }
    //convert to one big string
    std::string str = "p xnf "+std::to_string(get_const_opts()->num_vars)+" "+std::to_string(n_cls)+"\n";
    for(const auto &cls : xclss_str) {
        str += cls + " 0\n";
    }
    return str;
};


//overwrite to_str() func
std::string impl_graph::to_str() const noexcept {
    std::map<var_t, vec<var_t> > edges;
    for (var_t c_idx = 0; c_idx < no_v; ++c_idx) {
        //get color at c_idx
        const var_t c = L[c_idx];
        //add all out-neighbors to edges[c_idx]
        edges[c] = get_out_neighbour_vector( c );
    #ifdef USE_LHGR
        assert( edges[c].size() == CD_out[c] );
    #endif
        //sort color c_idx
        std::sort( std::execution::par, edges[c].begin(), edges[c].end() );
    }

    //generate string of edges with lexicographic ordering!
    vec< std::string > str_edges(0);
    for (const auto &[src,v] : edges) {
        vec<std::string> out_edges_str( v.size() );
        //construct strings!
        auto to_str = [&](const var_t dst) -> std::string {return "(" + vl.Vxlit(src).to_str() + "," + vl.Vxlit(dst).to_str() + ")";};
        std::transform(v.begin(), v.end(), out_edges_str.begin(), to_str);
        std::sort(std::execution::par, out_edges_str.begin(), out_edges_str.end());

        std::stringstream ss;
        std::copy(out_edges_str.begin(), out_edges_str.end(), std::ostream_iterator<std::string>(ss, " "));
        std::string out_edges = ss.str();
        if (!out_edges.empty()) {
            out_edges.resize(out_edges.length() - 1); // trim trailing space
            str_edges.emplace_back(out_edges);
        }
    }

    std::stringstream ss;
    std::sort(std::execution::par, str_edges.begin(), str_edges.end());
    std::copy(str_edges.begin(), str_edges.end(), std::ostream_iterator<std::string>(ss, "; "));
    std::string result = ss.str();
    if (!result.empty()) {
        result.resize(result.length() - 2); // trim trailing semicolon+space
        result += "\n";
    } else {
        result = "";
    }

    //append lin-eqs
    for(const auto & sys : Lsys) {
        result += sys.to_str() + "\n";
    }
    if(!xsys_stack.empty()) result.resize(result.length()-1);
    
    return result;
} 


std::string impl_graph::to_str_() const { return graph::to_str(); };

bool impl_graph::assert_data_structs() const noexcept {
    //make sure V and Vxlit have proper keys!
    for ([[maybe_unused]] auto &v : get_v_range()) {
        //std::cout << "check label of vert v=" <<  std::to_string(v) << "..." << std::endl;
        //if(v != vl.V( vl.Vxlit(v) )) {
        //    std::cout << "V and Vxlit do not agree on vert v=" << std::to_string(v) << " with label " << vl.Vxlit(v).to_str() << std::endl;
        //    std::cout << "Vxlit(v)=" << vl.Vxlit(v).to_str() << std::endl;
        //    lineral lit_v = vl.Vxlit(v);
        //    std::cout << "V( Vxlit(v) )=" << V( lit_v ) << std::endl;
        //}
        assert(v == vl.V( vl.Vxlit(v) ));
        //std::cout << "retrieved lit=" << Vxlit(v).to_str() << std::endl;
    }
    assert( graph::assert_data_structs() );
    return true;
};