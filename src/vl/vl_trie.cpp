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

#include "vl_trie.hpp"

#include <algorithm>

n_t curr_node, new_node_idx;

void vl_trie::backtrack(trie_repr&& r, [[maybe_unused]] const var_t dl) noexcept {
    v_node = std::move(r.v_node);
    assigned_vert.clear();
    for(const auto [v_idx,n_idx] : v_node) assigned_vert.insert( {n_idx,v_idx} );
    num_vs = std::move(r.num_vs);
    //prune trie
    prune(dl);
}

void vl_trie::prune(const var_t dl) noexcept {
    //std::cout << "pruning start! (trie-size = " << std::to_string(get_num_nodes()) << ", dl = " << std::to_string(dl) << ")" << std::endl;
    //loop through all nodes added in higher dls and remove them
    while((var_t) nodes_in_dl.size() > dl+1) {
        for (auto &n_idx : nodes_in_dl.top()) {
            remove_node( n_idx );
        }
        nodes_in_dl.pop();
    }
    //std::cout << "pruning done! (trie-size = " << std::to_string(get_num_nodes()) << ")" << std::endl;
}


const trie_insert_return_type vl_trie::insert(const var_t v, const lineral& lit, const var_t dl) {
    if( v_node.contains(v) ) return trie_insert_return_type(false, false, get_vert(v_node[v]));

    curr_node = ROOT;
    bool node_added = false; //as soon as one node was added, we know that we never have to check nodes[curr_node].children again!
    //ignores constant!
    for (auto it = lit.get_idxs_().rbegin(); it != lit.get_idxs_().rend(); ++it) {
        var_t ind = *it;
        const auto search = node_added ? nodes[curr_node].children.end() : nodes[curr_node].children.find(ind);
        if(search==nodes[curr_node].children.end()) {
            //add new node!
            new_node_idx = add_node(curr_node, ind, dl);
            curr_node = new_node_idx;
            node_added = true;
        } else {
            curr_node = search->second;
        }
    }

    if(lit.has_constant()) {
        //lit has constant
        if(has_assigned_vert(curr_node)) {
            return trie_insert_return_type(false, true, get_vert(curr_node));
        } else {
            //go down one more lvl!
            const auto search = nodes[curr_node].children.find( 0 );
            if(search==nodes[curr_node].children.end()) {
                //add new node!
                new_node_idx = add_node(curr_node, 0, dl);
                curr_node = new_node_idx;
            } else {
                curr_node = search->second;
            }
            //curr_node is at lit!
            if(has_assigned_vert(curr_node)){
                return trie_insert_return_type(false, false, get_vert(curr_node));
            } else {
                //assign node!
                assign_vert(curr_node, v);
                return trie_insert_return_type(true, false, get_vert(curr_node));
            }
        }
    } else {
        //lit has no constant -- i.e. we are at lit
        if(has_assigned_vert(curr_node)) {
            return trie_insert_return_type(false, false, get_vert(curr_node));
        } else {
            //check if lit with constant exists with label!
            const auto search = nodes[curr_node].children.find( 0 );
            if(search!=nodes[curr_node].children.end() && has_assigned_vert(nodes[curr_node].children.at(0)) ) {
                return trie_insert_return_type(false, true, get_vert(nodes[curr_node].children[0]) );
            } else {
                //assign curr_node to vert v!
                assign_vert(curr_node, v);
                return trie_insert_return_type(true, false, get_vert(curr_node));
            }
        }
    }
};

bool vl_trie::erase(const var_t v) {
    assert( contains(v) );
    //note: only removes assigned vert, does not change underlying data struct!
    if(contains(v)) {
        assigned_vert.erase( v_node.at(v) );
        v_node.erase(v);
        num_vs--;
        assert(!contains(v));
        return true;
    } else {
        return false;
    }
};


std::pair<var_t,bool> vl_trie::update(const var_t v, const lineral& l, const var_t dl) {
    assert(contains(v));
    [[maybe_unused]] bool erased = erase(v);
    assert(erased);
    auto ins = insert(v, l, dl);
    //prune(dl);
    return { ins.inserted ? v : ins.vert, ins.found_plus_one };
};


lineral vl_trie::operator[](const var_t v) const {
    vec<var_t> idxs(0);
    idxs.reserve( get_num_nodes()/num_vs + num_vars/10 );
    curr_node = v_node.at(v);
    while(curr_node!=ROOT) {
        idxs.push_back( nodes[curr_node].label );
        curr_node = nodes[curr_node].parent;
    }
    return lineral(std::move(idxs), true);
};

lineral vl_trie::at(const var_t v) const {
    //if(v_node.at(v)==ROOT && assigned_vert.at(v_node.at(v)) != v) throw std::out_of_range("Label of vertex " + std::to_string(v) + " not found in trie.");
    return lineral( vec<var_t>(begin(v),end()), true);
};

var_t vl_trie::operator[](const lineral& lit) const {
    curr_node = ROOT;
    for (auto it = lit.get_idxs_().rbegin(); it != lit.get_idxs_().rend(); ++it) {
        var_t ind = *it;
        const auto search = nodes[curr_node].children.find( ind );
        if(search != nodes[curr_node].children.end()) {
            curr_node = search->second;
        } else {
            return 0;
        };
    }
    if(lit.has_constant()) {
        const auto search = nodes[curr_node].children.find( 0 );
        if(search != nodes[curr_node].children.end()) {
            curr_node = search->second;
        } else {
            return 0;
        }
    }
    //check if final curr_node has a vertex
    return has_assigned_vert(curr_node) ? assigned_vert.at(curr_node) : 0;
};

var_t vl_trie::at(const lineral& lit) const {
    curr_node = ROOT;
    for (auto it = lit.get_idxs_().rbegin(); it != lit.get_idxs_().rend(); ++it) {
        var_t ind = *it;
        const auto search = nodes[curr_node].children.find( ind );
        if(search != nodes[curr_node].children.end()) {
            curr_node = search->second;
        } else {
            throw std::out_of_range("Vertex of label " + lit.to_str() + " not found in trie.");
        };
    }
    if(lit.has_constant()) {
        const auto search = nodes[curr_node].children.find( 0 );
        if(search != nodes[curr_node].children.end()) {
            curr_node = search->second;
        } else {
            throw std::out_of_range("Vertex of label " + lit.to_str() + " not found in trie.");
        }
    }

    const auto search = assigned_vert.find( curr_node );
    if(search == assigned_vert.end() ) throw std::out_of_range("Vertex of label " + lit.to_str() + " not found in trie.");
    return search->second;
};

std::pair<var_t,bool> vl_trie::at_(const lineral& lit) const {
    //iter down the trie!
    curr_node = ROOT;
    for (auto it = lit.get_idxs_().rbegin(); it != lit.get_idxs_().rend(); ++it) {
        var_t ind = *it;
        const auto search = nodes[curr_node].children.find( ind );
        //go down one more lvl -- if possible!
        if(search != nodes[curr_node].children.end()) {
            curr_node = search->second;
        } else {
            throw std::out_of_range("Vertex of label " + lit.to_str() + " not found in trie.");
        };
    }
    //we ignored constant so far!

    //check if there is an assigned vert!
    const auto search = assigned_vert.find( curr_node );
    if(search != assigned_vert.end() ) {
        return {search->second, !lit.has_constant()};
    } else {
        //lit without constant is not contained, check if lit with constant is there!
        //go down on more step!
        const auto search_ = nodes[curr_node].children.find( 0 );
        if(search_ != nodes[curr_node].children.end()) {
            curr_node = search_->second;
        } else {
            throw std::out_of_range("Vertex of label " + lit.to_str() + " not found in trie.");
        };

        //now curr_node points to lit with constant; check if it has an assigned vert!
        const auto search2 = assigned_vert.find( curr_node );
        if(search2 == assigned_vert.end() ) throw std::out_of_range("Vertex of label " + lit.to_str() + " not found in trie.");
        return {search2->second, lit.has_constant()};
    }
};

bool vl_trie::contains(const lineral& lit) const {
    //if( lit.has_constant() ) return false;
    curr_node = ROOT;
    for (auto it = lit.get_idxs_().rbegin(); it != lit.get_idxs_().rend(); ++it) {
        var_t ind = *it;
        const auto search = nodes[curr_node].children.find( ind );
        if(search != nodes[curr_node].children.end()) {
            curr_node = search->second;
        } else {
            return 0;
        };
    }
    if(lit.has_constant()) {
        const auto search = nodes[curr_node].children.find( 0 );
        if(search != nodes[curr_node].children.end()) {
            curr_node = search->second;
        } else {
            return 0;
        }
    }
    return has_assigned_vert(curr_node);
};
    
std::string vl_trie::to_str() const {
    std::string str = "";
    for(const auto& [v,n_idx] : v_node) {
        str += "(" + std::to_string(v) + "," + ( this->operator[](v) ).to_str() + ") ";
    }
    if(str.size() > 0) str.pop_back();
    return str;
};
    
vec<var_t> diff(0);
lineral vl_trie::sum(const var_t lhs, const var_t rhs) const {
  diff.clear();
  diff.reserve( get_num_nodes()/num_vs + num_vars/10 );
  std::set_symmetric_difference(std::execution::par, begin(lhs), end(), begin(rhs), end(), std::back_inserter(diff));
  //NOTE back_insterter might lead to repeated reallocations!
  return lineral(std::move(diff), true); //call ctor that does NOT sort diff
}