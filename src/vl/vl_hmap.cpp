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


#include "vl_hmap.hpp"


const vl_hmap_insert_return_type vl_hmap::insert(const var_t v, lineral&& lit, [[maybe_unused]] const var_t dl) {
    assert(!lit.has_constant());
    auto inserted = xl_to_v_stack.top().emplace( lit, v );
    if(inserted.second) {
        [[maybe_unused]] auto inserted2 = v_to_xl_stack.top().emplace( v, std::move(lit) ); 
        assert(inserted2.second);
    }
    return vl_hmap_insert_return_type(inserted.second, inserted.second ? v : (inserted.first)->second);
};

bool vl_hmap::erase(const var_t v) {
    bool erased = xl_to_v_stack.top().erase( v_to_xl_stack.top().at(v) );
    bool erased2 = v_to_xl_stack.top().erase( v );
    assert(erased && erased2);
    return erased && erased2;
};

std::pair<var_t,bool> vl_hmap::update(const var_t v, lineral&& l, const var_t dl) {
    [[maybe_unused]] bool erased = erase( v );
    assert(erased);
    bool found_plus_one = false;
    if(l.has_constant()) {
        l.add_one();
        found_plus_one = true;
    }
    auto ins = insert(v, std::move(l), dl);
    return { ins.inserted ? v : ins.vert, found_plus_one};
};


std::string vl_hmap::to_str() const {
    std::string str = "";
    for(const auto& [v,l] : v_to_xl_stack.top()) {
        str += "(" + std::to_string(v) + "," + l.to_str() + ") ";
    }
    if(str.size() > 0) str.pop_back();
    return str;
};