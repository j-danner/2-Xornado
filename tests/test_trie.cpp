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

//file to test implementation of trie
#include "../src/vl/vl_trie.hpp"

#include "../src/impl_graph.hpp"
#include "../src/solve.hpp"

#include <catch2/catch_all.hpp>


TEST_CASE( "trie creation, insert", "[trie]" ) {
    vl_trie tr = vl_trie(4,4);

    lineral l1 = lineral({0,1,2,3});
    var_t v1 = 0;

    lineral l2 = lineral({2,3});
    var_t v2 = 1;
    
    lineral l3 = lineral({0,1,2});
    var_t v3 = 2;

    auto ins1 = tr.insert(v1, l1, 0);
    CHECK( tr.size()==1 );
    CHECK( tr.get_num_nodes()==5 );
    CHECK( tr.to_str() == "(0,x1+x2+x3+1)" );

    auto ins2 = tr.insert(v2, l2, 0);
    CHECK( tr.size()==2 );
    CHECK( tr.get_num_nodes()==5 );
    CHECK( tr.to_str() == "(0,x1+x2+x3+1) (1,x2+x3)" );

    auto ins3 = tr.insert(v3, l3, 0);
    CHECK( tr.size()==3 );
    CHECK( tr.get_num_nodes()==8 );
    CHECK( tr.to_str() == "(0,x1+x2+x3+1) (1,x2+x3) (2,x1+x2+1)" );

    CHECK( ins1.inserted );
    CHECK( ins2.inserted );
    CHECK( ins3.inserted );

    CHECK( ins1.vert == v1 );
    CHECK( ins2.vert == v2 );
    CHECK( ins3.vert == v3 );

    lineral tmp = tr[v1];

    CHECK( tr[v1] == l1 );
    CHECK( tr[v2] == l2 );
    CHECK( tr[v3] == l3 );
    
    CHECK( tr[l1] == v1 );
    CHECK( tr[l2] == v2 );
    CHECK( tr[l3] == v3 );

    CHECK( tr.at(v1) == l1 );
    CHECK( tr.at(v2) == l2 );
    CHECK( tr.at(v3) == l3 );
    
    CHECK( tr.at(l1) == v1 );
    CHECK( tr.at(l2) == v2 );
    CHECK( tr.at(l3) == v3 );

    CHECK( tr.contains(l1) );
    CHECK( tr.contains(l2) );
    CHECK( tr.contains(l3) );
    
    CHECK( tr.contains(v1) );
    CHECK( tr.contains(v2) );
    CHECK( tr.contains(v3) );

    CHECK( !tr.contains(l1+l2) );
    CHECK( !tr.contains(l2+l3) );
    CHECK( !tr.contains(3) );
    CHECK( !tr.contains(4) );

    {
        lineral n = lineral({1});
        auto ins_fail1 = tr.insert(v1, n, 0);

        CHECK( !ins_fail1.inserted );
        CHECK( ins_fail1.vert == v1 );
        CHECK( tr.size()==3 );
        CHECK( tr.get_num_nodes()==8 );
    }
    
    {
        var_t n = 3;
        auto ins_fail2 = tr.insert(n, l1, 0);

        CHECK( !ins_fail2.inserted );
        CHECK( ins_fail2.vert == v1 );
        CHECK( tr.size()==3 );
        CHECK( tr.get_num_nodes()==8 );
    }

    //CHECK_THROWS( tr.erase(4) );

    CHECK( tr.erase(1) );
    CHECK( tr.to_str() == "(0,x1+x2+x3+1) (2,x1+x2+1)" );
    CHECK( tr.erase(0) );
    CHECK( tr.to_str() == "(2,x1+x2+1)" );
    CHECK( tr.erase(2) );
    CHECK( tr.to_str() == "" );

    lineral zero({},5);
    tr.insert(1,zero, 0);
    CHECK(tr[zero] == 1);
    CHECK(tr[1] == zero);
    CHECK(tr.at(zero) == 1);
    CHECK(tr.at(1) == zero);

    CHECK_THROWS( tr.at(2) );
}


TEST_CASE( "trie update/erase/insert", "[trie]" ) {
    vl_trie tr = vl_trie(4,4);
    lineral f = lineral({2,3,4});
    var_t vf = 0;
    lineral g = lineral({1,3,4});
    var_t vg = 1;
    lineral h = lineral({1,3});
    var_t vh = 2;
    
    auto ins1 = tr.insert(vf, f, 0);
    CHECK( tr.size()==1 );
    CHECK( tr.get_num_nodes()==4 );
    CHECK( tr.to_str() == "(0,x2+x3+x4)" );

    auto ins2 = tr.insert(vg, g, 0);
    CHECK( tr.size()==2 );
    CHECK( tr.get_num_nodes()==5 );
    CHECK( tr.to_str() == "(0,x2+x3+x4) (1,x1+x3+x4)" );

    auto ins3 = tr.insert(vh, h, 0);
    CHECK( tr.size()==3 );
    CHECK( tr.get_num_nodes()==7 );
    CHECK( tr.to_str() == "(0,x2+x3+x4) (1,x1+x3+x4) (2,x1+x3)" );

    lineral r = lineral({4});
    LinEqs L = LinEqs({r});
    CHECK( L.to_str() == "x4");

    tr.erase(vf);
    CHECK( tr.to_str() == "(1,x1+x3+x4) (2,x1+x3)" );

    tr.insert(vf,f, 0);
    CHECK( tr.to_str() == "(0,x2+x3+x4) (1,x1+x3+x4) (2,x1+x3)" );

    lineral l = lineral({1,2});
    auto [v,b] = tr.update(vf, l, 0);
    CHECK( !b );
    CHECK( v==vf );
    CHECK( tr.to_str() == "(0,x1+x2) (1,x1+x3+x4) (2,x1+x3)" );
    
    auto [v_,b_] = tr.update(vf, g, 0);
    CHECK( v_==vg );
    CHECK( !b_ );
    CHECK( tr.to_str() == "(1,x1+x3+x4) (2,x1+x3)" );

    auto [v__,b__] = tr.update(vg, h.plus_one(), 0);
    CHECK( v__==vh );
    CHECK( b__ );
    CHECK( tr.to_str() == "(2,x1+x3)" );
}

/*
TEST_CASE( "trie impl graph test", "[trie]" ) {
    //auto clss = parse_file("../../benchmarks/instances/2xnfs/rand/rand-10-20_2.xnf");
    //auto clss = parse_file("../../benchmarks/instances/2xnfs/bivium/bivium-t100-g40.xnf");
    //auto clss = parse_file("../../benchmarks/instances/2xnfs/mq/toyexamples/ToyExample-type1-n10-seed0.xnf");
    auto clss = parse_file("../../benchmarks/instances/2xnfs/mq/challenge-1-55-0.xnf");
    auto IG = impl_graph(clss);

    vl_trie tr = vl_trie(IG.get_no_v(), IG.get_opts()->num_vars);

    n_t sz = 0;
    n_t sum_idxs_sz = 0;
    for (const auto &[v,lit] : IG.get_Vxlit()) {
        const auto ins = tr.insert( v, lit , 0);
        sz++;
        sum_idxs_sz += lit.get_idxs().size();
        REQUIRE( ins.inserted );
        REQUIRE( tr.size() == sz );
    }
    
    for (const auto &[v,lit] : IG.get_Vxlit()) {
        REQUIRE( tr[v]==lit );
        REQUIRE( tr[lit]==v );
    }

    std::cout << "num vars     : " << IG.get_opts()->num_vars << std::endl;
    std::cout << "IG verts     : " << IG.get_no_v() << std::endl;
    std::cout << "lits stored  : " << tr.size() << std::endl;
    std::cout << "sum lit size : " << sum_idxs_sz << std::endl;
    std::cout << "trie nodes   : " << tr.nodes() << std::endl;
}
*/

TEST_CASE( "trie copy ctor", "[trie]" ) {
    vl_trie tr(10,5);

    tr.insert(1, lineral({1,2}), 0);

    CHECK( tr.to_str() == "(1,x1+x2)" );

    vl_trie tr_cpy(tr);
    
    CHECK( tr_cpy.to_str() == "(1,x1+x2)" );
    
    //SECTION("del tr, del tr_cpy") {
    //    tr.~vl_trie();

    //    REQUIRE( tr_cpy.to_str() == "(1,x1+x2)" );
    //    tr_cpy.~vl_trie();
    //}

    //SECTION("del tr") {
    //    tr.~vl_trie();

    //    CHECK( tr_cpy.to_str() == "(1,x1+x2)" );

    //    CHECK( tr_cpy[1].to_str()=="x1+x2+x3" );
    //    CHECK( tr_cpy[lineral({1,2},5)]==1 );
    //}

    SECTION("edit tr") {
        auto ins1 = tr.insert(4, lineral({2,4}), 0);
        CHECK( ins1.inserted );
        CHECK( tr_cpy.to_str() == "(1,x1+x2)" );

        auto ins1_cpy = tr_cpy.insert(4, lineral({2,4}), 0);
        CHECK( ins1_cpy.inserted );
    }
}

TEST_CASE( "trie iterator", "[trie]") {
    vl_trie tr = vl_trie(4,4);
    lineral f = lineral({2,3,4});
    var_t vf = 0;
    lineral g = lineral({1,3,4});
    var_t vg = 1;
    lineral h = lineral({1,3});
    var_t vh = 2;
    
    auto ins1 = tr.insert(vf, f, 0);
    auto ins2 = tr.insert(vg, g, 0);
    auto ins3 = tr.insert(vh, h, 0);

    //manually check one iterator
    auto it = tr.begin(0);
    CHECK(*it == 2);
    ++it;
    CHECK(*it == 3);
    ++it;
    CHECK(*it == 4);
    ++it;
    CHECK(it == tr.end());

    //check sum-func
    CHECK( tr.sum(0,1) == lineral({1,2}) );
    CHECK( tr.sum(0,2) == lineral({1,2,4}) );
    CHECK( tr.sum(1,2) == lineral({4}) );

}