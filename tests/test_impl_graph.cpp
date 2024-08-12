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

//file to test implementation of xcls
#include "../src/impl_graph.hpp"
#include "../src/solve.hpp"

#include <catch2/catch_all.hpp>


TEST_CASE( "implication graph construction", "[graph][impl-graph]" ) {
    //construct list of xor-clauses
    vec< vec<lineral> > clss;
    options opt(4);

    SECTION( "simple graph" ){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({1}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,2,3})), lineral(vec<var_t>({2,3}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1,2,4}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,4}))}) );
        //reduces only to linear and empty clauses!

        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "x1 x2+1 x4+1" );
    }
    
    SECTION( "single clause"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1)\n0" );
    }

    SECTION( "non-trivial graph"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({2})), lineral(vec<var_t>({3}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1) (x2+1,x2+x3+1) (x2+1,x3); (x2+x3,x2) (x2+x3,x3); (x3+1,x2) (x3+1,x2+x3+1) (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0" );
        
        
        //construct IG -- simple
        opt.ext = constr::simple;
        auto IG2 = impl_graph(clss, opt);
        CHECK( IG2.to_str() == "(x1+1,x2); (x2+1,x1) (x2+1,x3); (x3+1,x2) (x3+1,x4); (x4+1,x3)\n0" );
    }

}

TEST_CASE( "implication graph analysis - scc and to", "[graph][impl-graph][scc][to]" ) {
    vec< vec<lineral> > clss;
    options opt(4);
    opt.ext = constr::simple;

    SECTION("two symmetrical comps") {
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,2})), lineral(vec<var_t>({3}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,3})), lineral(vec<var_t>({1}))}) );
    
        //construct IG -- simple
        impl_graph IG(clss, opt);
        CHECK( IG.to_str() == "(x1+1,x3+1); (x1,x2); (x2+1,x1+1); (x2,x3); (x3+1,x2+1); (x3,x1)\n0" );

        LinEqs L = IG.scc_analysis();
        CHECK( L.to_str() == "x1+x3 x2+x3");
        CHECK( L.is_consistent() );
    }
    
    SECTION("a self-symmetrical comp") {
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,2})), lineral(vec<var_t>({0,1}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({0,3}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,4})), lineral(vec<var_t>({1}))}) );
    
        //construct IG -- simple
        impl_graph IG(clss, opt);
        CHECK( IG.to_str() == "(x1+1,x3+1) (x1+1,x4+1); (x1,x2) (x1,x2+1); (x2+1,x1+1); (x2,x1+1); (x3+1,x4); (x3,x1); (x4+1,x3); (x4,x1)\n0" );

        LinEqs L = IG.scc_analysis();
        //CHECK( L.to_str() == "x1+x4 x2+x4 x3+x4 1");
        CHECK(L.dim() == 4);
        CHECK( !L.is_consistent() );
    }
    
    SECTION("topological ordering of disconnected paths") {
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,2})), lineral(vec<var_t>({3}))}) );
    
        //construct IG -- simple
        impl_graph IG(clss, opt);
        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1); (x2,x3); (x3+1,x2+1)\n0" );

        auto to = IG.get_TO();
        CHECK( to == vec<var_t>({0,5,2,3,4,1}) ); //corr to: [x1, x3+1, x2, x2+1, x3, x1+1]
    }
}


TEST_CASE( "update implication graph", "[graph][impl-graph][update]" ) {
    //construct list of xor-clauses
    vec< vec<lineral> > clss;
    options opt(4);
    
    SECTION( "single clause update -- simple -- no new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({0,1}))};
        LinEqs L = LinEqs( linerals );
        
        CHECK( L.to_str() == "x1+1" );
        IG.add_new_xsys(L);

        LinEqs L_new = IG.update_graph( L );

        CHECK( L_new.to_str() == "0" );
    }
    
    SECTION( "single clause update -- simple -- new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral({1})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1" );

        LinEqs L_new = IG.update_graph( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "single clause - new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1)\n0" );
        
        vec<lineral> linerals = { lineral({1,2})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x2" );

        LinEqs L_new = IG.update_graph( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "merging two clauses1"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3)\n0" );
        
        vec<lineral> linerals = { lineral({1,3}), lineral({2,4})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x3 x2+x4" );

        LinEqs L_new = IG.update_graph( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x4); (x4+1,x3)\n0\nx1+x3 x2+x4" );
    }
    
    SECTION( "merging two clauses2"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1); (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0" );
        
        vec<lineral> xlits1 = { lineral(vec<var_t>({1,3}))};
        LinEqs L1 = LinEqs( xlits1 );
        IG.add_new_xsys(L1);
        
        CHECK( L1.to_str() == "x1+x3" );
        
        LinEqs L_new1 = IG.update_graph( L1 );
        IG.add_new_xsys(L_new1);

        CHECK( L_new1.to_str() == "0" );

        CHECK( IG.to_str() == "(x2+1,x2+x3+1) (x2+1,x3); (x2+x3,x2) (x2+x3,x3); (x3+1,x2) (x3+1,x2+x3+1) (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3\n0" );
        
        vec<lineral> xlits2 = {lineral(vec<var_t>({2,4}))};
        LinEqs L2 = LinEqs( xlits2 );
        IG.add_new_xsys(L2);
        
        CHECK( L2.to_str() == "x2+x4" );
        
        LinEqs L_new2 = IG.update_graph( L2 );

        CHECK( L_new1.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3\n0\nx2+x4" );
    }
    
    SECTION( "merging three clauses"){
        opt.num_vars = 6;

        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({5})), lineral(vec<var_t>({6}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3); (x5+1,x6); (x6+1,x5)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({1,3})), lineral(vec<var_t>({1,5}))};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys( L );
        
        CHECK( L.to_str() == "x1+x5 x3+x5" );

        LinEqs L_new = IG.update_graph( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x2+1,x5); (x4+1,x5); (x5+1,x2) (x5+1,x4) (x5+1,x6); (x6+1,x5)\n0\nx1+x5 x3+x5" );
    }
}

TEST_CASE( "update-hfd implication graph", "[graph][impl-graph][update-hfd]" ) {
    //construct list of xor-clauses
    vec< vec<lineral> > clss;
    options opt(4);
    
    SECTION( "single clause update -- simple -- no new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({0,1}))};
        LinEqs L = LinEqs( linerals );
        
        CHECK( L.to_str() == "x1+1" );

        LinEqs L_new = IG.update_graph_hash_fight_dev( L );

        CHECK( L_new.to_str() == "0" );
    }
    
    SECTION( "single clause update -- simple -- new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral({1})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1" );

        LinEqs L_new = IG.update_graph_hash_fight_dev( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "single clause - new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1)\n0" );
        
        vec<lineral> linerals = { lineral({1,2})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x2" );

        LinEqs L_new = IG.update_graph_hash_fight_dev( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "merging two clauses1"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3)\n0" );
        
        vec<lineral> linerals = { lineral({1,3}), lineral({2,4})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x3 x2+x4" );

        LinEqs L_new = IG.update_graph_hash_fight_dev( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x4); (x4+1,x3)\n0\nx1+x3 x2+x4" );
    }
    
    SECTION( "merging two clauses2"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1); (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0" );
        
        vec<lineral> xlits1 = { lineral(vec<var_t>({1,3}))};
        LinEqs L1 = LinEqs( xlits1 );
        IG.add_new_xsys(L1);
        
        CHECK( L1.to_str() == "x1+x3" );
        
        LinEqs L_new1 = IG.update_graph_hash_fight_dev( L1 );

        CHECK( L_new1.to_str() == "0" );

        CHECK( IG.to_str() == "(x2+1,x2+x3+1) (x2+1,x3); (x2+x3,x2) (x2+x3,x3); (x3+1,x2) (x3+1,x2+x3+1) (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3" );
        //CHECK( IG.to_str() == "(x2+1,x2+x3+1) (x2+1,x3); (x2+x3,x2) (x2+x3,x3); (x3+1,x2) (x3+1,x2+x3+1) (x3+1,x3+x4+1) (x3+1,x4); (x4+1,x3) (x4+1,x3+x4+1); (x3+x4,x3) (x3+x4,x4)\n0" );
        
        vec<lineral> xlits2 = {lineral(vec<var_t>({2,4}))};
        LinEqs L2 = LinEqs( xlits2 );
        IG.add_new_xsys(L2);
        
        CHECK( L2.to_str() == "x2+x4" );
        
        LinEqs L_new2 = IG.update_graph_hash_fight_dev( L2 );

        CHECK( L_new1.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3\nx2+x4" );
    }
    
    SECTION( "merging three clauses"){
        opt.num_vars = 6;

        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({5})), lineral(vec<var_t>({6}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3); (x5+1,x6); (x6+1,x5)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({1,3})), lineral(vec<var_t>({1,5}))};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x5 x3+x5" );

        LinEqs L_new = IG.update_graph_hash_fight_dev( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x2+1,x5); (x4+1,x5); (x5+1,x2) (x5+1,x4) (x5+1,x6); (x6+1,x5)\n0\nx1+x5 x3+x5" );
        //CHECK( IG.to_str() == "(x2+1,x5); (x4+1,x5); (x5+1,x2) (x5+1,x4) (x5+1,x6); (x6+1,x5)\n0" );
    }
}

TEST_CASE( "update-hf implication graph", "[graph][impl-graph][update-hf]" ) {
    //construct list of xor-clauses
    vec< vec<lineral> > clss;
    options opt(4);
    
    SECTION( "single clause update -- simple -- no new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({0,1}))};
        LinEqs L = LinEqs( linerals );
        
        CHECK( L.to_str() == "x1+1" );

        LinEqs L_new = IG.update_graph_hash_fight( L );

        CHECK( L_new.to_str() == "0" );
    }
    
    SECTION( "single clause update -- simple -- new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral({1})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1" );

        LinEqs L_new = IG.update_graph_hash_fight( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "single clause - new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1)\n0" );
        
        vec<lineral> linerals = { lineral({1,2})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x2" );

        LinEqs L_new = IG.update_graph_hash_fight( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "merging two clauses1"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3)\n0" );
        
        vec<lineral> linerals = { lineral({1,3}), lineral({2,4})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x3 x2+x4" );

        LinEqs L_new = IG.update_graph_hash_fight( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x4); (x4+1,x3)\n0\nx1+x3 x2+x4" );
    }
    
    SECTION( "merging two clauses2"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1); (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0" );
        
        vec<lineral> xlits1 = { lineral(vec<var_t>({1,3}))};
        LinEqs L1 = LinEqs( xlits1 );
        IG.add_new_xsys(L1);
        
        CHECK( L1.to_str() == "x1+x3" );
        
        LinEqs L_new1 = IG.update_graph_hash_fight( L1 );

        CHECK( L_new1.to_str() == "0" );

        CHECK( IG.to_str() == "(x2+1,x2+x3+1) (x2+1,x3); (x2+x3,x2) (x2+x3,x3); (x3+1,x2) (x3+1,x2+x3+1) (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3" );
        //CHECK( IG.to_str() == "(x2+1,x2+x3+1) (x2+1,x3); (x2+x3,x2) (x2+x3,x3); (x3+1,x2) (x3+1,x2+x3+1) (x3+1,x3+x4+1) (x3+1,x4); (x4+1,x3) (x4+1,x3+x4+1); (x3+x4,x3) (x3+x4,x4)\n0" );
        
        vec<lineral> xlits2 = {lineral(vec<var_t>({2,4}))};
        LinEqs L2 = LinEqs( xlits2 );
        IG.add_new_xsys(L2);
        
        CHECK( L2.to_str() == "x2+x4" );
        
        LinEqs L_new2 = IG.update_graph_hash_fight( L2 );

        CHECK( L_new1.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3\nx2+x4" );
    }
    
    SECTION( "merging three clauses"){
        opt.num_vars = 6;

        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({5})), lineral(vec<var_t>({6}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3); (x5+1,x6); (x6+1,x5)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({1,3})), lineral(vec<var_t>({1,5}))};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x5 x3+x5" );

        LinEqs L_new = IG.update_graph_hash_fight( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x2+1,x5); (x4+1,x5); (x5+1,x2) (x5+1,x4) (x5+1,x6); (x6+1,x5)\n0\nx1+x5 x3+x5" );
        //CHECK( IG.to_str() == "(x2+1,x5); (x4+1,x5); (x5+1,x2) (x5+1,x4) (x5+1,x6); (x6+1,x5)\n0" );
    }
}


TEST_CASE( "update-par implication graph", "[graph][impl-graph][update-par]" ) {
    //construct list of xor-clauses
    vec< vec<lineral> > clss;
    options opt(4);

    omp_set_num_threads( omp_get_max_threads() ); //activate parallelism!
    
    SECTION( "single clause update -- simple -- no new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({0,1}))};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+1" );

        LinEqs L_new = IG.update_graph_par( L );

        CHECK( L_new.to_str() == "0" );
    }
    
    SECTION( "single clause update -- simple -- new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1,x2); (x2+1,x1+1)\n0" );
        
        vec<lineral> linerals = { lineral({1})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1" );

        LinEqs L_new = IG.update_graph_par( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "single clause - new linerals"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1)\n0" );
        
        vec<lineral> linerals = { lineral({1,2})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x2" );

        LinEqs L_new = IG.update_graph_par( L );

        CHECK( L_new.to_str() == "x2" );
    }
    
    SECTION( "merging two clauses"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3)\n0" );
        
        vec<lineral> linerals = { lineral({1,3}), lineral({2,4})};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x3 x2+x4" );

        LinEqs L_new = IG.update_graph_par( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x4); (x4+1,x3)\n0\nx1+x3 x2+x4" );
    }
    
    SECTION( "merging two clauses"){
        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        
        //construct IG
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x1+x2+1) (x1+1,x2); (x1+x2,x1) (x1+x2,x2); (x2+1,x1) (x2+1,x1+x2+1); (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0" );
        
        vec<lineral> xlits1 = { lineral(vec<var_t>({1,3}))};
        LinEqs L1 = LinEqs( xlits1 );
        IG.add_new_xsys(L1);
        
        CHECK( L1.to_str() == "x1+x3" );
        
        LinEqs L_new1 = IG.update_graph_par( L1 );

        CHECK( L_new1.to_str() == "0" );

        CHECK( IG.to_str() == "(x2+1,x2+x3+1) (x2+1,x3); (x2+x3,x2) (x2+x3,x3); (x3+1,x2) (x3+1,x2+x3+1) (x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3" );
        
        vec<lineral> xlits2 = {lineral(vec<var_t>({2,4}))};
        LinEqs L2 = LinEqs( xlits2 );
        IG.add_new_xsys(L2);
        
        CHECK( L2.to_str() == "x2+x4" );
        
        LinEqs L_new2 = IG.update_graph_par( L2 );

        CHECK( L_new1.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x3+1,x3+x4+1) (x3+1,x4); (x3+x4,x3) (x3+x4,x4); (x4+1,x3) (x4+1,x3+x4+1)\n0\nx1+x3\nx2+x4" );
    }
    
    SECTION( "merging three clauses"){
        opt.num_vars = 6;

        clss.push_back( vec<lineral>({lineral(vec<var_t>({1})), lineral(vec<var_t>({2}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({3})), lineral(vec<var_t>({4}))}) );
        clss.push_back( vec<lineral>({lineral(vec<var_t>({5})), lineral(vec<var_t>({6}))}) );
        
        //construct IG
        opt.ext = constr::simple;
        impl_graph IG(clss, opt);

        CHECK( IG.to_str() == "(x1+1,x2); (x2+1,x1); (x3+1,x4); (x4+1,x3); (x5+1,x6); (x6+1,x5)\n0" );
        
        vec<lineral> linerals = { lineral(vec<var_t>({1,3})), lineral(vec<var_t>({1,5}))};
        LinEqs L = LinEqs( linerals );
        IG.add_new_xsys(L);
        
        CHECK( L.to_str() == "x1+x5 x3+x5" );

        LinEqs L_new = IG.update_graph_par( L );

        CHECK( L_new.to_str() == "0" );
        
        CHECK( IG.to_str() == "(x2+1,x5); (x4+1,x5); (x5+1,x2) (x5+1,x4) (x5+1,x6); (x6+1,x5)\n0\nx1+x5 x3+x5" );
    }
}

TEST_CASE( "parsing and solving instances/2xnfs test instances" , "[impl-graph][graph][parser][scc]" ) {
    SECTION( "test0.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test0.xnf");
        auto IG = impl_graph(clss);
        CHECK( IG.assert_data_structs() );
        //CHECK( IG.to_str() == "(x1,x2+1) (x1,x3); (x2+1,x3+1); (x2,x1+1); (x3+1,x1+1); (x3,x2)\n0" );
        auto L = IG.scc_analysis();
        CHECK( L.to_str() == "0" );
    
        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT
        //CHECK( s.sol == vec<bool>({true,false,true}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test1.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test1.xnf");
        auto IG = impl_graph(clss);
        CHECK( IG.assert_data_structs() );
        auto L = IG.scc_analysis();
        CHECK( L.to_str() == "x1+x5 x2+x5 x3+x5 x4+x5" );
    
        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT
        //CHECK( s.sol == vec<bool>({false,false,false,false,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test2.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test2.xnf");
        auto IG = impl_graph(clss);
    
        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //UNSAT
        //CHECK( s.sol == vec<bool>({false,false,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test3.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test3.xnf");
        auto IG = impl_graph(clss);
        //CHECK( IG.to_str() == "x1+x5 x2+x5 x3+x5 x4+x5 1" );
    
        stats s = IG.dpll_solve();
        CHECK( s.sat == false ); //UNSAT
    }
    
    SECTION( "test4.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test4.xnf");
        auto IG = impl_graph(clss);
        
        //SCCs lead to inconsistent eqs
        auto L = IG.scc_analysis();
        //CHECK( L.to_str() == "x1+x5+x7 x2+x5+x7 x3+x5+x7 x4+x5+x7 x6+x7 1" );

        stats s = IG.dpll_solve();
        CHECK( s.sat == false ); //UNSAT!
    }

    SECTION( "test5.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test5.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({false,false,false,false,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test6.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test6.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({false,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test6_.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test6_.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true );
        //CHECK( s.sol == vec<bool>({false,true}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test7.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test7.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true );
        //CHECK( s.sol == vec<bool>({true,false,false,true,true}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test8.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test8.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({true,true,true}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test9.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test9.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({true,false,true,false,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test10.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test10.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({true,true,false,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test11.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test11.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test12.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test12.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test13.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test13.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test14.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test14.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test15.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test15.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test16.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test16.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test17.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test17.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test18.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test18.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "test19.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test19.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test20.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test20.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test21.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test21.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test22.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test22.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test23.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test23.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test24.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test24.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "test25.xnf" ) {
        auto clss = parse_file("tests/2xnfs/test25.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }
}

TEST_CASE( "solving 2xnf test instances with gp" , "[solver][gp]" ) {
    SECTION( "test33.xnf" ) {
        reordering P;
        P.insert(2,1);
        P.insert(1,2);
        P.insert(4,3);
        P.insert(5,4);
        P.insert(3,5);
        auto clss_gp = parse_file_gp("tests/2xnfs/test33.xnf", P);
        auto slvr_gp = impl_graph(clss_gp);
        auto clss = parse_file("tests/2xnfs/test33.xnf");
        auto slvr = impl_graph(clss);

        stats s_gp = slvr_gp.dpll_solve();
        stats s = slvr.dpll_solve();

        CHECK( s.sat == true ); //SAT!
        CHECK( s_gp.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
        CHECK( check_sol(clss_gp.cls, s_gp.sol) );

        //check that after reordering sol of s_gp, we get a sol of original clss
        CHECK( !check_sol(clss.cls, s_gp.sol) );
        s_gp.reorder_sol(P);
        CHECK( check_sol(clss.cls, s_gp.sol) );
    }
    
    SECTION( "flat30-100.xnf" ) {
        reordering P;
        P.insert(2,1);
        P.insert(1,2);
        //auto P = parse_gp("gp");
        auto clss_gp = parse_file_gp("tests/2xnfs/flat30-100.xnf", P);
        auto slvr_gp = impl_graph(clss_gp);
        slvr_gp.get_opts()->dh = dec_heu::lex;
        auto clss = parse_file("tests/2xnfs/flat30-100.xnf");
        auto slvr = impl_graph(clss);
        slvr.get_opts()->dh = dec_heu::lex;

        stats s_gp = slvr_gp.dpll_solve();
        stats s = slvr.dpll_solve();

        CHECK( s.sat == true ); //SAT!
        CHECK( s_gp.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
        CHECK( check_sol(clss_gp.cls, s_gp.sol) );

        //check that after reordering sol of s_gp, we get a sol of original clss
        CHECK( !check_sol(clss.cls, s_gp.sol) );
        s_gp.reorder_sol(P);
        CHECK( check_sol(clss.cls, s_gp.sol) );
    }
}


TEST_CASE( "implication graph update + backtracking", "[graph][impl-graph][update][backtrack]" ) {
    //solve some first more 'serious' examples

    SECTION( "rand-3-6" ) {
        auto clss = parse_file("tests/2xnfs/rand-3-6.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({true,false,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "rand-5-10" ) {
        auto clss = parse_file("tests/2xnfs/rand-5-10.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({true,false,false,false,true}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "rand-10-20" ) {
        auto clss = parse_file("tests/2xnfs/rand-10-20.xnf");
        auto IG = impl_graph(clss);
        IG.get_opts()->dh = dec_heu::mr;
        IG.get_opts()->verb = 100;

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({false,false,false,false,true,false,false,false,true,true}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "rand-10-20_2" ) {
        auto clss = parse_file("tests/2xnfs/rand-10-20_2.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({true,false,true,true,true,false,false,false,true,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "rand-10-20_2" ) {
        auto clss = parse_file("tests/2xnfs/rand-10-20_2.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({true,false,true,true,true,false,false,false,true,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
 
    SECTION( "rand-10-30" ) {
        auto clss = parse_file("tests/2xnfs/rand-10-30.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        //CHECK( s.sol == vec<bool>({false,true,false,true,true,true,true,false,true,false}) );
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "rand-20-60" ) {
        auto clss = parse_file("tests/2xnfs/rand-20-60.xnf");
        auto IG = impl_graph(clss);

        stats s = IG.dpll_solve();
        CHECK( s.sat == true ); //SAT!
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    //SECTION( "bivium-t500-g177" ) {
    //    auto clss = parse_file("tests/2xnfs/bivium/bivium-t500-g177.xnf");
    //    auto IG = impl_graph(clss);
    //    stats s = IG.dpll_solve();
    //    CHECK( s.sat == true ); //SAT!
    //    CHECK( check_sol(clss.cls, s.sol) );
    //}
}
