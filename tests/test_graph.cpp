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

//file to test implementation of LHGR
#include <vector>
#include <set>

#include "../src/graph/graph.hpp"

#include <catch2/catch_all.hpp>


void CHECK_EQ_SET(vec<var_t> a, vec<var_t> b) {
    CHECK(std::set<var_t>(a.begin(),a.end()) == std::set<var_t>(b.begin(),b.end()));
}


TEST_CASE( "graph creation, edge-removal, and backtracking (small undirected graph, non-trivial sigma)", "[LHDGR]" ) {
    //cycle of length 3
    vec<var_t> sigma = {1,0,3,2};
    vec< std::pair<var_t,var_t>> E = {std::pair<var_t,var_t>(0,2),
                                              std::pair<var_t,var_t>(1,3)};

    graph G = graph(E,sigma.size());

    CHECK( G.assert_data_structs() );

    CHECK( G.to_str() == "(0,2); (1,3); (2,0); (3,1)");
    
    graph_repr G_orig = G.get_state();

    SECTION("remove_edge") {
        #if USE_LHGR
            G.remove_edge(0, 0); //remove first edge outgoing from 0, i.e., (0,1) OR (0,2)
        #else
            G.remove_edge(0, 2);
        #endif

        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,3); (2,0)");
    }
    
    SECTION("remove_all_edges") {
        G.remove_all_edges(0); //remove all edges outgoing from, i.e., (0,2)

        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,3); (2,0)");

    }
    
    G.backtrack( std::move(G_orig) );
        
    CHECK( G.to_str() == "(0,2); (1,3); (2,0); (3,1)");
}

/*
TEST_CASE( "graph creation, edge-removal, and backtracking (small undirected graph, trivial sigma)", "[LHDGR]" ) {
    //cycle of length 3
    vec<var_t> sigma = {0,1,2};
    vec< std::pair<var_t,var_t>> E = {std::pair<var_t,var_t>(0,1),
                                              std::pair<var_t,var_t>(1,2),
                                              std::pair<var_t,var_t>(2,0)};
                                            
    graph G = graph(E,3);

    CHECK( G.assert_data_structs() );
    CHECK( G.to_str() == "(0,1) (0,2); (1,0) (1,2); (2,0) (2,1)");

    G.assert_data_structs();

    graph_repr G_orig = G.get_state();

    #ifdef USE_LHGR
        G.remove_edge(0, 0); //remove first edge outgoing from 0, i.e., (0,1) OR (0,2)
    #else
        G.remove_edge(0, 1);
    #endif

    CHECK( G.assert_data_structs() );
    CHECK( G.to_str() == "(0,2); (1,2); (2,0) (2,1)");

    G.backtrack( std::move(G_orig) );
    
    CHECK( G.to_str() == "(0,1) (0,2); (1,0) (1,2); (2,0) (2,1)");
}
*/

TEST_CASE( "graph creation, edge-removal, vertex merging, and backtracking (small directed graph, non-trivial sigma)", "[LHDGR]" ) {
    //cycle of length 3
    vec<var_t> sigma = {1,0,3,2,5,4};
    vec< std::pair<var_t,var_t>> E = {std::pair<var_t,var_t>(5,3),
                                              std::pair<var_t,var_t>(4,0),
                                              std::pair<var_t,var_t>(1,2)};

    graph G = graph(E,6);

    CHECK( G.assert_data_structs() );

    CHECK( G.to_str() == "(1,2) (1,5); (2,4); (3,0); (4,0); (5,3)");
        
    graph_repr G_orig = G.get_state();
    
    CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
    CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({3,4}) );
    CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({5,2}) );
    CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
    CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({4}) );
    CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({1}) );
    CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({0}) );
    CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({5}) );
    CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({0}) );
    CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({2}) );
    CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({3}) );
    CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({1}) );


    SECTION("remove_edge") {
        #ifdef USE_LHGR
            G.remove_edge(2, 0); //remove first edge outgoing from 2, i.e., (2,4) (and symmetric edge (5,3))
        #else
            G.remove_edge(2, 4); //remove first edge outgoing from 2, i.e., (2,4) (and symmetric edge (5,3))
        #endif
        CHECK(G.assert_data_structs());

        CHECK(G.to_str() == "(1,2) (1,5); (3,0); (4,0)");

        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({4,3}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({5,2}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({1}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({1}) );
    };
    
    SECTION("remove_edge") {
        #ifdef USE_LHGR
            G.remove_edge(1, 0); //remove first edge outgoing from 1, i.e., (1,5) (and symmetric edge (4,0))
        #else
            G.remove_edge(1, 5); //remove first edge outgoing from 1, i.e., (1,5) (and symmetric edge (4,0))
        #endif
        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,2); (2,4); (3,0); (5,3)");
        
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({3}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({2}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({1}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({5}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({2}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({3}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({}) );
    };
    
    SECTION("remove_edge") {
        #ifdef USE_LHGR
            G.remove_edge(1, 1); //remove second edge outgoing from 1, i.e., (1,2) (and symmetric edge (3,0))
        #else
            G.remove_edge(1, 2); //remove second edge outgoing from 1, i.e., (1,2) (and symmetric edge (3,0))
        #endif
        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,5); (2,4); (4,0); (5,3)");

        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({5}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({5}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({2}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({3}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({1}) );
    };
    
    SECTION("remove_vert") {
        G.remove_vert(1); //remove vertex 1, and its symmetric vertex 0

        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(2,4); (5,3)");
        
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({5}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({2}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({3}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({}) );
    }
    
    SECTION("remove_vert") {
        G.remove_vert(0); //remove vertex 0, and its symmetric vertex 1

        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(2,4); (5,3)");
        
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({5}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({2}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({3}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({}) );
    }
    
    SECTION("remove_vert") {
        G.remove_vert(5); //remove vertex 5, and its symmetric vertex 4
        
        G.assert_data_structs();

        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,2); (3,0)");

        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({3}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({2}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({1}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({}) );
    }
    
    SECTION("vertex_merging") {
        G.merge_verts(1,2); //merge verts 1,2 (and 0,3)

        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,4) (1,5); (4,0); (5,0)");
        
        #ifdef USE_LHGR
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), G.get_out_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), G.get_in_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), G.get_out_neighbour_vector(2) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), G.get_in_neighbour_vector(2) );
        #endif
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({5,4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({4,5}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({1}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({1}) );
    }
    #ifdef USE_LHGR //merging symmetric verts not supported!
    SECTION("vertex_merging (symmetric pair!)") {
        G.merge_verts(4,5); //merge verts 4 and 5

        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,2) (1,4); (2,4); (3,0); (4,0) (4,3)" );
        
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({4,3}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({4,2}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(2), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(2), vec<var_t>({1}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(3), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(3), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), G.get_out_neighbour_vector(5) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), G.get_in_neighbour_vector(5) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({0,3}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({1,2}) );
    }
    #endif
    
    SECTION("double vertex merging (no-overlap)") {
        G.merge_verts(1,2); //merge verts 1 and 2

        CHECK(G.assert_data_structs());

        G.merge_verts(4,5);
        
        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,4); (4,0)");

        #ifdef USE_LHGR
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), G.get_out_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), G.get_in_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), G.get_out_neighbour_vector(2) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), G.get_in_neighbour_vector(2) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({4}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), G.get_out_neighbour_vector(5) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), G.get_in_neighbour_vector(5) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({1}) );
        #endif
    }
    

    SECTION("double vertex merging (enlarging non-trivial color) and backtracking") {
        G.merge_verts(1,2); //merge verts 1 and 2

        CHECK(G.assert_data_structs());
        
        CHECK( G.to_str() == "(1,4) (1,5); (4,0); (5,0)");

        graph_repr G_orig_ = G.get_state();

        G.merge_verts(1,4);
        
        CHECK(G.assert_data_structs());

        CHECK( G.to_str() == "(1,0)"); //reduces to a single self-symmetric edge!

        #ifdef USE_LHGR
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), G.get_out_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), G.get_in_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), G.get_out_neighbour_vector(2) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), G.get_in_neighbour_vector(2) );
        #endif
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({1}) );
        #ifdef USE_LHGR
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), G.get_out_neighbour_vector(2) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), G.get_in_neighbour_vector(2) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), G.get_out_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), G.get_in_neighbour_vector(3) );
        #endif

        //backtrack 1-level
        G.backtrack( std::move(G_orig_) );
        
        CHECK( G.to_str() == "(1,4) (1,5); (4,0); (5,0)");

        #ifdef USE_LHGR
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), G.get_out_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), G.get_in_neighbour_vector(3) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), G.get_out_neighbour_vector(2) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), G.get_in_neighbour_vector(2) );
        #endif
        CHECK_EQ_SET( G.get_out_neighbour_vector(1), vec<var_t>({5,4}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(1), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(0), vec<var_t>({}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(0), vec<var_t>({4,5}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(4), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(4), vec<var_t>({1}) );
        CHECK_EQ_SET( G.get_out_neighbour_vector(5), vec<var_t>({0}) );
        CHECK_EQ_SET( G.get_in_neighbour_vector(5), vec<var_t>({1}) );
    }
    
    G.backtrack( std::move(G_orig) );

    G.to_str();
        
    CHECK( G.to_str() == "(1,2) (1,5); (2,4); (3,0); (4,0); (5,3)");
}

TEST_CASE( "edge case in vertex-merging", "[LHDGR]" ) {
    vec<var_t> sigma = {1,0,3,2};
    vec< std::pair<var_t,var_t>> E = {std::pair<var_t,var_t>(0,2)};

    graph G = graph(E,4);
    
    CHECK( G.to_str() == "(0,2); (3,1)");

    graph_repr G_orig = G.get_state();

    SECTION("merge to cycle") {
        G.merge_verts(1,2);

        CHECK( G.to_str() == "(0,1)"); //symmetric edges now coincide!
    }
    
    SECTION("merge to cycle (with different label for cols)") {
        G.merge_verts(2,1);

        CHECK( G.to_str() == "(3,2)"); //symmetric edges now coincide!
    }
    
    SECTION("merge to disconnect") {
        G.merge_verts(1,3);

        CHECK( G.to_str() == ""); //edges become in-color, i.e., are removed
    }
    
    G.backtrack( std::move(G_orig) );
    
    CHECK( G.to_str() == "(0,2); (3,1)");
}

TEST_CASE( "graph creation, edge-removal, vertex merging, and backtracking (longer non-trivial example)", "[LHDGR]" ) {
    //cycle of length 3
    vec<var_t> sigma = {1,0,3,2,5,4,7,6,9,8,11,10};
    vec< std::pair<var_t,var_t>> E = {std::pair<var_t,var_t>(0,10),
                                              std::pair<var_t,var_t>(3,9),
                                              std::pair<var_t,var_t>(2,4),
                                              std::pair<var_t,var_t>(6,7),
                                              std::pair<var_t,var_t>(6,1),
                                              std::pair<var_t,var_t>(8,5)};

    auto state_str_stack = std::list< std::pair<graph_repr, std::string> >();

    graph G = graph(E,12);

    state_str_stack.push_back( 
        std::pair<graph_repr, std::string>(G.get_state(), "(0,7) (0,10); (2,4); (3,9); (4,9); (5,3); (6,1) (6,7); (8,2) (8,5); (11,1)")
        );
    CHECK( G.to_str() == state_str_stack.back().second );
    CHECK( G.assert_data_structs() );
    CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5,6,7,8,9,10,11}) );

    #ifdef USE_LHGR
    SECTION("edge cases...") {
        graph_repr G_orig = G.get_state();

        G.merge_verts(6,7);

        CHECK( G.to_str() == "(0,6) (0,10); (2,4); (3,9); (4,9); (5,3); (6,1); (8,2) (8,5); (11,1)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5,6,8,9,10,11}) );

        G.backtrack( std::move(G_orig) );
        
        CHECK( G.to_str() == "(0,7) (0,10); (2,4); (3,9); (4,9); (5,3); (6,1) (6,7); (8,2) (8,5); (11,1)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5,6,7,8,9,10,11}) );
    }
    #endif

    G.merge_verts(5,10); //merge 5,10 and 4,11
    state_str_stack.push_back( 
        std::pair<graph_repr, std::string>(G.get_state(), "(0,5) (0,7); (2,4); (3,9); (4,1) (4,9); (5,3); (6,1) (6,7); (8,2) (8,5)")
        );
    CHECK( G.get_no_e() == 11 );
    CHECK( G.to_str() == state_str_stack.back().second );
    CHECK( G.assert_data_structs() );
    CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5,6,7,8,9}) );
    
    G.merge_verts(1,9); //merge 1,9 and 0,8
    state_str_stack.push_back( 
        std::pair<graph_repr, std::string>(G.get_state(), "(0,2) (0,5) (0,7); (2,4); (3,1); (4,1); (5,3); (6,1) (6,7)")
        );
    CHECK( G.get_no_e() == 9 );
    CHECK( G.to_str() == state_str_stack.back().second );
    CHECK( G.assert_data_structs() );
    CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5,6,7}) );

    SECTION("in-between vertex-deletion + backtracking") {
        graph_repr G_orig = G.get_state();

        G.remove_vert(6); //removes 6 and 7

        CHECK( G.to_str() == "(0,2) (0,5); (2,4); (3,1); (4,1); (5,3)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5}) );

        G.backtrack( std::move(G_orig) );

        CHECK( G.to_str() == "(0,2) (0,5) (0,7); (2,4); (3,1); (4,1); (5,3); (6,1) (6,7)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5,6,7}) );
    }

    SECTION("in-between vertex-deletion + backtracking") {
        graph_repr G_orig = G.get_state();

        G.remove_vert(1); //removes 0 and 1

        CHECK( G.to_str() == "(2,4); (5,3); (6,7)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({2,3,4,5,6,7}) );

        G.backtrack( std::move(G_orig) );
        
        CHECK( G.to_str() == "(0,2) (0,5) (0,7); (2,4); (3,1); (4,1); (5,3); (6,1) (6,7)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5,6,7}) );
    }

    G.merge_verts(1,6); //merge 1,6 and 0,7
    state_str_stack.push_back( 
        std::pair<graph_repr, std::string>(G.get_state(), "(0,2) (0,5); (1,0); (2,4); (3,1); (4,1); (5,3)")
        );
    CHECK( G.get_no_e() == 7 );
    CHECK( G.to_str() == state_str_stack.back().second );
    CHECK( G.assert_data_structs() );
    CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5}) );


    SECTION("in-between vertex-deletion + backtracking") {
        graph_repr G_orig = G.get_state();

        G.remove_vert(1); //removes 0 and 1

        CHECK( G.to_str() == "(2,4); (5,3)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({2,3,4,5}) );

        G.backtrack( std::move(G_orig) );
        
        CHECK( G.to_str() == "(0,2) (0,5); (1,0); (2,4); (3,1); (4,1); (5,3)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({0,1,2,3,4,5}) );
    }

    G.merge_verts(5,1); //merge 1,5 and 0,4
    state_str_stack.push_back( 
        std::pair<graph_repr, std::string>(G.get_state(), "(2,4); (3,5); (4,2) (4,5); (5,3) (5,4)")
        );
    CHECK( G.get_no_e() == 6 );
    CHECK( G.to_str() == state_str_stack.back().second );
    CHECK( G.assert_data_structs() );
    CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({2,3,4,5}) );

    SECTION("in-between vertex-deletion + backtracking") {
        graph_repr G_orig = G.get_state();

        G.remove_vert(2);

        CHECK( G.to_str() == "(4,5); (5,4)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({4,5}) );

        G.backtrack( std::move(G_orig) );
        
        CHECK( G.to_str() == "(2,4); (3,5); (4,2) (4,5); (5,3) (5,4)" );
        CHECK( G.assert_data_structs() );
        CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({2,3,4,5}) );
    }

    G.merge_verts(3,5); //merge 3 and 5
    state_str_stack.push_back( 
        std::pair<graph_repr, std::string>(G.get_state(), "(2,3); (3,2)")
        );
    CHECK( G.get_no_e() == 2 );
    CHECK( G.to_str() == state_str_stack.back().second );
    CHECK( G.assert_data_structs() );
    CHECK_EQ_SET( G.get_v_vector(), vec<var_t>({2,3}) );


    //check backtracking!
    for (auto &&p : state_str_stack)
    {
        G.backtrack( std::move(p.first) );
        //check string repr
        CHECK( p.second == G.to_str() );
        //check data struct validity
        CHECK( G.assert_data_structs() );
    }

}