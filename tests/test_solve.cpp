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

//file to test implementation of solve
#include "../src/solve.hpp"

#include <catch2/catch_all.hpp>


TEST_CASE( "solving with different options" , "[impl-graph][graph][parser][solve]" ) {
    auto fname = GENERATE("tests/2xnfs/ToyExample-type1-n10-seed1.xnf", "tests/2xnfs/ToyExample-type1-n10-seed0.xnf", "tests/2xnfs/rand-3-6.xnf", "tests/2xnfs/rand-20-40.xnf");
    //auto fname = "tests/2xnfs/rand-10-20.xnf";
    //auto fname = "tests/2xnfs/rand-3-6.xnf";
    //auto fname = "tests/2xnfs/ToyExample-type1-n10-seed0.xnf";
    auto clss = parse_file( fname );
    auto xnf = clss.cls;
    var_t num_vars = clss.num_vars;
    var_t num_cls = clss.num_cls;

    SECTION( "dh:fv-fls:no-upd:ts" ) {
        options opts(num_vars, num_cls, dec_heu::fv, fls_alg::no, upd_alg::ts, 1, 0, 0);
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "dh:mr-fls:no-upd:ts" ) {
        options opts(num_vars, num_cls, dec_heu::mr, fls_alg::no, upd_alg::ts, 1, 0, 0);
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "dh:mr-fls:full-upd:hf" ) {
        options opts(num_vars, num_cls, dec_heu::mr, fls_alg::full, upd_alg::hf, 1, 0, 0);
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "dh:mp-fls:trivial-upd:hf" ) {
        options opts(num_vars, num_cls, dec_heu::mp, fls_alg::trivial, upd_alg::hf, 1, 0, 0);
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "dh:msp-fls:trivial-upd:hfd" ) {
        options opts(num_vars, num_cls, dec_heu::mp, fls_alg::trivial, upd_alg::hfd, 1, 0, 0);
        opts.score = sc::active;
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }

    SECTION( "dh:mbn-fls:trivial-upd:hfd" ) {
        options opts(num_vars, num_cls, dec_heu::mbn, fls_alg::trivial, upd_alg::hfd, 1, 0, 0);
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "dh:mp-fls:trivial-upd:hf -- 4 jobs" ) {
        options opts(num_vars, num_cls, dec_heu::mp, fls_alg::trivial, upd_alg::hf, 4, 0, 0);
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "dh:mp-fls:trivial-upd:par -- 4 jobs" ) {
        options opts(num_vars, num_cls, dec_heu::mp, fls_alg::trivial, upd_alg::par, 4, 0, 0);
        stats s = solve(xnf, opts);
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }
    
    SECTION( "dh:mp-fls:trivial-upd:hf -- terminate within timeout" ) {
        options opts(num_vars, num_cls, dec_heu::mp, fls_alg::trivial, upd_alg::hf, 1, 0, 0);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        stats s = solve(xnf, opts);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        float time_no_time_out = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0f;

        options opts2(num_vars, num_cls, dec_heu::mp, fls_alg::trivial, upd_alg::hf, 1, 0, time_no_time_out+3);
        begin = std::chrono::steady_clock::now();
        stats s2 = solve(xnf, opts2);
        end = std::chrono::steady_clock::now();
        float time_with_time_out = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0f;
        CHECK( time_no_time_out == Catch::Approx(time_with_time_out).margin(0.1) );

        CHECK( s.finished == true ); //SAT
        CHECK( s.sat == true ); //SAT
        CHECK( check_sol(clss.cls, s.sol) );
    }
}

TEST_CASE("solving with different options -- timeout", "[impl-graph][graph][parser][solve]")  {
    auto clss = parse_file("tests/2xnfs/ToyExample-type1-n20-seed0.xnf");
    auto xnf = clss.cls;
    auto num_vars = clss.num_vars;
    auto num_cls = clss.num_cls;

    options opts(num_vars, num_cls, dec_heu::mp, fls_alg::trivial, upd_alg::hf, 1, 0, 1);
    stats s = solve(xnf, opts);

    CHECK( s.finished == false );
}
