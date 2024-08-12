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

//file to test implementation of xsys
#include "../src/LA/lineqs.hpp"

#include <catch2/catch_all.hpp>


TEST_CASE( "xsys creation/reduction/addition", "[xsys]" ) {
    auto linerals = vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({1,2})), lineral(vec<var_t>({0,2}))});
    LinEqs ls = LinEqs( linerals );
    CHECK( ls.to_str() == "x1+1 x2+1" );
    
    linerals = vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({0,1})), lineral(vec<var_t>({0,1}))});
    ls = LinEqs( linerals );
    CHECK( ls.to_str() == "x1+1" );
    
    linerals = vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({1,2})), lineral(vec<var_t>({0,2}))});
    ls = LinEqs( linerals );
    CHECK( ls.to_str() == "x1+1 x2+1" );
    
    
    linerals = vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({0,1})), lineral(vec<var_t>({0,1}))});
    ls = LinEqs( linerals );
    CHECK( ls.to_str() == "x1+1" );
    
    linerals = vec<lineral>({lineral(vec<var_t>({0,1})), lineral(vec<var_t>({1,2})), lineral(vec<var_t>({2,3})), lineral(vec<var_t>({3,4})), lineral(vec<var_t>({4,5}))});
    ls = LinEqs( linerals );
    CHECK( ls.to_str() == "x1+1 x2+1 x3+1 x4+1 x5+1" );
    CHECK( ls.reduce( lineral(vec<var_t>({0})) ).to_str() == "1" );
    CHECK( ls.reduce( lineral(vec<var_t>({1})) ).to_str() == "1" );
    CHECK( ls.reduce( lineral(vec<var_t>({0,1})) ).to_str() == "0" );
    CHECK( ls.reduce( lineral(vec<var_t>({6})) ).to_str() == "x6" );
    
    linerals = vec<lineral>({lineral(vec<var_t>({0,1,2,3})), lineral(vec<var_t>({1,2,3,5})), lineral(vec<var_t>({3,4})), lineral(vec<var_t>({0,4})), lineral(vec<var_t>({0,5,6})) });
    ls = LinEqs( linerals );
    CHECK( linerals[0].to_str() == "x1+x2+x3+1" );
    CHECK( linerals[1].to_str() == "x1+x2+x3+x5" );
    CHECK( linerals[2].to_str() == "x3+x4" );
    CHECK( linerals[3].to_str() == "x4+1" );
    CHECK( linerals[4].to_str() == "x5+x6+1" );

    CHECK( ls.to_str() == "x1+x2 x3+1 x4+1 x5+1 x6" );
    CHECK( ls.reduce( lineral(vec<var_t>({2})) ).to_str() == "x2" );
    CHECK( ls.reduce( lineral(vec<var_t>({2,3,0})) ).to_str() == "x2" );
    CHECK( ls.reduce( lineral(vec<var_t>({2,3,4,6,0})) ).to_str() == "x2+1" );
    vec<bool> sol(6, false);
    ls.solve( sol );
    CHECK( sol == vec<bool>({0,0,1,1,1,0}) );
    vec<bool> sol_(6, true);
    ls.solve( sol_ );
    CHECK( sol_ == vec<bool>({1,1,1,1,1,0}) );

    CHECK( linerals[0].to_str() == "x1+x2+x3+1" );
    CHECK( linerals[1].to_str() == "x1+x2+x3+x5" );
    CHECK( linerals[2].to_str() == "x3+x4" );
    CHECK( linerals[3].to_str() == "x4+1" );
    CHECK( linerals[4].to_str() == "x5+x6+1" );
    
    std::reverse(linerals.begin(), linerals.end());
    REQUIRE( linerals[4].to_str() == "x1+x2+x3+1" );
    REQUIRE( linerals[3].to_str() == "x1+x2+x3+x5" );
    REQUIRE( linerals[2].to_str() == "x3+x4" );
    REQUIRE( linerals[1].to_str() == "x4+1" );
    REQUIRE( linerals[0].to_str() == "x5+x6+1" );

    ls = LinEqs( linerals );

    REQUIRE( linerals[4].to_str() == "x1+x2+x3+1" );
    REQUIRE( linerals[3].to_str() == "x1+x2+x3+x5" );
    REQUIRE( linerals[2].to_str() == "x3+x4" );
    REQUIRE( linerals[1].to_str() == "x4+1" );
    REQUIRE( linerals[0].to_str() == "x5+x6+1" );

    CHECK( ls.to_str() == "x1+x2 x3+1 x4+1 x5+1 x6" );
    CHECK( ls.reduce( lineral(vec<var_t>({2})) ).to_str() == "x2" );
    CHECK( ls.reduce( lineral(vec<var_t>({2,3,0})) ).to_str() == "x2" );
    CHECK( ls.reduce( lineral(vec<var_t>({2,3,4,6,0})) ).to_str() == "x2+1" );


    //test operator +=
    linerals = vec<lineral>({lineral(vec<var_t>({0,1,2,3})), lineral(vec<var_t>({1,2,3,5}))});
    auto xlits2 = vec<lineral>({lineral(vec<var_t>({3,4})), lineral(vec<var_t>({0,4})), lineral(vec<var_t>({0,5,6})) });

    auto sys1 = LinEqs(linerals);
    auto sys2 = LinEqs(xlits2);

    REQUIRE( sys1.to_str() == "x1+x2+x3+1 x5+1");
    REQUIRE( sys2.to_str() == "x3+1 x4+1 x5+x6+1");

    sys1 += sys2;
    CHECK( sys1.to_str() == "x1+x2 x3+1 x4+1 x5+1 x6" );

    //operator +
    sys1 = LinEqs(linerals);
    sys2 = LinEqs(xlits2);

    REQUIRE( sys1.to_str() == "x1+x2+x3+1 x5+1");
    REQUIRE( sys2.to_str() == "x3+1 x4+1 x5+x6+1");

    auto sys3 = sys1+sys2;
    CHECK( sys3.to_str() == "x1+x2 x3+1 x4+1 x5+1 x6" );
}