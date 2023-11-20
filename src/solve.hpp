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

#pragma once
//offers function to solve (and parse) xnf instances and guessing path files

#include <set>

#include "misc.hpp"
#include "LA/lineral.hpp"

/**
 * @brief solves xnf using provided opts
 * 
 * @param xnf vector of vector representing list of xor-clauses to be solved -- only works for 2-XNFs so far!
 * @param opts options specifying update alg, timeout, inprocessing settings etc
 * @param s stats to put statistics into
 * 
 * @return ret exit code
 */
int solve(const vec< vec<lineral> >& xnf, const options& opts, stats& s);

stats solve(const vec< vec<lineral> >& xnf, const options& opts);

/**
 * @brief uses the preprocessing of the class impl_graph and returns a string representation of the resulting system
 * 
 * @param xnf vector of vector representing list of xor-clauses to be solved -- only works for 2-XNFs so far!
 * @param opts options specifying update alg, timeout, inprocessing settings etc
 * @param s stats to put statistics into
 * @return std::string representation of pre-processed xnf instance
 */
std::string preprocess(const vec< vec<lineral> >& xnf, const options& opts, stats& s);

/**
 * @brief write string to file
 * 
 * @param fname file to write to
 * @param out string to write
 */
void write_str(const std::string& fname, const std::string& out);

struct parsed_xnf {
    var_t num_vars;
    var_t num_cls;
    vec< vec<lineral> > cls;

    parsed_xnf(var_t _num_vars, var_t _num_cls, vec< vec<lineral> > _cls) : num_vars(_num_vars), num_cls(_num_cls), cls(_cls) {};
    parsed_xnf(const parsed_xnf& o) : num_vars(o.num_vars), num_cls(o.num_cls), cls(o.cls) {};
};

/**
 * @brief parses file with name fname
 * 
 * @param fname guessing path file name; each line contains one index
 * @return reordering of variables, s.t. lex is the correct order
 */
reordering parse_gp(const std::string& fname);

/**
 * @brief parses file with name fname
 * 
 * @param fname xnf-file name
 * @param reordering permutation of indices
 * @return parsed_Xnf parsed num-cls, num-vars and parsed xlits
 */
parsed_xnf parse_file_gp(const std::string& fname, const reordering& P);

/**
 * @brief parses file with name fname
 * 
 * @param fname xnf-file name
 * @return parsed_Xnf parsed num-cls, num-vars and parsed xlits
 */
parsed_xnf parse_file(const std::string& fname);

/**
 * @brief print parsed xcls to string
 * 
 * @param xclss vector of xcls (repr as vector of linerals)
 * @return std::string string repr of the xclauses
 */
std::string to_str(const vec< vec<lineral> >& xclss);

/**
 * @brief checks whether sol is a solution of given xcls
 * 
 * @param clss vector of xcls (repr as vector of linerals)
 * @param sol solution to be checked
 * @return true iff sol is a solution of the xclauses
 */
inline bool check_sol(const vec< vec<lineral> >& clss, const vec<bool>& sol) {
    return std::all_of( clss.begin(), clss.end(), /* all clauses need to be satisfied */
                    [&sol] (vec<lineral> xcls) -> bool { 
                        return std::any_of(xcls.begin(), xcls.end(), [&sol](lineral l) { return l.eval(sol); } ); /* at least one lit of clause must be satisfied */
                        }
                    );
}
