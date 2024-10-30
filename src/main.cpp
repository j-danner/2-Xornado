// Copyright (c) 2022-2024 Julian Danner <julian.danner@uni-passau.de>
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

#include "solve.hpp"

#include "argparse/argparse.hpp"


//main -- parses args
int main(int argc, char const *argv[])
{
    argparse::ArgumentParser program(__PROJECT_NAME, __VERSION__, argparse::default_arguments::help);
    program.add_argument("-v", "--version")
      .action([=]([[maybe_unused]] const std::string& s) {
        std::stringstream out;
        out << "c 2xornado created by J. Danner (2022-2024)" << std::endl;
        out << "c version:           " << __PROJECT_VERSION << std::endl;
        out << "c compilation date:  " << __DATE__ << " at " << __TIME__ << std::endl;
        out << "c compiler:          " << __CMAKE_CXX_COMPILER_ID << " " << __CMAKE_CXX_COMPILER_VERSION << " using C++" << __CMAKE_CXX_STANDARD << std::endl;
        out << "c compilation flags:" << __CMAKE_CXX_FLAGS << std::endl;
        #ifdef FULL_REDUCTION
            out << "c crGCP uses NR for updating vertices" << std::endl;
        #else
            out << "c crGCP only reduces LTs of vertices" << std::endl;
        #endif
        #ifdef USE_TRIE
            out << "c using trie for mapping vertices to linerals" << std::endl;
        #else
            out << "c using hashmaps for mapping vertices to linerals" << std::endl;
        #endif
        #ifdef USE_LHGR
            out << "c using lhgr for graph representation" << std::endl;
        #else
            out << "c using AL for graph representation" << std::endl;
        #endif
        std::cout << out.str();
        std::exit(0);
      })
      .default_value(false)
      .help("shows version and compilation information")
      .implicit_value(true)
      .nargs(0);
    
    //add args:
    //fname
    program.add_argument("fname")
        .help("path to 2xnf-instance");
    //dec_heu
    program.add_argument("-dh","--decision-heuristic")
        .help("decision heuristic; 'mp' for MaxPath, 'mr' for MaxReach, 'mbn' for MaxBottleNeck, 'fv' for FirstVert")
        .default_value(std::string("mp"))
        .action([](const std::string& value) {
            static const vec<std::string> choices = { "fv", "mp", "mr", "mbn" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            //arg invalid!
            throw std::runtime_error("invalid argument passed for parameter -dh");
        });

    //fls opts
    program.add_argument("-fls","--failed-lineral-search")
        .help("failed lineral search; 'no' to deactivate, 'trivial' to only search for trivial, 'full' to search for all failed linerals.")
        .default_value(std::string("no"))
        .action([](const std::string& value) {
            static const vec<std::string> choices = { "no", "trivial", "trivial_cc", "full" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            //arg invalid!
            throw std::runtime_error("invalid argument passed for parameter -fls");
        });

    //fls_schedule
    program.add_argument("-flss","--fls-schedule")
        .help("number n s.t. every n-th crGCP we perform fls")
        .default_value(1)
        .scan<'i', int>();
        
    
    //upd opts
    //program.add_argument("-upd","--update-alg")
    //    .help("algorithm to use for update-graph function, 'ts' for alg in two steps (1. update all linerals, 2. merge verts); 'hf' for hash-fight based update; 'par' for parallel version.")
    //    .default_value(std::string("ts"))
    //    .action([](const std::string& value) {
    //        static const vec<std::string> choices = { "ts", "hf", "par", "hfd" };
    //        if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
    //            return value;
    //        }
    //        //arg invalid!
    //        throw std::runtime_error("invalid argument passed for parameter -upd");
    //    });
    
    //score opts
    program.add_argument("-sc","--score")
        .help("activate weighting of vars based on score (inspired by VSIDS)")
        .default_value(false)
        .implicit_value(true);
    
    //simple trivial graph
    program.add_argument("-simple")
        .help("construct the trivial IGS from the input 2-XNF instead of the extended trivial IGS")
        .default_value(false)
        .implicit_value(true);
    
    //score opts
    program.add_argument("-pp","--preprocess")
        .help("expects 'no', 'scc', 'fls_scc' (failed lineral search + SCC), or 'fls_scc_ee' (failed lineral search + SCC + edge extension)")
        .default_value(std::string("fls_scc"))
        .action([](const std::string& value) {
            static const vec<std::string> choices = { "no","scc","fls_scc","fls_scc_ee" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            //arg invalid!
            throw std::runtime_error("invalid argument passed for parameter -pp");
        });
    
    program.add_argument("-ppo","--preprocess-out")
        .help("path for output of xnf after pre-processing (input and output xnf are equivalent)");

            
    //guessing path input
    program.add_argument("-gp","--guessing-path")
        .help("path to file storing guessing path; each line contains exactly one number corr to the corresponding variable; USE WITH CAUTION!");

    
    //verbosity
    #ifdef VERBOSITY
        program.add_argument("-vb", "--verb")
            .help("verbosity (choose in 0-100)")
            .default_value(0)
            .scan<'i', int>();
    #endif
    
    //timeout
    program.add_argument("-t","--time-out")
        .help("timeout in seconds (negative to deactivate)")
        .default_value(-1)
        .scan<'i', int>();

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    //parse string-input to 
    auto fname = program.get<std::string>("fname");

    auto dh_str = program.get<std::string>("-dh");
    dec_heu dh = dec_heu::mp;
    if(dh_str=="fv") dh = dec_heu::fv;
    else if(dh_str=="mp") dh = dec_heu::mp;
    else if(dh_str=="mr") dh = dec_heu::mr;
    else if(dh_str=="mbn") dh = dec_heu::mbn;
    else if(dh_str=="lex") dh = dec_heu::lex;
    
    auto fls_str = program.get<std::string>("-fls");
    fls_alg fls = fls_alg::no;
    if(fls_str=="no") fls = fls_alg::no;
    else if(fls_str=="trivial") fls = fls_alg::trivial;
    else if(fls_str=="trivial_cc") fls = fls_alg::trivial_cc;
    else if(fls_str=="full") fls = fls_alg::full;

    auto fls_s = program.get<int>("-flss");
    
    //auto upd_str = program.get<std::string>("-upd");
    upd_alg upd = upd_alg::ts;
    //if(upd_str=="ts") upd = upd_alg::ts;
    //else if(upd_str=="hf") upd = upd_alg::hf;
    //else if(upd_str=="par") upd = upd_alg::par;
    //else if(upd_str=="hfd") upd = upd_alg::hfd;
    
    sc score = program.get<bool>("-sc") ? sc::active : sc::inactive;

    constr ext = program.get<bool>("-simple") ? constr::simple : constr::extended;
    
    auto pp_str = program.get<std::string>("-pp");
    preproc pp = preproc::fls_scc;
    if(pp_str=="no") pp = preproc::no;
    else if(pp_str=="scc") pp = preproc::scc;
    else if(pp_str=="fls_scc") pp = preproc::fls_scc;
    else if(pp_str=="fls_scc_ee") pp = preproc::fls_scc_ee;

    const bool only_preprocess = program.is_used("-ppo");
    const std::string pp_out = only_preprocess ? program.get<std::string>("-ppo") : "";
        
    const std::string gp_fname = program.is_used("-gp") ? program.get<std::string>("-gp") : "";
    if(program.is_used("-gp")) dh = dec_heu::lex;


    #ifdef VERBOSITY
        int verb = program.get<int>("-vb");
    #else
        int verb = 0;
    #endif
    
    auto time_out = program.get<int>("-t");

    stats s;
    s.begin = std::chrono::steady_clock::now();

    //parse file
    try {
        reordering P = parse_gp( gp_fname );
        parsed_xnf p_xnf = P.size()==0 ? parse_file( fname ) : parse_file_gp( fname, P );


        //init options
        options opts( p_xnf.num_vars, p_xnf.num_cls, dh, fls, fls_s, upd, score, ext, pp, verb, time_out, P );


        if(only_preprocess) {
            std::string out = preprocess(p_xnf.cls, opts, s);
            if(out.size()>0) {
                write_str(pp_out, out);
                return 0;
            }
            return 1; //preprocessing failed.
        }

        int ret = solve(p_xnf.cls, opts, s);
        if(ret != 0) throw std::runtime_error("solving failed!");
        //print stats
        s.end = std::chrono::steady_clock::now();
        s.print_final();

        if(s.finished && s.sat) { //check sol!
            const bool is_sol = check_sol(p_xnf.cls, s.sol);
            std::cout << "c solution " << ( is_sol ? "verified" : "INCORRECT!") << std::endl;
            return is_sol ? 0 : -1;
        }
        
        return 0;
    } catch (std::exception &ex) {
        std::cout << "s INDEFINITE" << std::endl;
        return 1;
    }
}
