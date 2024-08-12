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


#include "../src/solve.hpp"
#include "../src/impl_graph.hpp"

#include <benchmark/benchmark.h>

#define concat(first, second) first second


static void BM_dpll_solve(benchmark::State& state, std::string fname) {
    for (auto _ : state) {
        auto clss = parse_file(fname);
        auto IG = impl_graph(clss);
        stats s = IG.dpll_solve();
    }
}

BENCHMARK_CAPTURE(BM_dpll_solve, rand-20-60, concat(BENCH_FILES, "/2xnfs/rand-20-60.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);
BENCHMARK_CAPTURE(BM_dpll_solve, rand-10-20, concat(BENCH_FILES, "/2xnfs/rand-10-20.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);
BENCHMARK_CAPTURE(BM_dpll_solve, rand-5-10,  concat(BENCH_FILES, "/2xnfs/rand-5-10.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);
BENCHMARK_CAPTURE(BM_dpll_solve, mq-toyexample-type1-n15,  concat(BENCH_FILES, "/2xnfs/ToyExample-type1-n15-seed0.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);
BENCHMARK_CAPTURE(BM_dpll_solve, mq-toyexample-type1-n15,  concat(BENCH_FILES, "/2xnfs/ToyExample-type1-n15-seed1.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);
BENCHMARK_CAPTURE(BM_dpll_solve, mq-toyexample-type1-n15,  concat(BENCH_FILES, "/2xnfs/ToyExample-type1-n15-seed2.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);
BENCHMARK_CAPTURE(BM_dpll_solve, mq-toyexample-type1-n15,  concat(BENCH_FILES, "/2xnfs/ToyExample-type1-n15-seed3.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);
BENCHMARK_CAPTURE(BM_dpll_solve, mq-toyexample-type1-n15,  concat(BENCH_FILES, "/2xnfs/ToyExample-type1-n15-seed4.xnf") )->Unit(benchmark::kMillisecond)->MinTime(2);

int xlit_performance(var_t n, long k) {
    //compute k random lineral additions in n vars
    vec< lineral > linerals;
    linerals.reserve(2*k);

    vec<var_t> xlit_set;
    srand((unsigned)time(NULL));
    for(int j=0; j<2*k; j++) {
        for (int i=0; i < n; i++){
            if(rand() % 2) xlit_set.push_back(i);
        }
        linerals.push_back( lineral(xlit_set) );
        xlit_set.clear();
    }

    //performance analysis:
    auto start = std::chrono::steady_clock::now();
    vec<lineral> sums;
    sums.reserve(k);

    for (int i = 0; i < k; i++)
    {
        sums[i] = linerals[2*i] + linerals[2*i+1];
    }
    auto end = std::chrono::steady_clock::now();

    std::cout << k << " additions of random linerals in " << n << " inds took " << std::chrono::duration_cast<std::chrono::seconds> (end - start).count() << "s." << std::endl;

    return 1;
}

//    //xlit_performance(1000, 1000, 0.5);
//    //xlit_performance(1000, 1<<20);
//    //xlit_performance(100, 1<<15, 0.5);
//    //xlit_performance(1000, 1<<15);
//    //xlit_performance(1000, 1<<15, 0.5);
//    //xlit_performance(1000, 1<<15, 0.2);
//    ////xlit_performance(1000, 1<<15, 0.1);
//    //xlit_performance(1000, 1<<15, 0.01);
//    //xlit_performance(1000, 1<<15, 0.001);




BENCHMARK_MAIN();
