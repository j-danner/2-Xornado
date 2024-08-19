# 2-Xornado

> Graph-based DPLL-SAT Solver for propositional logic formulas in XOR-OR-AND normal form (XNF) implemented in C++20.

### Usage

Run `2xornado -h` to get the following help:

```
Positional arguments:
fname                        	path to 2xnf-instance

Optional arguments:
-h --help                    	shows help message and exits [default: false]
-v --version                 	shows version and compilation information [default: false]
-dh --decision-heuristic     	decision heuristic; 'mp' for MaxPath, 'mr' for MaxReach, 'mbn' for MaxBottleNeck, 'fv' for FirstVert [default: "mp"]
-fls --failed-lineral-search 	failed lineral search; 'no' to deactivate, 'trivial' to only search for trivial, 'full' to search for all failed linerals. [default: "no"]
-flss --fls-schedule         	number n s.t. every n-th crGCP we perform fls [default: 1]
-sc --score                  	activate weighting of vars based on score (inspired by VSIDS) [default: false]
-simple                      	construct the trivial IGS from the input 2-XNF instead of the extended trivial IGS [default: false]
-pp --preprocess             	expects 'no', 'scc', 'fls_scc' (failed lineral search + SCC), or 'fls_scc_ee' (failed lineral search + SCC + edge extension) [default: "fls_scc"]
-ppo --preprocess-out        	path for output of xnf after pre-processing (input and output xnf are equivalent)
-gp --guessing-path          	path to file storing guessing path; each line contains exactly one number corr to the corresponding variable; USE WITH CAUTION!
-vb --verb                   	verbosity (choose in 0-100) [default: 0]
-t --time-out                	timeout in seconds (negative to deactivate) [default: -1]
```


#### XNF Format

Encode your 2-XNF problem in a DIMACS-like format with header `p xnf n_vars n_cls` and where each clause is in a new line where linerals, a XOR of literals, are encoded as literals connected by `+` and the clause terminates with `0`.

The XNF clause $(\neg X_1 \oplus X_2 \oplus X_3) \vee (X_4\oplus X_5)$ is
encoded as `-1+2+3 4+5 0`.

#### ANF Input

To solve systems of (quadratic) polynomials, use our `anf_to_2xnf` conversion tool from [here](https://github.com/Wrazlmumfp/anf_to_2xnf.git) to encode the polynomials as an instance in 2-XNF.

### Build

On Ubuntu/Debian ensure that you have installed the packages `build-essential`, `cmake`, `libm4ri-dev`, and `libtbb-dev` (optionally `libjemalloc-dev`, `libbenchmark-dev`, `catch2`).
Then run `cmake .` and `make 2xornado`. When the build is finished you will have the executable `2xornado`.

If this does not work for you or you are running a different OS, you can also use the docker image `jdanner/2xornado` (download via `docker pull jdanner/2xornado:latest`) or build it yourself using the provided `Dockerfile` by `docker build -t jdanner/2xornado:latest .`. Then access 2-Xornado through `docker_2xornado`.

