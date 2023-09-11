# 2xnf_solver

> Graph-based DPLL-SAT Solver for propositional logic formulas in XOR-OR-AND normal form (XNF) implemented in C++20.

### Usage

Run `2xnf_solver -h` to get the following help:

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
-pp --preprocess             	expects 'no', 'scc', 'fls_scc', or 'fls_scc_ee' [default: "fls_scc"]
-ppo --preprocess-out        	fname for pre-processed input xnf; output xnf is equivalent to input xnf [default: ""]
-vb --verb                   	verbosity (choose in 0-100) [default: 0]
-t --time-out                	timeout in seconds (negative to deactivate) [default: -1]
```

Encode your 2-XNF problem in a DIMACS-like format with header `p xnf n_vars n_cls` and where each clause is in a new line where linerals, a XOR of literals, are encoded as literals connected by `+` and the clause terminates with `0`.

The XNF clause $(\neg X_1 \oplus X_2 \oplus X_3) \vee (X_4\oplus X_5)$ is
encoded as `-1+2+3 4+5 0`.


### Build

On Ubuntu/Debian ensure that you have installed the packages `build-essential`, `cmake`, and `libm4ri-dev` (optionally `libjemalloc-dev`).
Then run `cmake .` and `make 2xnf_solver`. When the build is finished you will have the executible `2xnf_solver`.

If this does not work for you or you are running a different OS, you can also use the provided `Dockerfile` to build your own docker-image by `docker build -t 2xnf_solver:1.0 .` and use the solver through `docker_2xnf_solver`.

