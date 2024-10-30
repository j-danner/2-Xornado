# 2-Xornado

> Graph-based DPLL-SAT Solver for propositional logic formulas in 2-XNF implemented in C++20.

This repository contains the source code of our implementation of Algorithm 8 (`G_2XNF_DPLL`) from our paper *SAT Solving Using XOR-OR-AND Normal Forms* called `2xornado`.

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


#### Input Format

###### 2-XNF

Informally, a XNF-formula is a CNF-formula, where literals are replaced by XORs of literals, called *linerals*.
The solver reads 2-XNF-formulas in a DIMACS-like format with header `p xnf n_vars n_cls`, where linerals, a XOR of literals, are encoded as literals connected by `+`; and clauses are terminated by `0`.
For example, the XNF clause $(\neg X_1 \oplus X_2 \oplus X_3) \vee (X_4\oplus X_5)$ is encoded as `-1+2+3 4+5 0`.

###### 2-CNF-XOR

A CNF-XOR formula consists of a CNF formula and XOR constraints on the variables. `2xornado` supports these constraints in the usual encoding as a line starting with `x` followed by the involved literals. For example the constraint $\neg X_1 \oplus X_2$ can be encoded as `x -1 2 0`, but also as the XNF unit clause `-1+2 0`.

###### ANF

To solve systems in algebraic normal form (ANF), use our [`anf_to_2xnf`](https://github.com/Wrazlmumfp/anf_to_2xnf.git) tool to encode the system of polynomials in 2-XNF first.


### Build

On Ubuntu/Debian ensure that the packages `build-essential`, `cmake`, and `libm4ri-dev` (optionally `libjemalloc-dev`, `libbenchmark-dev`, `catch2`/`libcatch2-dev`) are installed. On Arch-based systems the packages `base-devel`, `cmake`, and `m4ri` (optionally `jemalloc`, `benchmark`, `catch2`) are required.
Then run `cmake .` and `make 2xornado`.

Alternatively, use the docker image `jdanner/2xornado` (download via `docker pull jdanner/2xornado:latest`) or build it yourself using the provided [`Dockerfile`](Dockerfile) by `docker build -t jdanner/2xornado:latest .`. Then access 2-Xornado through [`docker_2xornado`](docker_2xornado).

