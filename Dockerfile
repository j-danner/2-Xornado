#dockerfile to build the solver and run benchmarks on it
FROM debian:bullseye-slim
COPY . /2xnf_solver
RUN apt update && apt install --assume-yes --no-install-recommends build-essential cmake libm4ri-dev libjemalloc-dev && cd 2xnf_solver && cmake . && make 2xnf_solver
ENTRYPOINT [ "2xnf_solver/2xnf_solver" ]

# build image with
#   docker build -t 2xnf_solver:1.0 .
# solve instance at path/fname with
#   docker run -v path/fname:path/fname 2xnf_solver_docker path/fname [args]
