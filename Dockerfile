#dockerfile to build the solver and run benchmarks on it
FROM debian:trixie-slim
COPY . /2xornado
RUN apt update && apt install --assume-yes --no-install-recommends build-essential cmake pkg-config libm4ri-dev libjemalloc-dev libbenchmark-dev libcatch2-dev && cd 2xornado && cmake . && make
ENTRYPOINT [ "2xornado/2xornado" ]

# build image with
#   docker build -t 2xornado:1.0 .
# solve instance at path/fname with
#   docker run -v path/fname:path/fname 2xornado_docker path/fname [args]
