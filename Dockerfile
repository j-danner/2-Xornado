#dockerfile to build the solver and run benchmarks on it
FROM debian:bullseye-slim
COPY . /2-Xornado
RUN apt update && apt install --assume-yes --no-install-recommends build-essential cmake libm4ri-dev libjemalloc-dev && cd 2xornado && cmake . && make 2xornado
ENTRYPOINT [ "2-Xornado/2xornado" ]

# build image with
#   docker build -t 2xornado:1.0 .
# solve instance at path/fname with
#   docker run -v path/fname:path/fname 2xornado_docker path/fname [args]
