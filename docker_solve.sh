#!/bin/bash

#run docker solver on file
docker run -v $(realpath $1):$(realpath $1) 2xnf_solver:1.0 $(realpath $1) ${@:2}
