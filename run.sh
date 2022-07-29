#!/bin/bash

make clean && make

data=$1

K=$2

./bin/degeneracy_cliques -i $data -t A -d 0 -k $K