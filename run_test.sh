#!/bin/bash

opts="-loop-unroll"

cd benchmarks/gcc_loops
make clean
make OPTFLAGS="$opts"
make run > "gcc_loops.csv"

cd ../lcals
make clean
make OPTFLAGS="$opts"
make run OUTDIR="lcals_results"

cd ../NPB-SER
make clean
make suite OPTFLAGS="$opts"
make run > "npb_results.txt"

cd ../../
python3 parse_results.py "$opts" "benchmarks/gcc_loops/gcc_loops.csv" "benchmarks/lcals/lcals_results" "benchmarks/NPB-SER/npb_results.txt"

