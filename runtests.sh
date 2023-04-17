#!/bin/bash

passes="passperms.txt"
tests="testpaths.txt"
results="results.txt"

while read p; do
    while read t; do
        clang-14 -S -c -emit-llvm ${t} -o ${t}.bc
        opt-14 -enable-new-pm=0 -mem2reg ${p} ${t}.bc -S -o ${t}.ll
        echo ${p} >> "$results"
        echo ${t} >> "$results"
        (time clang-14 ${t}.ll)2>"$results"
        (time ./a.out)2>"$results"
        du a.out
    done < "$tests"
done < "$passes"