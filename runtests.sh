#!/bin/bash

passes="passperms.txt"
tests="testpaths.txt"

while read p; do
    while read t; do
        clang-14 -S -c -emit-llvm ${t} -o ${t}.bc
        opt-14 -enable-new-pm=0 -mem2reg ${p} ${t}.bc -S -o ${t}.ll
        time clang-14 ${t}.ll
        time ./a.out
        du a.out
    done < "$tests"
done < "$passes"