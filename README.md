# CompilersFinal

### Running Tests
The `runtests.sh <start> <end>` file will run the benchmarks on the permutations start through end (inclusive). Permutation $i$ is the $i$ th line of `permutation.txt` (1-indexed). If permuation $i$ is $(i_1, i_2, i_3, i_4)$, then the optimizations run are those on lines $i_1, i_2, i_3$, and $i_4$ in `passes.txt`. Results are saved in the `results` directory and are named according to the permuation number. 

For example, the 0th permuation is (0, 1, 2, 3) and will run the below optimizations (in addition to the base optimizations):
```
-loop-deletion -loop-extract -loop-extract-single -loop-reduce
```

The followng optimizations are always run:
```
-mem2reg -loop-simplify -lcssa -loop-rotate -licm
```

### SLURM
To run using SLURM, use `submit_jobs.sh <start> <end>`. This will batch submit SLURM jobs to run permutations <start> through <end> inclusively where each job runs 91 permutations. For example, `submit_jobs.sh 1 200` will submit three jobs, 1 for permutations 1 through 91, 1 for permutations 92 through 182, and 1 for permutations 183 through 200. 

### Running Individually
To run the GCC loops, TSVC, or LCALS benchmark, use the following as a guideline:
```
# GCC Loops
cd benchmarks/gcc_loops   # GCC loops
cd benchmarks/TSVC_2/src  # TSVC
cd benchmarks/lcals       # LCALS

make OPTFLAGS="{passes to run}"
make size
make run
```
To run the NPB benchmark, use the following as a guideline
```
cd benchmarks/NPB-SER

make suite OPTFLAGS="{passes to run}"
make size
make run
```
