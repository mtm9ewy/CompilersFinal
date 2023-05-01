# CompilersFinal

## Approach
We identified 16 loop optimizations passes, which give rise to 16! permutations. As such, it is infeasible to test all possible orderings. Instead, we use an iterative approach. We first test all permuations of length 4 of the 16 passes, which yields 43,680 unique passes. Then we take the top five results and generate all permuations of length 4 with the remaning 12 passes, which yields $5 \cdot 11,880 = 59,400$ unique passe. This process is continued


\# Results from Previous Iteration | # Remaining Passes | # Unique passes
----------------------------------|--------------------|------------------
0 | 16 |  43,680
5 | 12 | 59,40
25 | 8 | 42,000
1,000 | 4 | 24,000


After this, we take the top result out of the $1 + 1 + 1 + 24,000 = 24,003$ choices.

## Reproducing Results
We use the cortado machines via SLURM for all tests.

### Generating Optimization Pass Permutations
To generate the permutations for $n$ remaining passes, run 
```
python3 create_optimizations_list <path to file containing best optimizations from previous iteration>
```
For the case of 16 remaning passes, there is no previous iteration, so pass an empty file. Results are saved as `optimization_permutations_{n}.txt`.

## Running Permuations via SLURM
The `runtests.sh <start> <end>` script will run the benchmarks on the permutations start through end (inclusive). Permutation $i$ is the $i$ th line of `optimization_permutations_{n}.txt` (1-indexed). Note that `run_test.sh` must be manually modified to use the correct `optimization_permutations_{n}.txt` file.  


To run using SLURM, use `submit_jobs.sh <start> <end>`. This will batch submit SLURM jobs to run permutations <start> through <end> inclusively where each job runs a fixed number of permuations. Again, this file must be manually modified to set the correct bounds. 

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

## Evaluating Permuations
All results are saved in the `temp_results` folder. Run `python3 aggregate_results.py` to aggregate them into a single `results.csv` file. Then manually modify `score_results.py` to take the top $n$ results along with the proper file pathes for the results file and the `optimization_permutations_{n}.txt` file. Run `python3 score_results.py` and the top permutations will be saved into the filename set via `best_file`. 
