#!/bin/bash

permutation_number="$1"
GCC_DIR="$2"
LCALS_DIR="$3"
NPB_DIR="$4"
TSVC_DIR="$5"

# Get the permutation on nth line from permuations.txt (1 indexed)
permutation=$(sed "${permutation_number}q;d" permutations.txt)
# Get the passes corresonding to this permutation
opts=$(python3 get_optimizations.py "$permutation")

# Run the GCC Loops benchmark
cd "benchmarks/${GCC_DIR}"
make clean
{ time make OPTFLAGS="$opts"; } 2> "../../results/${permutation_number}_gcc_compile.txt"
make run > "gcc_loops_${permutation_number}.csv"
make size > "../../results/${permutation_number}_gcc_size.txt"


# Run the LCALS benchmark
cd "../${LCALS_DIR}"
make clean
{ time make OPTFLAGS="$opts"; } 2> "../../results/${permutation_number}_lcals_compile.txt"
make run OUTDIR="lcals_results_${permutation_number}"
make size > "../../results/${permutation_number}_lcals_size.txt"


# Run the NPB benchmark
cd "../${NPB_DIR}"
make clean
{ time make suite OPTFLAGS="$opts"; } 2> "../../results/${permutation_number}_npb_compile.txt"
make run > "npb_results_${permutation_number}.txt"
make size > "../../results/${permutation_number}_npb_size.txt"


# Run the TSVC benchmark
cd "../${TSVC_DIR}"
make clean
{ time make OPTFLAGS="$opts"; } 2>"../../results/${permutation_number}_tsvc_compile.txt"
make run > "tsvc_results_${permutation_number}.txt"
make size > "../../results/${permutation_number}_tsvc_size.txt"


# Parse and record the timing results
cd ../../
python3 parse_results.py "$permutation_number" "benchmarks/${GCC_DIR}/gcc_loops_${permutation_number}.csv" "benchmarks/${LCALS_DIR}/lcals_results_${permutation_number}" "benchmarks/${NPB_DIR}/npb_results_${permutation_number}.txt" "benchmarks/${TSVC_DIR}/tsvc_results_${permutation_number}.txt"

exit 0

