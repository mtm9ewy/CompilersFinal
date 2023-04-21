#!/bin/bash

#SBATCH --job-name="SLURM Loop Optimization Ordering"
#SBATCH --error="loop_pass.err"
#SBATCH --output="loop_pass.output"

source /etc/profile.d/modules.sh
module load python3
module load clang-llvm-14.0.6

start="$1"
stop="$2"

mkdir -p results

GCC_DIR="gcc_loops_${start}"
LCALS_DIR="lcals_${start}"
NPB_DIR="NPB-SER_${start}"
TSVC_DIR="TSVC_2_${start}"

# Create working directories
cp -r benchmarks/gcc_loops "benchmarks/${GCC_DIR}"
cp -r benchmarks/lcals "benchmarks/${LCALS_DIR}"
cp -r benchmarks/NPB-SER "benchmarks/${NPB_DIR}"
cp -r benchmarks/TSVC_2/src "benchmarks/${TSVC_DIR}"

for i in $(seq "${start}" "${stop}")
do
    echo "$i"
    bash run_test.sh "$i" "$GCC_DIR" "$LCALS_DIR" "$NPB_DIR" "$TSVC_DIR"
done

# Remove working directories
rm -r "benchmarks/${GCC_DIR}"
rm -r "benchmarks/${LCALS_DIR}"
rm -r "benchmarks/${NPB_DIR}"
rm -r "benchmarks/${TSVC_DIR}"
