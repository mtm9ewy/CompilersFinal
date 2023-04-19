# CompilersFinal

### Running
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
