SHELL=/bin/sh
BENCHMARK=mg
BENCHMARKU=MG

include ../config/make.def

OBJS = mg.o ${COMMON}/c_print_results.o  \
       ${COMMON}/c_${RAND}.o ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

mg.o:		mg.cpp npbparams.hpp 
	${CC} -S -c -emit-llvm -o mg.bc mg.cpp
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o mg.ll mg.bc
	$(CC) -c mg.ll

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.hpp core
