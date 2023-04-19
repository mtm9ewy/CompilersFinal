SHELL=/bin/sh
BENCHMARK=ep
BENCHMARKU=EP

include ../config/make.def

OBJS = ep.o ${COMMON}/c_print_results.o ${COMMON}/c_${RAND}.o \
       ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

ep.o:		ep.cpp npbparams.hpp
	${CC} -S -c -emit-llvm -o ep.bc ep.cpp
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o ep.ll ep.bc
	$(CC) -c ep.ll

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.hpp core
