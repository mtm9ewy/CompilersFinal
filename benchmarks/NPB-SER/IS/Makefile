SHELL=/bin/sh
BENCHMARK=is
BENCHMARKU=IS

include ../config/make.def

OBJS = is.o ${COMMON}/c_print_results.o ${COMMON}/c_${RAND}.o \
       ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}


is.o:		is.cpp npbparams.hpp
	${CC} -S -c -emit-llvm -o is.bc is.cpp
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o is.ll is.bc
	$(CC) -c is.ll

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.hpp core
