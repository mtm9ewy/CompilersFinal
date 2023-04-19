SHELL=/bin/sh
BENCHMARK=bt
BENCHMARKU=BT

include ../config/make.def

OBJS =	bt.o \
	${COMMON}/c_print_results.o ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

bt.o:             bt.cpp npbparams.hpp
	${CC} -S -c -emit-llvm -o bt.bc bt.cpp
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o bt.ll bt.bc
	$(CC) -c bt.ll

clean:
	- rm -f *.o *~ mputil*
	- rm -f  npbparams.hpp core
