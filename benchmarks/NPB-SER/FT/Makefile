SHELL=/bin/sh
BENCHMARK=ft
BENCHMARKU=FT

include ../config/make.def

OBJS = ft.o ${COMMON}/c_${RAND}.o ${COMMON}/c_print_results.o \
       ${COMMON}/c_timers.o ${COMMON}/c_wtime.o #../omp-prof.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

ft.o:             ft.cpp npbparams.hpp
	${CC} -S -c -emit-llvm -o ft.bc ft.cpp
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o ft.ll ft.bc
	$(CC) -c ft.ll

clean:
	- rm -f *.o *~ mputil*
	- rm -f ft npbparams.hpp core
