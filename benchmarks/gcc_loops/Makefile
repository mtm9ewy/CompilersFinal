.PHONY: all clean

CC=clang++-14
OPT=opt-14

all: gcc-loops

gcc-loops: FORCE
	$(CC) -Wno-align-mismatch -S -c -emit-llvm -o gcc-loops.bc gcc-loops.cpp
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o gcc-loops.ll gcc-loops.bc
	$(CC) -c -o gcc-loops.o gcc-loops.ll
	$(CC) gcc-loops.o -o gcc-loops

clean :
	rm -f *.o *.s *.ll *.bc

run:
	./gcc-loops

size:
	du gcc-loops



FORCE: ;

