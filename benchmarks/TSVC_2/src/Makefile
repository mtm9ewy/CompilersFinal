.PHONY: all clean

libs=-lm
CC=clang -O -Xclang -disable-llvm-passes
OPT=opt

all: tsvc

tsvc: tsvc.o dummy.o common.o
	$(CC) $^ $(libs) -o tsvc

%.o : %.c 
	$(CC) -S -c -emit-llvm -o $*.bc $<
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o $*.ll $*.bc
	$(CC) -c -o $@ $*.ll

clean :
	rm -f *.o *.s *.ll *.bc

run :
	./tsvc

size :
	du tsvc
