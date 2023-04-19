SHELL=/bin/sh
CLASS=S
SFILE=config/suite.def

default: header
	@ $(SHELL) sys/print_instructions

	       
BT: bt
bt: header
	cd BT; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)
		       
SP: sp		       
sp: header	       
	cd SP; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)
		       
LU: lu		       
lu: header	       
	cd LU; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)
		       
MG: mg		       
mg: header	       
	cd MG; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)
		       
FT: ft		       
ft: header	       
	cd FT; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)
		       
IS: is		       
is: header	       
	cd IS; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)
		       
CG: cg		       
cg: header	       
	cd CG; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)
		       
EP: ep		       
ep: header	       
	cd EP; $(MAKE) CLASS=$(CLASS) OPTFLAGS=$(OPTFLAGS)


# Awk script courtesy cmg@cray.com
suite:
	@ awk '{ if ($$1 !~ /^#/ &&  NF > 0)                              \
	printf "make %s CLASS=%s\n", $$1, $$2 }' $(SFILE)  \
	| $(SHELL)

run: 
	./bin/bt.W
	./bin/cg.W
	./bin/ep.W
	./bin/ft.W
	./bin/is.W
	./bin/lu.W
	./bin/mg.W
	./bin/sp.W

size: 
	du bin/*.W


# It would be nice to make clean in each subdirectory (the targets
# are defined) but on a really clean system this will won't work
# because those makefiles need config/make.def
clean:
	- rm -f core 
	- rm -f *~ */core */*~ */*.o */npbparams.hpp */*.obj */*.exe
	- rm -f sys/setparams sys/makesuite sys/setparams.hpp

cleanall: clean
	- rm -r bin/*

veryclean: clean
	- rm config/make.def config/suite.def Part*
	- rm bin/sp.* bin/lu.* bin/mg.* bin/ft.* bin/bt.* bin/is.* bin/ep.* bin/cg.*

header:
	@ $(SHELL) sys/print_header

kit: 
	- makekit -s100k -k30 * */* */*/*

