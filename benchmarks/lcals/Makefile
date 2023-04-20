#  
# See README-LCALS_license.txt for access and distribution restrictions
#

#
#  Modify stuff in this Makfile to fit your machine and compiler.
#
#  Uncomment ONE AND ONLY ONE opt/debug option and "LCALS_ARCH" 
#  variable line for the platform you want. 
#
# IMPORTANT:  Make sure CXXFLAGS are what you want. They are used in 
#             the source code, LCALSParams.hxx in particular, to
#             set code compilation options. 
# 
#             Also, make sure the options set in the LCALS_rules.mk file are 
#             what you want.  
#
#             Finally, make sure the VERSION of the compiler is what you want!
#

OPT_DEBUG = opt
LCALS_ARCH = bgq_clang

##
## The LCALS_rules.mk file defines macro variables that specify LCALS behavior.
## To change the rules, the file in the LCALS include directory can be edited
## or it can be replaced with a custom version here.
##
include ./LCALS_rules.mk


#
#  By default, only the code in the source files associated with the
#  LCALS benchmark will be compiled.  To run other loop variants,
#  define BUILD_MISC compile loop variants in MISC_STUDIES directory.
#
#  This is necessary to build an executable if other variants are 
#  specified to run in main.cxx.
#
#BUILD_MISC = yes


#
#  Clang C++ compiler
# 
ifeq ($(LCALS_ARCH),bgq_clang)
#CXX 		= /usr/apps/gnu/clang/bin/mpiclang++11
#CXX 		= /usr/apps/gnu/clang/bin/bgclang++11

#Specific versions
#CXX            = /usr/apps/gnu/clang/r176829-20130309/bin/bgclang++11
#CXX             = /usr/apps/gnu/clang/r176751-20130307/bin/bgclang++11
#CXX             = /usr/apps/gnu/clang/r181589-20130510/bin/bgclang++11
#CXX             = /usr/apps/gnu/clang/r189357-20130827/bin/bgclang++11
#CXX             = /usr/apps/gnu/clang/r190771-20130914/bin/bgclang++11
#CXX             = /usr/apps/gnu/clang/r192411-20131010/bin/bgclang++11
CXX             = clang++-14 --std=c++14
OPT             = opt-14



ifeq ($(OPT_DEBUG),opt)
#CXX_COMPILE = $(CXX) -O3 -finline-functions -finline-limit=10000 -fomit-frame-pointer -minline-all-stringops -malign-double -ftree-vectorize -floop-block -ansi-alias -std=c++0x
#Opt 3
#CXX_COMPILE = $(CXX) -O3 -finline-functions -finline-limit=10000 -malign-double -std=c++0x
#Opt 2
CXX_COMPILE = $(CXX) -O3 -finline-functions  -ffast-math -std=c++0x
#Opt 1
#CXX_COMPILE = $(CXX) -O3 -finline-functions  -std=c++0x
#Opt 0
#CXX_COMPILE = $(CXX) -O0 -finline-functions  -std=c++0x
LDFLAGS	=
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -std=c++0x
LDFLAGS	=
endif

CXXFLAGS 	= -DLCALS_PLATFORM_BGQ -DLCALS_COMPILER_CLANG
LDPATH		=

endif 


OBJS := $(patsubst %.cxx,%.o,$(wildcard ./*.cxx))

ifdef BUILD_MISC
MISC_OBJS := $(patsubst %.cxx,%.o,$(wildcard ./MISC_STUDIES/*.cxx))
MISCCXXFLAG = -DLCALS_DO_MISC 
else
MISC_OBJS = 
MISCCXXFLAG = 
endif

OBJFILES = $(OBJS) $(MISC_OBJS)

CXX_FLAGS_BUILD = -I. $(CXXFLAGS) $(LCALS_RULES) $(MISCCXXFLAG)

lcals: $(OBJFILES)
	@echo "Creating ./lcalsversioninfo.txt file" && ./lcals_version \
	"$(CXX)"                                                  \
	"CXX_COMPILE      info   $(CXX_COMPILE)"                  \
	"CXX_FLAGS_BUILD  info   $(CXX_FLAGS_BUILD)"              \
	"LDFLAGS          info   $(LDFLAGS)"
	$(CXX_COMPILE) $(OBJFILES) $(LDFLAGS) $(LDPATH) -o  $@

# Uncomment this line to generate an executable to run...
%.o : %.cxx 
	# $(CXX_COMPILE) -c -o $@ $< $(CXX_FLAGS_BUILD)
	$(CXX) -S -c -emit-llvm -o $*.bc $< $(CXX_FLAGS_BUILD)
	$(OPT) -enable-new-pm=0 $(OPTFLAGS) -S -o $*.ll $*.bc
	$(CXX) -c -o $@ $*.ll

run: 
	./lcals $(OUTDIR)

size: 
	du lcals

## The following can be used to generate vectorization report and 
## assmbly code output for Intel compilers...
#%.o : %.cxx ; $(CXX_COMPILE) -S -o $@ $< $(CXX_FLAGS_BUILD)
#%.o : %.cxx ; $(CXX_COMPILE) -S -vec-report3 -o $@ $< $(CXX_FLAGS_BUILD)
#%.o : %.cxx ; $(CXX_COMPILE) -S -fcode-asm -vec-report3 -o $@ $< $(CXX_FLAGS_BUILD)

## The following can be used to generate vectorization report and 
## assmbly code output for IBM XLC compilers...
## See compiler flag options in the appropriate section above.
#%.o : %.cxx ; $(CXX_COMPILE) -c -o $@ $< $(CXX_FLAGS_BUILD)

clean-obj: 
	rm -rf *.o MISC_STUDIES/*.o

clean: clean-obj
	rm -rf *.s *.lst *.exe lcalsversioninfo.txt*


