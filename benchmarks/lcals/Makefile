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
#OPT_DEBUG = debug

#rzalastor
#LCALS_ARCH = x86_sse_icc
#LCALS_ARCH = x86_sse_gnu

#rzmerl
#LCALS_ARCH = x86_avx_icc
#LCALS_ARCH = x86_avx_gnu

#rzdawndev
## This is still here for posterity...
#LCALS_ARCH = bgp_xlc9

#rzuseq
#LCALS_ARCH = bgq_xlc12
LCALS_ARCH = bgq_clang
#LCALS_ARCH = bgq_gnu

#rzmic
#LCALS_ARCH = MIC


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
#  MIC
# 
ifeq ($(LCALS_ARCH),MIC)
CXX 		= icc

ifeq ($(OPT_DEBUG),opt)
CXX_COMPILE = $(CXX) -g -O3 -mmic -vec-report3  -inline-max-total-size=10000 -inline-forceinline -ansi-alias -std=c++0x -openmp
LDFLAGS	=
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -mmic  -O0 -std=c++0x -openmp
LDFLAGS	=
endif

CXXFLAGS	= -DLCALS_PLATFORM_X86_SSE -DLCALS_COMPILER_ICC
LDPATH		=

endif 


#
#  INTEL COMPILER (SSE 4.1) 
# 
ifeq ($(LCALS_ARCH),x86_sse_icc)
#CXX 		= /usr/local/tools/icc-12.1.339/bin/icc
#CXX 		= /usr/local/tools/ic-13.0.117/bin/icpc
#CXX 		= /usr/local/tools/ic-13.0.146/bin/icc
#CXX 		= /usr/local/tools/ic-13.1.163/bin/icpc
#CXX 		= /usr/local/tools/ic-14.0.080/bin/icc
#CXX 		= /usr/local/tools/ic-14.0.097/bin/icpc
#CXX 		= /usr/local/tools/ic-14.0.106/bin/icpc

ifeq ($(OPT_DEBUG),opt)
#CXX_COMPILE = $(CXX) -O3 -msse4.1 -inline-max-total-size=10000 -inline-forceinline -opt-streaming-stores always -ansi-alias -std=c++0x -openmp
CXX_COMPILE = $(CXX) -O3 -msse4.1 -inline-max-total-size=10000 -inline-forceinline -ansi-alias -std=c++0x -openmp
#CXX_COMPILE = $(CXX) -O3 -msse4.1 -inline-max-total-size=10000 -inline-forceinline -ansi-alias -std=c++0x 
LDFLAGS	=
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -std=c++0x -openmp
LDFLAGS	=
endif

CXXFLAGS	= -DLCALS_PLATFORM_X86_SSE -DLCALS_COMPILER_ICC
LDPATH		=

endif 


#
#  GNU  COMPILER (SSE 4.1)
# 
ifeq ($(LCALS_ARCH),x86_sse_gnu)
#CXX 		= /usr/local/bin/g++-4.4.6
#CXX 		= /usr/local/bin/g++-4.6.1
#CXX 		= /usr/apps/gnu/4.7.1/bin/g++
#CXX 		= /usr/apps/gnu/4.8.0/bin/g++

ifeq ($(OPT_DEBUG),opt)
#
# Use this with GNU 4.7X and later
CXX_COMPILE = $(CXX) -Ofast -msse4.1 -finline-functions -finline-limit=10000 -std=c++11 -openmp
#CXX_COMPILE = $(CXX) -O3 -msse4.1 -finline-functions -finline-limit=10000 -std=c++0x -openmp
##CXX_COMPILE = $(CXX) -O3 -msse4.1  -ansi-alias -std=c++0x
## inline flags...
#CXX_COMPILE = $(CXX) -O3 -msse4.1  -finline-functions -finline-limit=10000 -ansi-alias -std=c++0x
## inline flags + others...
##CXX_COMPILE = $(CXX) -O3 -msse4.1  -finline-functions -finline-limit=10000 -fomit-frame-pointer -minline-all-stringops -malign-double -ftree-vectorize -floop-block -ansi-alias -std=c++0x -openmp
##CXX_COMPILE = $(CXX) -O3 -msse4.1  -finline-functions -finline-limit=10000 -fomit-frame-pointer -malign-double -ftree-vectorize -floop-block -ansi-alias -std=c++0x
LDFLAGS	=
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -std=c++0x -openmp
LDFLAGS	=
endif

CXXFLAGS 	= -DLCALS_PLATFORM_X86_SSE -DLCALS_COMPILER_GNU
LDPATH		=

endif 


#
#  INTEL COMPILER (AVX 2)
# 
ifeq ($(LCALS_ARCH),x86_avx_icc)
#CXX 		= /usr/local/tools/icc-12.1.339/bin/icpc
#CXX 		= /usr/local/tools/ic-13.0.117/bin/icpc
#CXX 		= /usr/local/tools/ic-13.0.146/bin/icpc
#CXX 		= /usr/local/tools/ic-13.1.163/bin/icpc
#CXX 		= /usr/local/tools/ic-14.0.080/bin/icpc
#CXX 		= /usr/local/tools/ic-14.0.097/bin/icpc
CXX 		= /usr/local/tools/ic-14.0.106/bin/icpc

ifeq ($(OPT_DEBUG),opt)
#CXX_COMPILE = $(CXX) -O3 -mavx -inline-max-total-size=10000 -inline-forceinline -opt-streaming-stores always -ansi-alias -std=c++0x -openmp
#CXX_COMPILE = $(CXX) -O3 -mavx -inline-max-total-size=10000 -inline-forceinline -ansi-alias -std=c++0x -openmp -static-intel
CXX_COMPILE = $(CXX) -O3 -mavx -inline-max-total-size=10000 -inline-forceinline -ansi-alias -std=c++0x -openmp
#CXX_COMPILE = $(CXX) -O3 -mavx -inline-max-total-size=10000 -inline-forceinline -ansi-alias -std=c++0x 
LDFLAGS	=
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -std=c++0x -openmp
LDFLAGS	=
endif

CXXFLAGS 	= -DLCALS_PLATFORM_X86_AVX -DLCALS_COMPILER_ICC
LDPATH		=

endif


#
#  GNU COMPILER (AVX 2)
# 
ifeq ($(LCALS_ARCH),x86_avx_gnu)
#CXX 		= /usr/local/bin/g++-4.6.1
#CXX 		= /usr/apps/gnu/4.7.1/bin/g++
#CXX 		= /usr/apps/gnu/4.8.0/bin/g++

ifeq ($(OPT_DEBUG),opt)
#
# Use this with GNU 4.7X and later
CXX_COMPILE = $(CXX) -Ofast -mavx -finline-functions -finline-limit=10000 -std=c++11 -openmp
#
# These should work with older compiler versions...
#CXX_COMPILE = $(CXX) -O3 -mavx -ansi-alias -std=c++0x -openmp
#CXX_COMPILE = $(CXX) -O3 -mavx -std=c++0x
##CXX_COMPILE = $(CXX) -O3 -mavx -ansi-alias -std=c++0x
##CXX_COMPILE = $(CXX) -O3 -mavx -finline-functions -finline-limit=10000 -fomit-frame-pointer -malign-double -ftree-vectorize -floop-block -ansi-alias -std=c++0x
LDFLAGS	=
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -std=c++0x -openmp
LDFLAGS	=
endif

CXXFLAGS 	= -DLCALS_PLATFORM_X86_AVX -DLCALS_COMPILER_GNU
LDPATH		=

endif


#
#  XLC 9 COMPILER  (BG/P)
# 
ifeq ($(LCALS_ARCH),bgp_xlc9)
##
## Other options to try "-O4 -qhot=level=0", "-qhot=vector -qhot=simd" 
##                                           (automatic woth -O4 or -O5) 
##
CXX 		= mpixlcxx

ifeq ($(OPT_DEBUG),opt)
CXX_COMPILE = $(CXX) -O3 -qarch=450d -qtune=450 -qalias=allp -qhot -qsmp=omp
LDFLAGS	= -O3 -Wl,--relax -qsmp=omp -qnostaticlink -export-dynamic -Wl,--relax
##CXX_COMPILE = $(CXX) -O3 -qarch=450d -qtune=450 -qalias=allp -qhot 
##LDFLAGS	= -O3 -Wl,--relax  -qnostaticlink -export-dynamic -Wl,--relax
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -qarch=450d -qtune=450 -qalias=allp -qsmp=omp
LDFLAGS	= -g -O0 -Wl,--relax -qsmp=omp -qnostaticlink -export-dynamic -Wl,--relax
endif

CXXFLAGS 	= -DLCALS_PLATFORM_BGP -DLCALS_COMPILER_XLC9 -I/usr/gapps/bdiv/sles_10_ppc64/opt/platform/include
LDPATH		= -L/usr/local/tools/lib

endif


#
#  XLC 12 COMPILER  (BG/Q)
# 
ifeq ($(LCALS_ARCH),bgq_xlc12)
#CXX 		= /usr/local/tools/compilers/ibm/bgxlc++-12.1.0.2a
#CXX 		= /usr/local/tools/compilers/ibm/bgxlc++-12.1.0.3
CXX 		= mpixlcxx_r
#CXX 		= /usr/local/tools/compilers/ibm/mpixlcxx-12.1.0.1d
#CXX 		= /usr/local/tools/compilers/ibm/mpixlcxx-12.1.0.2
#CXX 		= /usr/local/tools/compilers/ibm/mpixlcxx-12.1.0.3

ifeq ($(OPT_DEBUG),opt)
#CXX_COMPILE = $(CXX) -O3 -qarch=qp -qhot=novector -qsimd=auto -qlanglvl=extended0x -qstrict -qinline=10000 -qsmp=omp
CXX_COMPILE = $(CXX) -O3 -qarch=qp -qhot=novector -qsimd=auto -qlanglvl=extended0x -qnostrict -qinline=10000 -qsmp=omp
##
## USE THESE LINE TO GENERATE VECTORIZATION REPORT INFO...
#CXX_COMPILE = $(CXX) -O3 -qarch=qp -qhot=novector -qsimd=auto -qlanglvl=extended0x -qstrict -qinline=10000  -qlist -qsource -qlistopt -qreport
#CXX_COMPILE = $(CXX) -O3 -qarch=qp -qhot=novector -qsimd=auto -qlanglvl=extended0x -qnostrict -qinline=10000  -qlist -qsource -qlistopt -qreport
#CXX_COMPILE = $(CXX) -O3 -qarch=qp -qhot=novector -qsimd=auto -qlanglvl=extended0x -qinline=100 -qlistfmt=html=inlines
LDFLAGS	= 
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -qarch=qp -qlanglvl=extended0x -qsmp=omp
LDFLAGS	= 
endif

CXXFLAGS 	= -DLCALS_PLATFORM_BGQ -DLCALS_COMPILER_XLC12
LDPATH		=

endif


#
#  Clang C++ compiler for BG/Q
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
CXX             = clang++


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


#
#  GNU compiler for BG/Q
#
ifeq ($(LCALS_ARCH),bgq_gnu)
#CXX             = /bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/bin/g++
CXX             = /usr/local/tools/compilers/ibm/mpicxx-4.7.2

#Previous versions

ifeq ($(OPT_DEBUG),opt)
CXX_COMPILE = $(CXX) -O3 -finline-functions -finline-limit=10000 -std=c++0x -fopenmp

LDFLAGS =
endif

ifeq ($(OPT_DEBUG),debug)
CXX_COMPILE = $(CXX) -g -O0 -std=c++0x -openmp
LDFLAGS =
endif

CXXFLAGS        = -DLCALS_PLATFORM_BGQ -DLCALS_COMPILER_GNU
LDPATH          =

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
	$(CXX_COMPILE) $(OBJFILES) $(LDFLAGS) $(LDPATH) -o  $@.exe

## Uncomment this line to generate an executable to run...
%.o : %.cxx ; $(CXX_COMPILE) -c -o $@ $< $(CXX_FLAGS_BUILD)

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


