##
## File defining macro variables that specify (some of) LCALS compilation 
## options and behavior.
##
## This should be included in the makefile to build code that uses LCALS
## or replace with an equivalent substitute to customize behavior.  Many
## of these items can also be set manually in the LCALSParams.hxx file.
## 
## \author  Rich Hornung, Center for Applied Scientific Computing, LLNL
##

##
## Available options for LCALS scalar floating point types are:
##
## -DLCALS_USE_DOUBLE
## -DLCALS_USE_FLOAT
##
## Exactly one must be defined!!!
##
LCALS_FP_TYPE	= -DLCALS_USE_DOUBLE
#LCALS_FP_TYPE	= -DLCALS_USE_FLOAT


##
## Available options for LCALS floating point pointer types are:
##
## -DLCALS_USE_BARE_PTR
## -DLCALS_USE_RESTRICT_PTR
## -DLCALS_USE_RESTRICT_ALIGNED_PTR
## -DLCALS_USE_PTR_CLASS
##
## Exactly one must be defined!!!
##
#LCALS_FPPTR_TYPE	= -DLCALS_USE_BARE_PTR
LCALS_FPPTR_TYPE	= -DLCALS_USE_RESTRICT_PTR
#LCALS_FPPTR_TYPE	= -DLCALS_USE_RESTRICT_ALIGNED_PTR
#LCALS_FPPTR_TYPE	= -DLCALS_USE_PTR_CLASS


##
## Enable various LCALS checksum options:
##
## -DLCALS_VERIFY_CHECKSUM  - report checksums for all loop results
## -DLCALS_VERIFY_CHECKSUM_ABBREVIATED  - same as above, but only run a
##                                        few loop samples 
##
#LCALS_CHECKSUM_TYPE	= -DLCALS_VERIFY_CHECKSUM
LCALS_CHECKSUM_TYPE	= -DLCALS_VERIFY_CHECKSUM_ABBREVIATED


##
## Available options for LCALS loop execution timing mechanism:
##
## -DLCALS_USE_CYCLE  - use high-resolution cycle.h option
## -DLCALS_USE_CLOCK  - use std library time.h mechanism
##
## Exactly one must be defined!!!
##
#LCALS_TIMER_TYPE	= -DLCALS_USE_CYCLE
LCALS_TIMER_TYPE	= -DLCALS_USE_CLOCK

 
##
## Do OMP loop variants only.
##
## Must be defined or not defined, but not both!!!
##
#LCALS_DO_OMP_ONLY       = -DLCALS_DO_OMP_ONLY
LCALS_DO_OMP_ONLY       =



##
## The following defines -D flags on the compilation line...
##
LCALS_RULES	= $(LCALS_FP_TYPE) $(LCALS_FPPTR_TYPE) $(LCALS_FT_OPT) \
	$(LCALS_CHECKSUM_TYPE) $(LCALS_TIMER_TYPE) $(LCALS_DO_OMP_ONLY) \
	$(LCALS_LOOP_LENGTH_SETTINGS)
