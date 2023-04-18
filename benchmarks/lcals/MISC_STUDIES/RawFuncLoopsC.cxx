//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "C" subset raw loops implemented in
// individual subroutines.
//

#include "LCALSSuite.hxx"

#include "RawFuncLoops.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>

#if defined(COMPILE_HYDRO_1D)
void doHYDRO_1D(LoopTimer& ltimer, int num_samples, Index_type len,
                Real_ptr x, Real_ptr y, Real_ptr z)
{
   LoopData& loop_data = getLoopData();

   const Real_type q = loop_data.scalar_Real[0];
   const Real_type r = loop_data.scalar_Real[1];
   const Real_type t = loop_data.scalar_Real[2];

   TIMER_START(ltimer);
   for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

      for (Index_type k=0 ; k<len ; k++ ) {
        x[k] = q + y[k]*( r*z[k+10] + t*z[k+11] );
      }

   }
   TIMER_STOP(ltimer);
}
#endif
