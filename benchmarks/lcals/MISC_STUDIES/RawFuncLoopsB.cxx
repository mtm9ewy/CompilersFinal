//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "B" subset raw loops implemented in
// individual subroutines.
//

#include "LCALSSuite.hxx"

#include "RawFuncLoops.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>

#if defined(COMPILE_INIT3)
void doINIT3(LoopTimer& ltimer, int num_samples, Index_type len,
             Real_ptr out1, Real_ptr out2, Real_ptr out3, 
             Real_ptr in1, Real_ptr in2)
{
   TIMER_START(ltimer);
   for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

      for (Index_type i=0 ; i<len ; i++ ) {
        out1[i] = out2[i] = out3[i] = - in1[i] - in2[i];
      }

   }
   TIMER_STOP(ltimer);
}
#endif


#if defined(COMPILE_MULADDSUB)
void doMULADDSUB(LoopTimer& ltimer, int num_samples, Index_type len,
                 Real_ptr out1, Real_ptr out2, Real_ptr out3, 
                 Real_ptr in1, Real_ptr in2)
                 
{
   TIMER_START(ltimer);
   for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

      for (Index_type i=0 ; i<len ; i++ ) {
#if defined(XLC_ALIGN_DECORATE)
        __alignx(LCALS_DATA_ALIGN, out1);
        __alignx(LCALS_DATA_ALIGN, out2);
        __alignx(LCALS_DATA_ALIGN, out3);
        __alignx(LCALS_DATA_ALIGN, in1);
        __alignx(LCALS_DATA_ALIGN, in2);
#endif
        out1[i] = in1[i] * in2[i] ;
        out2[i] = in1[i] + in2[i] ;
        out3[i] = in1[i] - in2[i] ;
      }

   }
   TIMER_STOP(ltimer);
}
#endif


#if defined(COMPILE_IF_QUAD)
void doIF_QUAD(LoopTimer& ltimer, int num_samples, Index_type len,
               Real_ptr a, Real_ptr b, Real_ptr c,
               Real_ptr x1, Real_ptr x2)
               
{
   TIMER_START(ltimer);
   for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

      for (Index_type i=0 ; i<len ; i++ ) {
         Real_type s = b[i]*b[i] - 4.0*a[i]*c[i];
         if ( s >= 0 ) {
            s = sqrt(s);
            x2[i] = (-b[i]+s)/(2.0*a[i]);
            x1[i] = (-b[i]-s)/(2.0*a[i]);
         } else {
            x2[i] = 0.0;
            x1[i] = 0.0;
         }
      }

   }
   TIMER_STOP(ltimer);
}
#endif
