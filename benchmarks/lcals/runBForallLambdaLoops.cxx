//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "B" subset forall lambda loops
//

#include "LCALSSuite.hxx"
#include "LCALSTraversalMethods.hxx"
#include "SubsetDataB.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runBForallLambdaLoops( std::vector<LoopStat>& loop_stats,
                            bool run_loop[],
                            LoopLength ilength )
{
   LoopSuiteRunInfo& loop_suite_run_info = getLoopSuiteRunInfo();
   LoopData& loop_data = getLoopData();

#if defined(COMPILE_LAMBDA_VARIANTS)

   for (unsigned iloop = 0; iloop < loop_suite_run_info.num_loops; ++iloop) {

      if ( run_loop[iloop] ) {

         LoopStat& stat = loop_stats[iloop];
         Index_type len = stat.loop_length[ilength];
         int num_samples = stat.samples_per_pass[ilength];
#if defined(LCALS_VERIFY_CHECKSUM_ABBREVIATED)
         num_samples = num_checksum_samples;
#endif

         LoopTimer ltimer;

         switch ( iloop ) {

#if defined(COMPILE_INIT3)
          case INIT3 : {

            loopInit(iloop, stat);

            Real_ptr out1 = loop_data.array_1D_Real[0];
            Real_ptr out2 = loop_data.array_1D_Real[1];
            Real_ptr out3 = loop_data.array_1D_Real[2];
            Real_ptr in1 = loop_data.array_1D_Real[3];
            Real_ptr in2 = loop_data.array_1D_Real[4];

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len,
               [&] (Index_type i) {
                  out1[i] = out2[i] = out3[i] = - in1[i] - in2[i];
               } );
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_MULADDSUB)
          case MULADDSUB : {

            loopInit(iloop, stat);

            Real_ptr out1 = loop_data.array_1D_Real[0];
            Real_ptr out2 = loop_data.array_1D_Real[1];
            Real_ptr out3 = loop_data.array_1D_Real[2];
            Real_ptr in1 = loop_data.array_1D_Real[3];
            Real_ptr in2 = loop_data.array_1D_Real[4]; 

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, 
               [&] (Index_type i) {
                  out1[i] = in1[i] * in2[i] ;
                  out2[i] = in1[i] + in2[i] ;
                  out3[i] = in1[i] - in2[i] ;
               } );
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_IF_QUAD)
          case IF_QUAD : {

            loopInit(iloop, stat);

            Real_ptr a = loop_data.array_1D_Real[0];
            Real_ptr b = loop_data.array_1D_Real[1];
            Real_ptr c = loop_data.array_1D_Real[2];
            Real_ptr x1 = loop_data.array_1D_Real[3];
            Real_ptr x2 = loop_data.array_1D_Real[4];

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<exec_policy>(0, len, 
               [&] (Index_type i) {
                  Real_type s = b[i]*b[i] - 4.0*a[i]*c[i];
                  if ( s >= 0 ) {
                     s = sqrt(s); 
                     x2[i] = (-b[i]+s)/(2.0*a[i]);
                     x1[i] = (-b[i]-s)/(2.0*a[i]);
                  } else {
                     x2[i] = 0.0;
                     x1[i] = 0.0;
                  }
               } );

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_TRAP_INT)
          case TRAP_INT : {

            loopInit(iloop, stat);

            Real_type xn = loop_data.scalar_Real[0];
            Real_type x0 = loop_data.scalar_Real[1];
            Real_type xp = loop_data.scalar_Real[2];
            Real_type y = loop_data.scalar_Real[3];
            Real_type yp = loop_data.scalar_Real[4];

            Index_type nx = loop_data.array_1D_Indx[0][0] + 1;

            const Real_type h = (xn - x0) / nx;
            Real_type sumx = 0.5*( trap_int_func(x0, y, xp, yp) +
                                   trap_int_func(xn, y, xp, yp) );

            Real_type val = 0;

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<exec_policy>(0, len, 
               [&] (Index_type i) {
                  Real_type x = x0 + i*h;
                  sumx += trap_int_func(x, y, xp, yp);
               } );

               val = sumx * h;

            }
            TIMER_STOP(ltimer);

            //
            // RDH put this here to prevent compiler from optimizing anything
            // out since result of loop (sumx) would not be used otherwise.
            //
            loop_data.scalar_Real[0] =
               (val + 0.00123) / (val - 0.00123);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

          default: {
//          std::cout << "\n Unknown loop id = " << iloop << std::endl;
          }

         } // switch on loop id

         copyTimer(stat, ilength, ltimer);

      }  // if loop with id should be run 

   }  // for loop over loops

#endif  // if COMPILE_LAMBDA_VARIANTS

}
