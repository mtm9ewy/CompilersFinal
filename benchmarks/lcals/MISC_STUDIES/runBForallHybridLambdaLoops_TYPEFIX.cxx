//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "B" subset forall hybrid lambda loops
//    with type information "fix"
//

#include "LCALSSuite.hxx"
#include "LCALSTraversalMethods.hxx"
#include "SubsetDataB.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runBForallHybridLambdaLoops_TYPEFIX( std::vector<LoopStat>& loop_stats,
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his,
               [&] (Index_type i) {
                  Real_ptr tout1 = out1;
                  Real_ptr tout2 = out2;
                  Real_ptr tout3 = out3;
                  tout1[i] = tout2[i] = tout3[i] = - in1[i] - in2[i];
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type i) {
                  Real_ptr tout1 = out1;
                  Real_ptr tout2 = out2;
                  Real_ptr tout3 = out3;
                  tout1[i] = in1[i] * in2[i] ;
                  tout2[i] = in1[i] + in2[i] ;
                  tout3[i] = in1[i] - in2[i] ;
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<exec_policy>(his,
               [&] (Index_type i) {
                 Real_ptr ta = a;
                  Real_ptr tb = b;
                  Real_ptr tc = c;
                  Real_ptr tx1 = x1;
                  Real_ptr tx2 = x2;

                  Real_type s = tb[i]*tb[i] - 4.0*ta[i]*tc[i];
                  if ( s >= 0 ) {
                     s = sqrt(s);
                     tx2[i] = (-tb[i]+s)/(2.0*ta[i]);
                     tx1[i] = (-tb[i]-s)/(2.0*ta[i]);
                  } else {
                     tx2[i] = 0.0;
                     tx1[i] = 0.0;
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

            stat.loop_is_run = false;

//
//  NOTE: There is no change to this loop for type fix.
//
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

