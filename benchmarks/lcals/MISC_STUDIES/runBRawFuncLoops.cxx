//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "B" subset raw loops
//

#include "LCALSSuite.hxx"

#include "RawFuncLoops.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runBRawFuncLoops( std::vector<LoopStat>& loop_stats,
                       bool run_loop[],
                       LoopLength ilength )
{
   LoopSuiteRunInfo& loop_suite_run_info = getLoopSuiteRunInfo();
   LoopData& loop_data = getLoopData();

#if defined(COMPILE_RAW_VARIANTS)

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

            doINIT3(ltimer, num_samples, len,
                    out1, out2, out3, 
                    in1, in2);

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

            doMULADDSUB(ltimer, num_samples, len,
                        out1, out2, out3, 
                        in1, in2);

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

            doIF_QUAD(ltimer, num_samples, len,
                      a, b, c, 
                      x1, x2);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_TRAP_INT)
          case TRAP_INT : {


            stat.loop_is_run = false;

//
//  NOTE: There is no change to this loop for subroutine form.
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

#endif  // if COMPILE_RAW_VARIANTS

}
