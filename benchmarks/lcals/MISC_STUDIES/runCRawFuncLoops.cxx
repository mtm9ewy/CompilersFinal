//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "C" subset raw loops
//

#include "LCALSSuite.hxx"

#include "RawFuncLoops.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runCRawFuncLoops( std::vector<LoopStat>& loop_stats,
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

#if defined(COMPILE_HYDRO_1D)
          case HYDRO_1D : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr z = loop_data.array_1D_Real[2];

            doHYDRO_1D(ltimer, num_samples, len,
                       x, y, z);

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

#endif  // if COMPILE_RAW_VARIANTS

}
