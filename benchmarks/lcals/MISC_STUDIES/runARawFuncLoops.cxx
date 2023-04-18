//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "A" subset raw loops
//

#include "LCALSSuite.hxx"

#include "SubsetDataA.hxx"
#include "RawFuncLoops.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runARawFuncLoops( std::vector<LoopStat>& loop_stats,
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

#if defined(COMPILE_PRESSURE_CALC)
          case PRESSURE_CALC : {

            loopInit(iloop, stat);

            Real_ptr compression = loop_data.array_1D_Real[0];  
            Real_ptr bvc = loop_data.array_1D_Real[1];  
            Real_ptr p_new = loop_data.array_1D_Real[2];  
            Real_ptr e_old = loop_data.array_1D_Real[3];  
            Real_ptr vnewc = loop_data.array_1D_Real[4];  

            doPRESSURE_CALC(ltimer, num_samples, len,
                            compression, bvc, p_new, 
                            e_old, vnewc); 

            loopFinalize(iloop, stat, ilength);
     
            break;
          }
#endif

#if defined(COMPILE_PRESSURE_CALC_ALT)
          case PRESSURE_CALC_ALT : {

            stat.loop_is_run = false;

//
//  NOTE: Alternative version of loop kernel is same as above for this variant.
//
            break;
          }
#endif

#if defined(COMPILE_ENERGY_CALC)
          case ENERGY_CALC : {

            loopInit(iloop, stat);

            Real_ptr e_new = loop_data.array_1D_Real[0];
            Real_ptr e_old = loop_data.array_1D_Real[1];
            Real_ptr delvc = loop_data.array_1D_Real[2];
            Real_ptr p_new = loop_data.array_1D_Real[3];
            Real_ptr p_old = loop_data.array_1D_Real[4];
            Real_ptr q_new = loop_data.array_1D_Real[5];
            Real_ptr q_old = loop_data.array_1D_Real[6];
            Real_ptr work = loop_data.array_1D_Real[7];
            Real_ptr compHalfStep = loop_data.array_1D_Real[8];
            Real_ptr pHalfStep = loop_data.array_1D_Real[9];
            Real_ptr bvc = loop_data.array_1D_Real[10];
            Real_ptr pbvc = loop_data.array_1D_Real[11];
            Real_ptr ql_old = loop_data.array_1D_Real[12];
            Real_ptr qq_old = loop_data.array_1D_Real[13];
            Real_ptr vnewc = loop_data.array_1D_Real[14];

            doENERGY_CALC(ltimer, num_samples, len,
                          e_new, e_old, delvc,
                          p_new, p_old, q_new,
                          q_old, work, compHalfStep,
                          pHalfStep, bvc, pbvc,
                          ql_old, qq_old, vnewc);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_ENERGY_CALC_ALT)
          case ENERGY_CALC_ALT : {

            stat.loop_is_run = false;

//
//  NOTE: Alternative version of loop kernel is same as above for this variant.
//
            break;
          }
#endif

#if defined(COMPILE_VOL3D_CALC)
          case VOL3D_CALC : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr z = loop_data.array_1D_Real[2];
            Real_ptr vol = loop_data.array_1D_Real[3];

            ADomain domain(ilength, /* ndims = */ 3);

            doVOL3D_CALC(ltimer, num_samples, domain,
                         x, y, z, vol);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_DEL_DOT_VEC_2D)
          case DEL_DOT_VEC_2D : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr xdot = loop_data.array_1D_Real[2];
            Real_ptr ydot = loop_data.array_1D_Real[3];
            Real_ptr div = loop_data.array_1D_Real[4];

            ADomain domain(ilength, /* ndims = */ 2);

            doDEL_DOT_VEC_2D(ltimer, num_samples, domain,
                             x, y, xdot, ydot, div);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_COUPLE)
          case COUPLE : {

            loopInit(iloop, stat);

            Complex_ptr t0 = loop_data.array_1D_Complex[0];
            Complex_ptr t1 = loop_data.array_1D_Complex[1];
            Complex_ptr t2 = loop_data.array_1D_Complex[2];
            Complex_ptr denac = loop_data.array_1D_Complex[3];
            Complex_ptr denlw = loop_data.array_1D_Complex[4];

            ADomain domain(ilength, /* ndims = */ 3);

            doCOUPLE(ltimer, num_samples, domain,
                     t0, t1, t2, denac, denlw);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_FIR)
          case FIR : {

            loopInit(iloop, stat);

            Real_ptr out = loop_data.array_1D_Real[0];
            Real_ptr in = loop_data.array_1D_Real[1];

            Index_type val = 
            doFIR(ltimer, num_samples, len, 
                  out, in);

            //
            // RDH added this. Without it compiler may optimize out outer 
            // sampling loop since each sample pass results in identical output.
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

#endif  // if COMPILE_RAW_VARIANTS

}
