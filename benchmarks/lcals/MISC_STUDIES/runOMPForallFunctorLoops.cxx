//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS OpenMP forall functor loops
//
// NOTE: These loops are a subset of other loops.  Which subset is indicated
//       at top of each loop case.
//

#include "LCALSSuite.hxx"
#include "LCALSTraversalMethods.hxx"

#include "SubsetDataA.hxx"
#include "SubsetDataB.hxx"
#include "FunctorKernelsA.hxx"
#include "FunctorKernelsB.hxx"
#include "FunctorKernelsC.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runOMPForallFunctorLoops( std::vector<LoopStat>& loop_stats,
                               bool run_loop[],
                               LoopLength ilength )
{
   LoopSuiteRunInfo& loop_suite_run_info = getLoopSuiteRunInfo();
   LoopData& loop_data = getLoopData();

#if defined(COMPILE_FUNCTOR_VARIANTS) && defined(COMPILE_OMP_VARIANTS) 

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
          case PRESSURE_CALC : {  // --> from subset "A"

            loopInit(iloop, stat);

            Real_ptr compression = loop_data.array_1D_Real[0];
            Real_ptr bvc = loop_data.array_1D_Real[1];
            Real_ptr p_new = loop_data.array_1D_Real[2];
            Real_ptr e_old = loop_data.array_1D_Real[3];
            Real_ptr vnewc = loop_data.array_1D_Real[4];

            const Real_type cls = loop_data.scalar_Real[0];
            const Real_type p_cut = loop_data.scalar_Real[1];
            const Real_type pmin = loop_data.scalar_Real[2];
            const Real_type eosvmax = loop_data.scalar_Real[3];

            PRESSURE_CALC_FunctorA kernela(bvc, compression, cls);

            PRESSURE_CALC_FunctorB kernelb(p_new, bvc, e_old, vnewc,
                                           p_cut, eosvmax, pmin);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

             #pragma omp parallel
              {
               forall<omp_for_nowait_exec>(0, len, kernela);

               forall<omp_for_nowait_exec>(0, len, kernelb);
              }  // omp parallel

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_PRESSURE_CALC_ALT)
          case PRESSURE_CALC_ALT : {  // --> from subset "A"

//
// NOTE: This is a slight variation on the OpenMP execution pattern
//       in the kernel above.
//

            loopInit(iloop, stat);

            Real_ptr compression = loop_data.array_1D_Real[0];
            Real_ptr bvc = loop_data.array_1D_Real[1];
            Real_ptr p_new = loop_data.array_1D_Real[2];
            Real_ptr e_old = loop_data.array_1D_Real[3];
            Real_ptr vnewc = loop_data.array_1D_Real[4];

            const Real_type cls = loop_data.scalar_Real[0];
            const Real_type p_cut = loop_data.scalar_Real[1];
            const Real_type pmin = loop_data.scalar_Real[2];
            const Real_type eosvmax = loop_data.scalar_Real[3];

            PRESSURE_CALC_FunctorA kernela(bvc, compression, cls);

            PRESSURE_CALC_FunctorB kernelb(p_new, bvc, e_old, vnewc,
                                           p_cut, eosvmax, pmin);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<omp_parallel_for_exec>(0, len, kernela);

               forall<omp_parallel_for_exec>(0, len, kernelb);

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_ENERGY_CALC)
          case ENERGY_CALC : {  // --> from subset "A"

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

            const Real_type rho0 = loop_data.scalar_Real[0];
            const Real_type e_cut = loop_data.scalar_Real[1];
            const Real_type emin = loop_data.scalar_Real[2];
            const Real_type q_cut = loop_data.scalar_Real[3];

            ENERGY_CALC_FunctorA kernela(e_new,
                                         e_old, delvc, p_old, q_old, work);

            ENERGY_CALC_FunctorB kernelb(q_new,
                                         delvc, compHalfStep, pHalfStep,
                                         e_new, bvc, pbvc,
                                         ql_old, qq_old,
                                         rho0);

            ENERGY_CALC_FunctorC kernelc(e_new,
                                         delvc, pHalfStep,
                                         p_old, q_old, q_new); 

            ENERGY_CALC_FunctorD kerneld(e_new,
                                         work, 
                                         e_cut, emin);

            ENERGY_CALC_FunctorE kernele(q_new,
                                         delvc, 
                                         e_new, p_new, vnewc,
                                         bvc, pbvc,
                                         ql_old, qq_old,
                                         p_old, q_old, pHalfStep,
                                         rho0, e_cut, emin);

            ENERGY_CALC_FunctorF kernelf(q_new,
                                         delvc, 
                                         e_new, p_new, vnewc,
                                         bvc, pbvc,
                                         ql_old, qq_old,
                                         rho0, q_cut);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

              #pragma omp parallel
              {
               forall<omp_for_nowait_exec>(0, len, kernela);

               forall<omp_for_nowait_exec>(0, len, kernelb);

               forall<omp_for_nowait_exec>(0, len, kernelc);

               forall<omp_for_nowait_exec>(0, len, kerneld);

               forall<omp_for_nowait_exec>(0, len, kernele);

               forall<omp_for_nowait_exec>(0, len, kernelf);
              } // omp parallel

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_ENERGY_CALC_ALT)
          case ENERGY_CALC_ALT : {  // --> from subset "A"

//
// NOTE: This is a slight variation on the OpenMP execution pattern
//       in the kernel above.
//

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

            const Real_type rho0 = loop_data.scalar_Real[0];
            const Real_type e_cut = loop_data.scalar_Real[1];
            const Real_type emin = loop_data.scalar_Real[2];
            const Real_type q_cut = loop_data.scalar_Real[3];

            ENERGY_CALC_FunctorA kernela(e_new,
                                         e_old, delvc, p_old, q_old, work);

            ENERGY_CALC_FunctorB kernelb(q_new,
                                         delvc, compHalfStep, pHalfStep,
                                         e_new, bvc, pbvc,
                                         ql_old, qq_old,
                                         rho0);

            ENERGY_CALC_FunctorC kernelc(e_new,
                                         delvc, pHalfStep,
                                         p_old, q_old, q_new); 

            ENERGY_CALC_FunctorD kerneld(e_new,
                                         work, 
                                         e_cut, emin);

            ENERGY_CALC_FunctorE kernele(q_new,
                                         delvc, 
                                         e_new, p_new, vnewc,
                                         bvc, pbvc,
                                         ql_old, qq_old,
                                         p_old, q_old, pHalfStep,
                                         rho0, e_cut, emin);

            ENERGY_CALC_FunctorF kernelf(q_new,
                                         delvc, 
                                         e_new, p_new, vnewc,
                                         bvc, pbvc,
                                         ql_old, qq_old,
                                         rho0, q_cut);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<omp_parallel_for_exec>(0, len, kernela);

               forall<omp_parallel_for_exec>(0, len, kernelb);

               forall<omp_parallel_for_exec>(0, len, kernelc);

               forall<omp_parallel_for_exec>(0, len, kerneld);

               forall<omp_parallel_for_exec>(0, len, kernele);

               forall<omp_parallel_for_exec>(0, len, kernelf);

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_VOL3D_CALC)
          case VOL3D_CALC : {   // --> from subset "A"

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr z = loop_data.array_1D_Real[2];
            Real_ptr vol = loop_data.array_1D_Real[3];

            ADomain domain(ilength, /* ndims = */ 3);

            UnalignedReal_ptr x0,x1,x2,x3,x4,x5,x6,x7 ;
            UnalignedReal_ptr y0,y1,y2,y3,y4,y5,y6,y7 ;
            UnalignedReal_ptr z0,z1,z2,z3,z4,z5,z6,z7 ;

            NDPTRSET(x,x0,x1,x2,x3,x4,x5,x6,x7) ;
            NDPTRSET(y,y0,y1,y2,y3,y4,y5,y6,y7) ;
            NDPTRSET(z,z0,z1,z2,z3,z4,z5,z6,z7) ;

            const Real_type vnormq = 0.083333333333333333; /* vnormq = 1/12 */

            VOL3D_CALC_Functor kernel(vol,
                                      x0,x1,x2,x3,x4,x5,x6,x7,
                                      y0,y1,y2,y3,y4,y5,y6,y7,
                                      z0,z1,z2,z3,z4,z5,z6,z7,
                                      vnormq); 


            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<omp_parallel_for_exec>(domain.fpz, domain.lpz + 1, 
                                                  kernel);

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_DEL_DOT_VEC_2D)
          case DEL_DOT_VEC_2D : {   // --> from subset "A"

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr xdot = loop_data.array_1D_Real[2];
            Real_ptr ydot = loop_data.array_1D_Real[3];
            Real_ptr div = loop_data.array_1D_Real[4];

            ADomain domain(ilength, /* ndims = */ 2);

            UnalignedReal_ptr x1,x2,x3,x4 ;
            UnalignedReal_ptr y1,y2,y3,y4 ;
            UnalignedReal_ptr fx1,fx2,fx3,fx4 ;
            UnalignedReal_ptr fy1,fy2,fy3,fy4 ;

            NDSET2D(x,x1,x2,x3,x4) ;
            NDSET2D(y,y1,y2,y3,y4) ;
            NDSET2D(xdot,fx1,fx2,fx3,fx4) ;
            NDSET2D(ydot,fy1,fy2,fy3,fy4) ;

            const Real_type ptiny = 1.0e-20;

            DEL_DOT_VEC_2D_Functor kernel(domain.real_zones,
                                    div,
                                    x1,x2,x3,x4,
                                    y1,y2,y3,y4,
                                    fx1,fx2,fx3,fx4,
                                    fy1,fy2,fy3,fy4,
                                    ptiny);


            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<omp_parallel_for_exec>(0, domain.n_real_zones,
                                                  kernel);

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_COUPLE)
          case COUPLE : {  // --> from subset "A"

            loopInit(iloop, stat);

            Complex_ptr t0 = loop_data.array_1D_Complex[0];
            Complex_ptr t1 = loop_data.array_1D_Complex[1];
            Complex_ptr t2 = loop_data.array_1D_Complex[2];
            Complex_ptr denac = loop_data.array_1D_Complex[3];
            Complex_ptr denlw = loop_data.array_1D_Complex[4];


            ADomain domain(ilength, /* ndims = */ 3);

            Index_type imin = domain.imin;
            Index_type imax = domain.imax;
            Index_type jmin = domain.jmin;
            Index_type jmax = domain.jmax;
            Index_type kmin = domain.kmin;
            Index_type kmax = domain.kmax;

            const Real_type clight=3.e+10;
            const Real_type csound=3.09e+7;
            const Real_type omega0= 0.9;
            const Real_type omegar= 0.9;
            const Real_type dt= 0.208;
            const Real_type c10 = 0.25 * (clight / csound);
            const Real_type fratio = sqrt(omegar / omega0);
            const Real_type r_fratio = 1.0/fratio;
            const Real_type c20 = 0.25 * (clight / csound) * r_fratio;
            const Complex_type ireal(0.0, 1.0);

            COUPLE_Functor kernel(imin, imax, jmin, jmax, 
                                  t0, t1, t2,
                                  denac, denlw,
                                  dt, c10, fratio, r_fratio, c20, ireal);


            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<omp_parallel_for_exec>(kmin, kmax, kernel);

            } 
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_FIR)
          case FIR : {  // --> from subset "A"

            loopInit(iloop, stat);

            Real_ptr out = loop_data.array_1D_Real[0];
            Real_ptr in = loop_data.array_1D_Real[1];

            const Index_type coefflen = 16;
            Real_type coeff[coefflen] = { 3.0, -1.0, -1.0, -1.0,
                                          -1.0, 3.0, -1.0, -1.0,
                                          -1.0, -1.0, 3.0, -1.0,
                                          -1.0, -1.0, -1.0, 3.0 };
            const Index_type len_minus_coeff = len - coefflen;

            Index_type val = 0;

            FIR_Functor kernel(in, out, coeff, coefflen);


            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<omp_parallel_for_exec>(0, len_minus_coeff, kernel);

               val = isamp;

            }
            TIMER_STOP(ltimer);

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

#if defined(COMPILE_INIT3)
          case INIT3 : {   // --> from subset "B"

            loopInit(iloop, stat);

            Real_ptr out1 = loop_data.array_1D_Real[0];
            Real_ptr out2 = loop_data.array_1D_Real[1];
            Real_ptr out3 = loop_data.array_1D_Real[2];
            Real_ptr in1 = loop_data.array_1D_Real[3];
            Real_ptr in2 = loop_data.array_1D_Real[4];

            INIT3_Functor kernel(out1, out2, out3,
                                 in1, in2);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<omp_parallel_for_exec>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_MULADDSUB)
          case MULADDSUB : {   // --> from subset "B"

            loopInit(iloop, stat);

            Real_ptr out1 = loop_data.array_1D_Real[0];
            Real_ptr out2 = loop_data.array_1D_Real[1];
            Real_ptr out3 = loop_data.array_1D_Real[2];
            Real_ptr in1 = loop_data.array_1D_Real[3];
            Real_ptr in2 = loop_data.array_1D_Real[4];

            MULADDSUB_Functor kernel(out1, out2, out3,
                                     in1, in2);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<omp_parallel_for_exec>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_IF_QUAD)
          case IF_QUAD : {  // --> from subset "B"

            loopInit(iloop, stat);

            Real_ptr a = loop_data.array_1D_Real[0];
            Real_ptr b = loop_data.array_1D_Real[1];
            Real_ptr c = loop_data.array_1D_Real[2];
            Real_ptr x1 = loop_data.array_1D_Real[3];
            Real_ptr x2 = loop_data.array_1D_Real[4];

            IF_QUAD_Functor kernel(a, b, c, x1, x2);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<omp_parallel_for_exec>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_TRAP_INT)
          case TRAP_INT : {   // --> from subset "B"

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

            TRAP_INT_Functor_OMP kernel(&sumx, x0, h, y, xp, yp);

            Real_type val = 0;

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<omp_parallel_for_exec>(0, len, kernel);

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

#if defined(COMPILE_PIC_2D)
          case PIC_2D : {  // --> from subset "C"

            loopInit(iloop, stat);

            //
            // RDH modified kernel to fix zero-based indexing error.
            //

            Real_ptr* p = loop_data.array_2D_Nx25_Real[0];
            Real_ptr* b = loop_data.array_2D_Nx25_Real[1];
            Real_ptr* c = loop_data.array_2D_Nx25_Real[2];

            Real_ptr y = loop_data.array_1D_Real[0];
            Real_ptr z = loop_data.array_1D_Real[1];

            Index_type* e = loop_data.array_1D_Indx[0];
            Index_type* f = loop_data.array_1D_Indx[1];

            Real_ptr* h = loop_data.array_2D_64x64_Real[0];

            PIC_2D_Functor_OMP kernel(p, b, c,
                                      y, z,
                                      e, f,
                                      h);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<omp_parallel_for_exec>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

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

#endif  // if COMPILE_FUNCTOR_VARIANTS && COMPILE_OMP_VARIANTS

}
