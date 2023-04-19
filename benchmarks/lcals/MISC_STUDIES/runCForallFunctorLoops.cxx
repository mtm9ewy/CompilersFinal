//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "C" subset forall functor loops
//

#include "LCALSSuite.hxx"
#include "LCALSTraversalMethods.hxx"

#include "FunctorKernelsC.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runCForallFunctorLoops( std::vector<LoopStat>& loop_stats,
                             bool run_loop[],
                             LoopLength ilength )
{
   LoopSuiteRunInfo& loop_suite_run_info = getLoopSuiteRunInfo();
   LoopData& loop_data = getLoopData();

#if defined(COMPILE_FUNCTOR_VARIANTS)

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

            const Real_type q = loop_data.scalar_Real[0];
            const Real_type r = loop_data.scalar_Real[1];
            const Real_type t = loop_data.scalar_Real[2];

            HYDRO_1D_Functor kernel(x, y, z, 
                                    q, r, t);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_ICCG)
          case ICCG : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Nx4_Real[0];
            Real_ptr v = loop_data.array_1D_Nx4_Real[1];

            Index_type ii, ipnt, ipntp, i;  

            ICCG_Functor kernel(x, v, &i); 

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

              ii = len;
              ipntp = 0;
              do {
                ipnt = ipntp;
                ipntp += ii;
                ii /= 2;
                i = ipntp;
                forall<exec_policy>(ipnt+1, ipntp, 2, kernel);
              } while ( ii>0 );

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }

#endif

#if defined(COMPILE_INNER_PROD)
          case INNER_PROD : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr z = loop_data.array_1D_Real[1];

            Real_type q = 0.0;
            Real_type val = 0.0;

            INNER_PROD_Functor kernel(x, z, &q);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               q = 0.0;
               forall<exec_policy>(0, len, kernel); 

               val = q*isamp;
            }
            TIMER_STOP(ltimer);

            //
            // RDH added this. Without it compiler may optimize out
            // outer sampling loop because value of q was not used.
            //
            loop_data.scalar_Real[0] =
               (val + 0.00123) / (val - 0.00123);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_BAND_LIN_EQ)
          case BAND_LIN_EQ : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];

            Real_type temp;
            Index_type lw;

            BAND_LIN_EQ_Functor kernel(x, y, &temp, &lw);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

              Index_type m = ( 1001-7 )/2;
              for ( Index_type k=6 ; k<1001 ; k=k+m ) {
                 lw = k - 6;
                 temp = x[k-1];
                 forall<exec_policy>(4, len, 5, kernel);
                 x[k-1] = y[4]*temp;
              }

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_TRIDIAG_ELIM)
          case TRIDIAG_ELIM : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr z = loop_data.array_1D_Real[2];

            TRIDIAG_ELIM_Functor kernel(x, y, z);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(1, len, kernel); 
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_EOS)
          case EOS : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];  
            Real_ptr y = loop_data.array_1D_Real[1];  
            Real_ptr z = loop_data.array_1D_Real[2];  
            Real_ptr u = loop_data.array_1D_Real[3];  

            const Real_type q = loop_data.scalar_Real[0];
            const Real_type r = loop_data.scalar_Real[1];
            const Real_type t = loop_data.scalar_Real[2];

            EOS_Functor kernel(x, y, z, u,
                               q, r, t);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_ADI)
          case ADI : {

            loopInit(iloop, stat);

            Real_ptr du1 = loop_data.array_1D_Real[0];
            Real_ptr du2 = loop_data.array_1D_Real[1];
            Real_ptr du3 = loop_data.array_1D_Real[2];

            Real_ptr** u1 = loop_data.array_3D_2xNx4_Real[0];
            Real_ptr** u2 = loop_data.array_3D_2xNx4_Real[1];
            Real_ptr** u3 = loop_data.array_3D_2xNx4_Real[2];

            const Real_type sig = loop_data.scalar_Real[0];
            const Real_type a11 = loop_data.scalar_Real[1];
            const Real_type a12 = loop_data.scalar_Real[2];
            const Real_type a13 = loop_data.scalar_Real[3];
            const Real_type a21 = loop_data.scalar_Real[4];
            const Real_type a22 = loop_data.scalar_Real[5];
            const Real_type a23 = loop_data.scalar_Real[6];
            const Real_type a31 = loop_data.scalar_Real[7];
            const Real_type a32 = loop_data.scalar_Real[8];
            const Real_type a33 = loop_data.scalar_Real[9];

            Index_type kx;

            ADI_Functor kernel(du1, du2, du3,
                               u1, u2, u3, 
                               sig,
                               a11, a12, a13,  
                               a21, a22, a23,  
                               a31, a32, a33,
                               &kx);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for ( kx=1 ; kx<3 ; kx++ ) {
                  forall<exec_policy>(1, len, kernel);
               }

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_INT_PREDICT)
          case INT_PREDICT : {

            loopInit(iloop, stat);

            Real_ptr* px = loop_data.array_2D_Nx25_Real[0];

            const Real_type dm22 = loop_data.scalar_Real[0];
            const Real_type dm23 = loop_data.scalar_Real[1];
            const Real_type dm24 = loop_data.scalar_Real[2];
            const Real_type dm25 = loop_data.scalar_Real[3];
            const Real_type dm26 = loop_data.scalar_Real[4];
            const Real_type dm27 = loop_data.scalar_Real[5];
            const Real_type dm28 = loop_data.scalar_Real[6];
            const Real_type c0 = loop_data.scalar_Real[7];

            INT_PREDICT_Functor kernel(px,
                                       dm22, dm23, dm24,
                                       dm25, dm26, dm27,
                                       dm28, c0);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_DIFF_PREDICT)
          case DIFF_PREDICT : {

            loopInit(iloop, stat);

            Real_ptr* px = loop_data.array_2D_Nx25_Real[0];
            Real_ptr* cx = loop_data.array_2D_Nx25_Real[1];

            DIFF_PREDICT_Functor kernel(px, cx);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_FIRST_SUM)
          case FIRST_SUM : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];

            FIRST_SUM_Functor kernel(x, y);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               x[0] = y[0];
               forall<exec_policy>(1, len, kernel);

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_FIRST_DIFF)
          case FIRST_DIFF : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];

            FIRST_DIFF_Functor kernel(x, y);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_PIC_2D)
          case PIC_2D : {

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

            PIC_2D_Functor kernel(p, b, c,
                                  y, z,
                                  e, f,
                                  h);
    
            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_PIC_1D)
          case PIC_1D : {

            loopInit(iloop, stat);

            Real_ptr vx = loop_data.array_1D_Real[0];
            Real_ptr xx = loop_data.array_1D_Real[1];
            Real_ptr xi = loop_data.array_1D_Real[2];
            Real_ptr ex = loop_data.array_1D_Real[3];
            Real_ptr ex1 = loop_data.array_1D_Real[4];
            Real_ptr dex = loop_data.array_1D_Real[5];
            Real_ptr dex1 = loop_data.array_1D_Real[6];
            Real_ptr rh = loop_data.array_1D_Real[7];
            Real_ptr rx = loop_data.array_1D_Real[8];

            const Real_type flx = loop_data.scalar_Real[0];

            Index_type* ix = loop_data.array_1D_Indx[2];
            Index_type* ir = loop_data.array_1D_Indx[3];
            Index_type* grd = loop_data.array_1D_Indx[4];


            PIC_1D_FunctorA kernela(vx, xx, xi,
                                    ex, ex1,
                                    dex, dex1, 
                                    ix, grd);

            PIC_1D_FunctorB kernelb(vx, xx, xi,
                                    ex1,
                                    dex1, 
                                    rx, 
                                    ir, 
                                    flx);

            PIC_1D_FunctorC kernelc(rh, rx,
                                    ir);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<exec_policy>(0, len, kernela);
               forall<exec_policy>(0, len, kernelb);
               forall<exec_policy>(0, len, kernelc);

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_HYDRO_2D)
          case HYDRO_2D : {

            loopInit(iloop, stat);

            Real_ptr* za = loop_data.array_2D_7xN_Real[0];
            Real_ptr* zb = loop_data.array_2D_7xN_Real[1];
            Real_ptr* zm = loop_data.array_2D_7xN_Real[2];
            Real_ptr* zp = loop_data.array_2D_7xN_Real[3];
            Real_ptr* zq = loop_data.array_2D_7xN_Real[4];
            Real_ptr* zr = loop_data.array_2D_7xN_Real[5];
            Real_ptr* zu = loop_data.array_2D_7xN_Real[6];
            Real_ptr* zv = loop_data.array_2D_7xN_Real[7];
            Real_ptr* zz = loop_data.array_2D_7xN_Real[8];

            Real_ptr* zrout = loop_data.array_2D_7xN_Real[9];
            Real_ptr* zzout = loop_data.array_2D_7xN_Real[10];

            const Real_type t = 0.0037;
            const Real_type s = 0.0041;

            Index_type kn = 6;
            Index_type jn = len;
            Index_type k;

            HYDRO_2D_FunctorA kernela(za, zb, zm, zp, zq, zr, &k);

            HYDRO_2D_FunctorB kernelb(zu, zv, za, zb, zr, zz, s, &k); 

            HYDRO_2D_FunctorC kernelc(zrout, zzout, zr, zz, zu, zv, t, &k); 

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for ( k=1 ; k<kn ; k++ ) {
                  forall<exec_policy>(1, jn, kernela);
               }
              
               for ( k=1 ; k<kn ; k++ ) {
                  forall<exec_policy>(1, jn, kernelb);
               }

               for ( k=1 ; k<kn ; k++ ) {
                  forall<exec_policy>(1, jn, kernelc);
               }

            } 
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_GEN_LIN_RECUR)
          case GEN_LIN_RECUR : {

            loopInit(iloop, stat);

            Real_ptr b5 = loop_data.array_1D_Real[0];
            Real_ptr sa = loop_data.array_1D_Real[1];
            Real_ptr sb = loop_data.array_1D_Real[2];
            Real_type stb5 = loop_data.scalar_Real[0];

            GEN_LIN_RECUR_FunctorA kernela(b5, sa, sb, &stb5);

            GEN_LIN_RECUR_FunctorB kernelb(b5, sa, sb, &stb5, len); 
            
            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<exec_policy>(0, len, kernela);
               forall<exec_policy>(1, len+1, kernelb);

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_DISC_ORD)
          case DISC_ORD : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr z = loop_data.array_1D_Real[2];
            Real_ptr u = loop_data.array_1D_Real[3];
            Real_ptr v = loop_data.array_1D_Real[4];
            Real_ptr w = loop_data.array_1D_Real[5];
            Real_ptr g = loop_data.array_1D_Real[6];
            Real_ptr xx = loop_data.array_1D_Real[7];
            Real_ptr vx = loop_data.array_1D_Real[9];
            const Real_type s = loop_data.scalar_Real[0];
            const Real_type t = loop_data.scalar_Real[1];
            const Real_type dk = loop_data.scalar_Real[2];

            DISC_ORD_Functor kernel(x, y, z, 
                                    u, v, w, g,
                                    xx, vx,
                                    s, t, dk);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_MAT_X_MAT)
          case MAT_X_MAT : {

            loopInit(iloop, stat);

            Real_ptr* px = loop_data.array_2D_Nx25_Real[0];
            Real_ptr* cx = loop_data.array_2D_Nx25_Real[1];
            Real_ptr* vy = loop_data.array_2D_64x64_Real[0];

            Index_type k, i;

            MAT_X_MAT_Functor kernel(px, cx, vy, &k, &i);        

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for ( k=0 ; k<25 ; k++ ) {
                  for ( i=0 ; i<25 ; i++ ) {
                     forall<exec_policy>(0, len, kernel);
                  }
               }

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_PLANCKIAN)
          case PLANCKIAN : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];
            Real_ptr y = loop_data.array_1D_Real[1];
            Real_ptr u = loop_data.array_1D_Real[2];
            Real_ptr v = loop_data.array_1D_Real[3];
            Real_ptr w = loop_data.array_1D_Real[4];

            Real_type expmax = 20.0;
            u[len-1] = 0.99*expmax*v[len-1];

            PLANCKIAN_Functor kernel(x, y, u, v, w);            

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(0, len, kernel);
            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_IMP_HYDRO_2D)
          case IMP_HYDRO_2D : {

            loopInit(iloop, stat);

            Real_ptr* za = loop_data.array_2D_7xN_Real[0];
            Real_ptr* zb = loop_data.array_2D_7xN_Real[1];
            Real_ptr* zr = loop_data.array_2D_7xN_Real[2];
            Real_ptr* zu = loop_data.array_2D_7xN_Real[3];
            Real_ptr* zv = loop_data.array_2D_7xN_Real[4];
            Real_ptr* zz = loop_data.array_2D_7xN_Real[5];

            Index_type j; 

            IMP_HYDRO_2D_Functor kernel(za, zb, zr, zu, zv, zz, &j);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for ( j=1 ; j<6 ; j++ ) {
                  forall<exec_policy>(1, len, kernel);
               }

            }
            TIMER_STOP(ltimer);

            loopFinalize(iloop, stat, ilength);

            break;
          }
#endif

#if defined(COMPILE_FIND_FIRST_MIN)
          case FIND_FIRST_MIN : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];

            Index_type m = 0;
            Index_type val = 0;

            FIND_FIRST_MIN_Functor kernel(x, &m);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               x[len/2] = (-1.0e+10)*isamp;
               m = 0;
               forall<exec_policy>(1, len, kernel);

               val = isamp;
            }
            TIMER_STOP(ltimer);

            //
            // RDH added this. Without it compiler optimized loop out
            // because value of m was not used.
            //
            loop_data.scalar_Real[0] = 0.123*val*m;

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

#endif  // if COMPILE_FUNCTOR_VARIANTS

}
