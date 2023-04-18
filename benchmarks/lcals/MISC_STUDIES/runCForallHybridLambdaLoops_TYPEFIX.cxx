//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Source file containing LCALS "C" subset forall hybrid lambda loops
//     with type information "fix" 
//

#include "LCALSSuite.hxx"
#include "LCALSTraversalMethods.hxx"

#include<cstdlib>
#include<iostream>
#include<cmath>


void runCForallHybridLambdaLoops_TYPEFIX( std::vector<LoopStat>& loop_stats,
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

#if defined(COMPILE_HYDRO_1D)
          case HYDRO_1D : {

            loopInit(iloop, stat);

            Real_ptr x = loop_data.array_1D_Real[0];  
            Real_ptr y = loop_data.array_1D_Real[1];  
            Real_ptr z = loop_data.array_1D_Real[2];  

            const Real_type q = loop_data.scalar_Real[0];
            const Real_type r = loop_data.scalar_Real[1];
            const Real_type t = loop_data.scalar_Real[2];

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  Real_ptr ty = y;
                  Real_ptr tz = z;
                  tx[k] = q + ty[k]*( r*tz[k+10] + t*tz[k+11] );
               } ); 
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
            unsigned is_count = 0;

//  Loop execution timing starts here...
            TIMER_START(ltimer);

            ii = len;
            do {
               ii /= 2;
               is_count++;     
            } while ( ii>0 ); 

            std::vector<HybridIndexSet> his(is_count);

            is_count = 0;
            ii = len;
            ipntp = 0;
            do {
              ipnt = ipntp;
              ipntp += ii;
              ii /= 2;
              his[is_count++].addRangeStrideIndices(ipnt+1, ipntp, 2);
            } while ( ii>0 );

//  Loop execution starts here...
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

              is_count = 0;
              ii = len;
              ipntp = 0;
              do {
                ipnt = ipntp;
                ipntp += ii;
                ii /= 2;
                i = ipntp ;
                forall<exec_policy>(his[is_count++], 
                [&] (Index_type k) {
                   i++;
                   Real_ptr tx = x;
                   Real_ptr tv = v;
                   tx[i] = tx[k] - tv[k  ]*tx[k-1] - tv[k+1]*tx[k+1];
                } );
              } while ( ii>0 );

            }
//  Loop execution ends here...

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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               q = 0.0;
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  Real_ptr tz = z;
                  q += tz[k]*tx[k];
               } );

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

            Index_type lw;
            Real_type temp;

            HybridIndexSet his;
            his.addRangeStrideIndices(4, len, 5);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

              Index_type m = ( 1001-7 )/2;
              for ( Index_type k=6 ; k<1001 ; k=k+m ) {
                 lw = k - 6;
                 temp = x[k-1];
                 forall<exec_policy>(his, 
                 [&] (Index_type j) {
                    Real_ptr tx = x;
                    Real_ptr ty = y;
                    temp -= tx[lw]*ty[j];
                    lw++;
                 } );
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

            HybridIndexSet his;
            his.addRangeIndices(1, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type i) {
                  Real_ptr tx = x;
                  Real_ptr ty = y;
                  Real_ptr tz = z;
                  tx[i] = tz[i]*( ty[i] - tx[i-1] );
               } );
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  Real_ptr ty = y;
                  Real_ptr tz = z;
                  Real_ptr tu = u;
                  tx[k] = tu[k] + r*( tz[k] + r*ty[k] ) +
                        t*( tu[k+3] + r*( tu[k+2] + r*tu[k+1] ) +
                           t*( tu[k+6] + q*( tu[k+5] + q*tu[k+4] ) ) );
               } );
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

            Index_type nl1 = 0;
            Index_type nl2 = 1;
            Index_type kx;

            HybridIndexSet his;
            his.addRangeIndices(1, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for (kx=1 ; kx<3 ; kx++ ) {
                  forall<exec_policy>(his, 
                  [&] (Index_type ky) {
                     Real_ptr tdu1 = du1;
                     Real_ptr tdu2 = du2;
                     Real_ptr tdu3 = du3;
                     Real_ptr** tu1 = u1;
                     Real_ptr** tu2 = u2;
                     Real_ptr** tu3 = u3;
                     tdu1[ky] = tu1[nl1][ky+1][kx] - tu1[nl1][ky-1][kx];
                     tdu2[ky] = tu2[nl1][ky+1][kx] - tu2[nl1][ky-1][kx];
                     tdu3[ky] = tu3[nl1][ky+1][kx] - tu3[nl1][ky-1][kx];
                     tu1[nl2][ky][kx]=
                        tu1[nl1][ky][kx]+a11*tdu1[ky]+a12*tdu2[ky]+a13*tdu3[ky] + sig*
                        (tu1[nl1][ky][kx+1]-2.0*tu1[nl1][ky][kx]+tu1[nl1][ky][kx-1]);
                        tu2[nl2][ky][kx]=
                        tu2[nl1][ky][kx]+a21*tdu1[ky]+a22*tdu2[ky]+a23*tdu3[ky] + sig*
                        (tu2[nl1][ky][kx+1]-2.0*tu2[nl1][ky][kx]+tu2[nl1][ky][kx-1]);
                     tu3[nl2][ky][kx]=
                        tu3[nl1][ky][kx]+a31*tdu1[ky]+a32*tdu2[ky]+a33*tdu3[ky] + sig*
                        (tu3[nl1][ky][kx+1]-2.0*tu3[nl1][ky][kx]+tu3[nl1][ky][kx-1]);
                  } );
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type i) {
                  Real_ptr* tpx = px;
                  tpx[i][0] = dm28*tpx[i][12] + dm27*tpx[i][11] + dm26*tpx[i][10] +
                       dm25*tpx[i][ 9] + dm24*tpx[i][ 8] + dm23*tpx[i][ 7] +
                       dm22*tpx[i][ 6] + c0*( tpx[i][ 4] + tpx[i][ 5]) + tpx[i][ 2];
               } );
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type i) {
                  Real_ptr* tpx = px;
                  Real_ptr* tcx = cx;
                  Real_type ar, br, cr;
                  ar        =      tcx[i][ 4];
                  br        = ar - tpx[i][ 4];
                  tpx[i][ 4] = ar;
                  cr        = br - tpx[i][ 5];
                  tpx[i][ 5] = br;
                  ar        = cr - tpx[i][ 6];
                  tpx[i][ 6] = cr;
                  br        = ar - tpx[i][ 7];
                  tpx[i][ 7] = ar;
                  cr        = br - tpx[i][ 8];
                  tpx[i][ 8] = br;
                  ar        = cr - tpx[i][ 9];
                  tpx[i][ 9] = cr;
                  br        = ar - tpx[i][10];
                  tpx[i][10] = ar;
                  cr        = br - tpx[i][11];
                  tpx[i][11] = br;
                  tpx[i][13] = cr - tpx[i][12];
                  tpx[i][12] = cr;
               } );
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

            HybridIndexSet his;
            his.addRangeIndices(1, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               x[0] = y[0];
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  Real_ptr ty = y;
                  tx[k] = tx[k-1] + ty[k];
               } );

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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  Real_ptr ty = y;
                  tx[k] = ty[k+1] - ty[k];
               } );
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
    
            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type ip) {
                  Real_ptr* tp = p;
                  Real_ptr* tb = b;
                  Real_ptr* tc = c;
                  Real_ptr ty = y;
                  Real_ptr tz = z;
                  Real_ptr* th = h;
                  Index_type i1, j1, i2, j2;
                  i1 = (Index_type) tp[ip][0];
                  j1 = (Index_type) tp[ip][1];
                  i1 &= 64-1;
                  j1 &= 64-1;
                  tp[ip][2] += tb[j1][i1];
                  tp[ip][3] += tc[j1][i1];
                  tp[ip][0] += tp[ip][2];
                  tp[ip][1] += tp[ip][3];
                  i2 = (Index_type) tp[ip][0];
                  j2 = (Index_type) tp[ip][1];
                  i2 = ( i2 & 64-1 ) ;
                  j2 = ( j2 & 64-1 ) ;
                  tp[ip][0] += ty[i2+32];
                  tp[ip][1] += tz[j2+32];
                  i2 += e[i2+32];
                  j2 += f[j2+32];
                  th[j2][i2] += 1.0;
               } );
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


            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<exec_policy>(his, 
               [&] (Index_type k) {
                   Real_ptr tvx = vx;
                   Real_ptr txx = xx;
                   Real_ptr txi = xi;
                   Real_ptr tex = ex;
                   Real_ptr tex1 = ex1;
                   Real_ptr tdex = dex;
                   Real_ptr tdex1 = dex1;
                   tvx[k] = 0.0;
                   txx[k] = 0.0;
                   ix[k] = (Index_type) grd[k];
                   txi[k] = (Real_type) ix[k];
                   tex1[k] = tex[ ix[k] - 1 ];
                   tdex1[k] = tdex[ ix[k] - 1 ];
               } );

               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tvx = vx;
                  Real_ptr txx = xx;
                  Real_ptr txi = xi;
                  Real_ptr tex = ex;
                  Real_ptr tex1 = ex1;
                  Real_ptr tdex = dex;
                  Real_ptr tdex1 = dex1;
                  Real_ptr trx = rx;
                  tvx[k] = tvx[k] + tex1[k] + ( txx[k] - txi[k] )*tdex1[k];
                  txx[k] = txx[k] + tvx[k]  + flx;
                  ir[k] = (Index_type) xx[k];
                  trx[k] = txx[k] - ir[k];
                  ir[k] = ( ir[k] & (2048-1) ) + 1;
                  txx[k] = trx[k] + ir[k];
               } );

               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr trh = rh;
                  Real_ptr trx = rx;
                  trh[ ir[k]-1 ] += 1.0 - trx[k];
                  trh[ ir[k]   ] += trx[k];
               } );

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

            HybridIndexSet his;
            his.addRangeIndices(1, jn);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for ( k=1 ; k<kn ; k++ ) {
                  forall<exec_policy>(his,
                  [&] (Index_type j) {
                     Real_ptr* tza = za;
                     Real_ptr* tzb = zb;
                     Real_ptr* tzm = zm;
                     Real_ptr* tzp = zp;
                     Real_ptr* tzq = zq;
                     Real_ptr* tzr = zr;
                     tza[k][j] = ( tzp[k+1][j-1] +tzq[k+1][j-1] -tzp[k][j-1] -tzq[k][j-1] )*
                                 ( tzr[k][j] +tzr[k][j-1] ) / ( tzm[k][j-1] +tzm[k+1][j-1]);
                     tzb[k][j] = ( tzp[k][j-1] +tzq[k][j-1] -tzp[k][j] -tzq[k][j] ) *
                                 ( tzr[k][j] +tzr[k-1][j] ) / ( tzm[k][j] +tzm[k][j-1]);
                  } );
               }

               for ( k=1 ; k<kn ; k++ ) {
                  forall<exec_policy>(his,
                  [&] (Index_type j) {
                     Real_ptr* tzu = zu;
                     Real_ptr* tzv = zv;
                     Real_ptr* tza = za;
                     Real_ptr* tzb = zb;
                     Real_ptr* tzr = zr;
                     Real_ptr* tzz = zz;
                     tzu[k][j] += s*( tza[k][j]   *( tzz[k][j] - tzz[k][j+1] ) -
                                     tza[k][j-1] *( tzz[k][j] - tzz[k][j-1] ) -
                                     tzb[k][j]   *( tzz[k][j] - tzz[k-1][j] ) +
                                     tzb[k+1][j] *( tzz[k][j] - tzz[k+1][j] ) );
                     tzv[k][j] += s*( tza[k][j]   *( tzr[k][j] - tzr[k][j+1] ) -
                                     tza[k][j-1] *( tzr[k][j] - tzr[k][j-1] ) -
                                     tzb[k][j]   *( tzr[k][j] - tzr[k-1][j] ) +
                                     tzb[k+1][j] *( tzr[k][j] - tzr[k+1][j] ) );
                  } );
               }

               for ( k=1 ; k<kn ; k++ ) {
                  forall<exec_policy>(his,
                  [&] (Index_type j) {
                     Real_ptr* tzr = zr;
                     Real_ptr* tzz = zz;
                     Real_ptr* tzu = zu;
                     Real_ptr* tzv = zv;
                     Real_ptr* tzrout = zrout;
                     Real_ptr* tzzout = zzout;
                     tzrout[k][j] = tzr[k][j] + t*tzu[k][j];
                     tzzout[k][j] = tzz[k][j] + t*tzv[k][j];
                  } );
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

            Index_type kb5i = 0;

            HybridIndexSet hisa;
            hisa.addRangeIndices(0, len);
            HybridIndexSet hisb;
            hisb.addRangeIndices(1, len+1);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               forall<exec_policy>(hisa,
               [&] (Index_type k) {
                  Real_ptr tb5 = b5;
                  Real_ptr tsa = sa;
                  Real_ptr tsb = sb;
                  tb5[k+kb5i] = tsa[k] + stb5*tsb[k];
                  stb5 = tb5[k+kb5i] - stb5;
               } );

               forall<exec_policy>(hisb,
               [&] (Index_type i) {
                  Real_ptr tb5 = b5;
                  Real_ptr tsa = sa;
                  Real_ptr tsb = sb;
                  Index_type k = len - i ;
                  tb5[k+kb5i] = tsa[k] + stb5*tsb[k];
                  stb5 = tb5[k+kb5i] - stb5;
               } );

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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  Real_ptr ty = y;
                  Real_ptr tz = z;
                  Real_ptr tu = u;
                  Real_ptr tv = v;
                  Real_ptr tw = w;
                  Real_ptr tg = g;
                  Real_ptr txx = xx;
                  Real_ptr tvx = vx;
                  Real_type di = ty[k] - tg[k] / ( txx[k] + dk );
                  Real_type dn = 0.2;
                  if ( di ) {
                     dn = tz[k]/di ;
                     if ( t < dn ) dn = t;
                     if ( s > dn ) dn = s;
                  }
                  tx[k] = ((tw[k] + tv[k]*dn )* txx[k] + tu[k]) / (tvx[k] + tv[k]*dn);
                  txx[k+1] = ( tx[k] - txx[k] )* dn + txx[k];
               } );
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for ( k=0 ; k<25 ; k++ ) {
                  for ( i=0 ; i<25 ; i++ ) {
                     forall<exec_policy>(his, 
                     [&] (Index_type j) {
                        Real_ptr* tpx = px;
                        Real_ptr* tcx = cx;
                        Real_ptr* tvy = vy;
                        tpx[j][i] += tvy[k][i] * tcx[j][k];
                     } );
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

            HybridIndexSet his;
            his.addRangeIndices(0, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  Real_ptr ty = y;
                  Real_ptr tu = u;
                  Real_ptr tv = v;
                  Real_ptr tw = w;
                  ty[k] = tu[k] / tv[k];
                  tw[k] = tx[k] / ( exp( ty[k] ) -1.0 );
               } ); 
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

            HybridIndexSet his;
            his.addRangeIndices(1, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               for ( j=1 ; j<6 ; j++ ) {
                  forall<exec_policy>(his, 
                  [&] (Index_type k) {
                     Real_ptr* tza = za;
                     Real_ptr* tzb = zb;
                     Real_ptr* tzr = zr;
                     Real_ptr* tzu = zu;
                     Real_ptr* tzv = zv;
                     Real_ptr* tzz = zz;
                     Real_type qa = tza[j+1][k]*tzr[j][k] + tza[j-1][k]*tzb[j][k] +
                          tza[j][k+1]*tzu[j][k] + tza[j][k-1]*tzv[j][k] + tzz[j][k];
                     tza[j][k] += 0.175*( qa - tza[j][k] );
                  } );
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

            HybridIndexSet his;
            his.addRangeIndices(1, len);

            TIMER_START(ltimer);
            for (SampIndex_type isamp = 0; isamp < num_samples; ++isamp) {

               x[len/2] = (-1.0e+10)*isamp;
               m = 0;
               forall<exec_policy>(his, 
               [&] (Index_type k) {
                  Real_ptr tx = x;
                  if ( tx[k] < tx[m] ) m = k;
               } );

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

#endif  // if COMPILE_LAMBDA_VARIANTS

}
