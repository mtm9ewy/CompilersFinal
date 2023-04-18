//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Header file with functors for LCALS "C" subset
//

#ifndef LCALSFunctorKernelsC_HXX
#define LCALSFunctorKernelsC_HXX

#include "LCALSParams.hxx"

#include<cmath>


#if defined(COMPILE_HYDRO_1D)
/*
 *******************************************************************
 *   Kernel 1 -- hydro fragment
 *******************************************************************
 *       DO 1 L = 1,Loop
 *       DO 1 k = 1,n
 *  1       X(k)= Q + Y(k)*(R*ZX(k+10) + T*ZX(k+11))
 */
class HYDRO_1D_Functor
{
public:
   HYDRO_1D_Functor(Real_ptr x, Real_ptr y, Real_ptr z,
                    const Real_type q, const Real_type r, const Real_type t)
   : m_x(x), m_y(y), m_z(z), 
     m_q(q), m_r(r), m_t(t) { ; }

   void operator() (Index_type k)
   {
      m_x[k] = m_q + m_y[k]*( m_r*m_z[k+10] + m_t*m_z[k+11] );
   }

   Real_ptr m_x; 
   Real_ptr m_y; 
   Real_ptr m_z;
   const Real_type m_q; 
   const Real_type m_r; 
   const Real_type m_t;
};
#endif

#if defined(COMPILE_ICCG)
/*
 *******************************************************************
 *   Kernel 2 -- ICCG excerpt (Incomplete Cholesky Conj. Gradient)
 *******************************************************************
 *    DO 200  L= 1,Loop
 *        II= n
 *     IPNTP= 0
 *222   IPNT= IPNTP
 *     IPNTP= IPNTP+II
 *        II= II/2
 *         i= IPNTP+1
 CDIR$ IVDEP
 *    DO 2 k= IPNT+2,IPNTP,2
 *         i= i+1
 *  2   X(i)= X(k) - V(k)*X(k-1) - V(k+1)*X(k+1)
 *        IF( II.GT.1) GO TO 222
 *200 CONTINUE
 */
class ICCG_Functor
{
public:
   ICCG_Functor(Real_ptr x, Real_ptr v,
                Index_type* i)
   : m_x(x), m_v(v), 
     m_i(i) { ; }

   void operator() (Index_type k)
   {
      (*m_i)++;
      m_x[(*m_i)] = m_x[k] - m_v[k  ]*m_x[k-1] - m_v[k+1]*m_x[k+1];
   }

   Real_ptr m_x;
   Real_ptr m_v;
   Index_type* m_i;
};
#endif

#if defined(COMPILE_INNER_PROD)
/*
 *******************************************************************
 *   Kernel 3 -- inner product
 *******************************************************************
 *    DO 3 L= 1,Loop
 *         Q= 0.0
 *    DO 3 k= 1,n
 *  3      Q= Q + Z(k)*X(k)
 */
class INNER_PROD_Functor
{
public:
   INNER_PROD_Functor(Real_ptr x, Real_ptr z,
                      Real_type* q)
   : m_x(x), m_z(z),
     m_q(q) { ; }

   void operator() (Index_type k)
   {
      *m_q += m_z[k]*m_x[k];
   }

   Real_ptr m_x;
   Real_ptr m_z;
   Real_type* m_q;
};
#endif

#if defined(COMPILE_BAND_LIN_EQ)
/*
 *******************************************************************
 *   Kernel 4 -- banded linear equations
 *******************************************************************
 *            m= (1001-7)/2
 *    DO 444  L= 1,Loop
 *    DO 444  k= 7,1001,m
 *           lw= k-6
 *         temp= X(k-1)
 CDIR$ IVDEP
 *    DO   4  j= 5,n,5
 *       temp  = temp   - XZ(lw)*Y(j)
 *  4        lw= lw+1
 *       X(k-1)= Y(5)*temp
 *444 CONTINUE
 */
class BAND_LIN_EQ_Functor
{
public:
   BAND_LIN_EQ_Functor(Real_ptr x, Real_ptr y,
                       Real_type* temp, Index_type* lw)
   : m_x(x), m_y(y),
     m_temp(temp),
     m_lw(lw) { ; }

   void operator() (Index_type j)
   {
      *m_temp -= m_x[*m_lw]*m_y[j];
      (*m_lw)++;
   }

   Real_ptr m_x;
   Real_ptr m_y;
   Real_type* m_temp;
   Index_type* m_lw;
};
#endif

#if defined(COMPILE_TRIDIAG_ELIM)
/*
 *******************************************************************
 *   Kernel 5 -- tri-diagonal elimination, below diagonal
 *******************************************************************
 *    DO 5 L = 1,Loop
 *    DO 5 i = 2,n
 *  5    X(i)= Z(i)*(Y(i) - X(i-1))
 */
class TRIDIAG_ELIM_Functor
{
public:
   TRIDIAG_ELIM_Functor(Real_ptr x, Real_ptr y, Real_ptr z)
   : m_x(x), m_y(y), m_z(z) { ; }

   void operator() (Index_type i)
   {
      m_x[i] = m_z[i]*( m_y[i] - m_x[i-1] );
   }

   Real_ptr m_x;
   Real_ptr m_y;
   Real_ptr m_z;
};
#endif

#if defined(COMPILE_EOS)
/*
 *******************************************************************
 *   Kernel 7 -- equation of state fragment
 *******************************************************************
 *    DO 7 L= 1,Loop
 *    DO 7 k= 1,n
 *      X(k)=     U(k  ) + R*( Z(k  ) + R*Y(k  )) +
 *   .        T*( U(k+3) + R*( U(k+2) + R*U(k+1)) +
 *   .        T*( U(k+6) + Q*( U(k+5) + Q*U(k+4))))
 *  7 CONTINUE
 */
class EOS_Functor
{
public:
   EOS_Functor(Real_ptr x, Real_ptr y, Real_ptr z, Real_ptr u,
               const Real_type q, const Real_type r, const Real_type t)
   : m_x(x), m_y(y), m_z(z), m_u(u),
     m_q(q), m_r(r), m_t(t) { ; }

   void operator() (Index_type k)
   {
      m_x[k] = m_u[k] + m_r*( m_z[k] + m_r*m_y[k] ) +
           m_t*( m_u[k+3] + m_r*( m_u[k+2] + m_r*m_u[k+1] ) +
              m_t*( m_u[k+6] + m_q*( m_u[k+5] + m_q*m_u[k+4] ) ) );
   }

   Real_ptr m_x;
   Real_ptr m_y;
   Real_ptr m_z;
   Real_ptr m_u;
   const Real_type m_q;
   const Real_type m_r;
   const Real_type m_t;
};
#endif

#if defined(COMPILE_ADI)
/*
 *******************************************************************
 *   Kernel 8 -- ADI integration
 *******************************************************************
 *    DO  8      L = 1,Loop
 *             nl1 = 1
 *             nl2 = 2
 *    DO  8     kx = 2,3
 CDIR$ IVDEP
 *    DO  8     ky = 2,n
 *          DU1(ky)=U1(kx,ky+1,nl1)  -  U1(kx,ky-1,nl1)
 *          DU2(ky)=U2(kx,ky+1,nl1)  -  U2(kx,ky-1,nl1)
 *          DU3(ky)=U3(kx,ky+1,nl1)  -  U3(kx,ky-1,nl1)
 *    U1(kx,ky,nl2)=U1(kx,ky,nl1) +A11*DU1(ky) +A12*DU2(ky) +A13*DU3(ky)
 *   .       + SIG*(U1(kx+1,ky,nl1) -2.*U1(kx,ky,nl1) +U1(kx-1,ky,nl1))
 *    U2(kx,ky,nl2)=U2(kx,ky,nl1) +A21*DU1(ky) +A22*DU2(ky) +A23*DU3(ky)
 *   .       + SIG*(U2(kx+1,ky,nl1) -2.*U2(kx,ky,nl1) +U2(kx-1,ky,nl1))
 *    U3(kx,ky,nl2)=U3(kx,ky,nl1) +A31*DU1(ky) +A32*DU2(ky) +A33*DU3(ky)
 *   .       + SIG*(U3(kx+1,ky,nl1) -2.*U3(kx,ky,nl1) +U3(kx-1,ky,nl1))
 *  8 CONTINUE
 */
class ADI_Functor
{
public:
   ADI_Functor(Real_ptr du1, Real_ptr du2, Real_ptr du3,
               Real_ptr** u1, Real_ptr** u2, Real_ptr** u3,
               const Real_type sig, 
               const Real_type a11, const Real_type a12, const Real_type a13,
               const Real_type a21, const Real_type a22, const Real_type a23,
               const Real_type a31, const Real_type a32, const Real_type a33,
               Index_type* kx)
   : m_du1(du1), m_du2(du2), m_du3(du3),
     m_u1(u1), m_u2(u2), m_u3(u3), 
     m_sig(sig), 
     m_a11(a11), m_a12(a12), m_a13(a13),
     m_a21(a21), m_a22(a22), m_a23(a23),
     m_a31(a11), m_a32(a32), m_a33(a33),
     m_kx(kx),
     m_nl1(0), m_nl2(1) { ; }

   void operator() (Index_type ky)
   {
      Index_type tkx = *m_kx; 
      m_du1[ky] = m_u1[m_nl1][ky+1][tkx] - m_u1[m_nl1][ky-1][tkx];
      m_du2[ky] = m_u2[m_nl1][ky+1][tkx] - m_u2[m_nl1][ky-1][tkx];
      m_du3[ky] = m_u3[m_nl1][ky+1][tkx] - m_u3[m_nl1][ky-1][tkx];
      m_u1[m_nl2][ky][tkx]=
         m_u1[m_nl1][ky][tkx]+m_a11*m_du1[ky]+m_a12*m_du2[ky]+m_a13*m_du3[ky] +
         m_sig*
         (m_u1[m_nl1][ky][tkx+1]-2.0*m_u1[m_nl1][ky][tkx]+m_u1[m_nl1][ky][tkx-1]);
      m_u2[m_nl2][ky][tkx]=
         m_u2[m_nl1][ky][tkx]+m_a21*m_du1[ky]+m_a22*m_du2[ky]+m_a23*m_du3[ky] +
         m_sig*
         (m_u2[m_nl1][ky][tkx+1]-2.0*m_u2[m_nl1][ky][tkx]+m_u2[m_nl1][ky][tkx-1]);
      m_u3[m_nl2][ky][tkx]=
         m_u3[m_nl1][ky][tkx]+m_a31*m_du1[ky]+m_a32*m_du2[ky]+m_a33*m_du3[ky] +         m_sig*
         (m_u3[m_nl1][ky][tkx+1]-2.0*m_u3[m_nl1][ky][tkx]+m_u3[m_nl1][ky][tkx-1]); 
   }

   Real_ptr m_du1;
   Real_ptr m_du2;
   Real_ptr m_du3;
   Real_ptr** m_u1;
   Real_ptr** m_u2;
   Real_ptr** m_u3;
   const Real_type m_sig;
   const Real_type m_a11;
   const Real_type m_a12;
   const Real_type m_a13;
   const Real_type m_a21;
   const Real_type m_a22;
   const Real_type m_a23;
   const Real_type m_a31;
   const Real_type m_a32;
   const Real_type m_a33;
   Index_type* m_kx;
   Index_type m_nl1;
   Index_type m_nl2;
};
#endif

#if defined(COMPILE_INT_PREDICT)
/*
 *******************************************************************
 *   Kernel 9 -- integrate predictors
 *******************************************************************
 *    DO 9  L = 1,Loop
 *    DO 9  i = 1,n
 *    PX( 1,i)= DM28*PX(13,i) + DM27*PX(12,i) + DM26*PX(11,i) +
 *   .          DM25*PX(10,i) + DM24*PX( 9,i) + DM23*PX( 8,i) +
 *   .          DM22*PX( 7,i) +  C0*(PX( 5,i) +      PX( 6,i))+ PX( 3,i)
 *  9 CONTINUE
 */
class INT_PREDICT_Functor
{
public:
   INT_PREDICT_Functor(Real_ptr* px,
                       const Real_type dm22, const Real_type dm23, 
                       const Real_type dm24, const Real_type dm25, 
                       const Real_type dm26, const Real_type dm27,
                       const Real_type dm28, const Real_type c0)
   : m_px(px),
     m_dm22(dm22), m_dm23(dm23), 
     m_dm24(dm24), m_dm25(dm25), 
     m_dm26(dm26), m_dm27(dm27),
     m_dm28(dm28),m_c0(c0) { ; }

   void operator() (Index_type i)
   {
      m_px[i][0] = m_dm28*m_px[i][12] + m_dm27*m_px[i][11] + 
                   m_dm26*m_px[i][10] + m_dm25*m_px[i][ 9] + 
                   m_dm24*m_px[i][ 8] + m_dm23*m_px[i][ 7] +
                   m_dm22*m_px[i][ 6] + 
                   m_c0*( m_px[i][ 4] + m_px[i][ 5]) + m_px[i][ 2];
   }

   Real_ptr* m_px;
   const Real_type m_dm22;
   const Real_type m_dm23;
   const Real_type m_dm24;
   const Real_type m_dm25;
   const Real_type m_dm26;
   const Real_type m_dm27;
   const Real_type m_dm28;
   const Real_type m_c0;
};
#endif

#if defined(COMPILE_DIFF_PREDICT)
/*
 *******************************************************************
 *   Kernel 10 -- difference predictors
 *******************************************************************
 *    DO 10  L= 1,Loop
 *    DO 10  i= 1,n
 *    AR      =      CX(5,i)
 *    BR      = AR - PX(5,i)
 *    PX(5,i) = AR
 *    CR      = BR - PX(6,i)
 *    PX(6,i) = BR
 *    AR      = CR - PX(7,i)
 *    PX(7,i) = CR
 *    BR      = AR - PX(8,i)
 *    PX(8,i) = AR
 *    CR      = BR - PX(9,i)
 *    PX(9,i) = BR
 *    AR      = CR - PX(10,i)
 *    PX(10,i)= CR
 *    BR      = AR - PX(11,i)
 *    PX(11,i)= AR
 *    CR      = BR - PX(12,i)
 *    PX(12,i)= BR
 *    PX(14,i)= CR - PX(13,i)
 *    PX(13,i)= CR
 * 10 CONTINUE
 */
class DIFF_PREDICT_Functor
{
public:
   DIFF_PREDICT_Functor(Real_ptr* px, Real_ptr* cx)
   : m_px(px), m_cx(cx) { ; }

   void operator() (Index_type i)
   {
      Real_type ar, br, cr;
      ar        =      m_cx[i][ 4];
      br        = ar - m_px[i][ 4];
      m_px[i][ 4] = ar;
      cr        = br - m_px[i][ 5];
      m_px[i][ 5] = br;
      ar        = cr - m_px[i][ 6];
      m_px[i][ 6] = cr;
      br        = ar - m_px[i][ 7];
      m_px[i][ 7] = ar;
      cr        = br - m_px[i][ 8];
      m_px[i][ 8] = br;
      ar        = cr - m_px[i][ 9];
      m_px[i][ 9] = cr;
      br        = ar - m_px[i][10];
      m_px[i][10] = ar;
      cr        = br - m_px[i][11];
      m_px[i][11] = br;
      m_px[i][13] = cr - m_px[i][12];
      m_px[i][12] = cr;
   }

   Real_ptr* m_px;
   Real_ptr* m_cx;
};
#endif

#if defined(COMPILE_FIRST_SUM)
/*
 *******************************************************************
 *   Kernel 11 -- first sum
 *******************************************************************
 *    DO 11 L = 1,Loop
 *        X(1)= Y(1)
 *    DO 11 k = 2,n
 * 11     X(k)= X(k-1) + Y(k)
 */
class FIRST_SUM_Functor
{
public:
   FIRST_SUM_Functor(Real_ptr x, Real_ptr y)
   : m_x(x), m_y(y) { ; }

   void operator() (Index_type k)
   {
      m_x[k] = m_x[k-1] + m_y[k];
   }

   Real_ptr m_x;
   Real_ptr m_y;
};
#endif

#if defined(COMPILE_FIRST_DIFF)
/*
 *******************************************************************
 *   Kernel 12 -- first difference
 *******************************************************************
 *    DO 12 L = 1,Loop
 *    DO 12 k = 1,n
 * 12     X(k)= Y(k+1) - Y(k)
 */
class FIRST_DIFF_Functor
{
public:
   FIRST_DIFF_Functor(Real_ptr x, Real_ptr y)
   : m_x(x), m_y(y) { ; }

   void operator() (Index_type k)
   {
      m_x[k] = m_y[k+1] - m_y[k];
   }

   Real_ptr m_x;
   Real_ptr m_y;
};
#endif

#if defined(COMPILE_PIC_2D)
/*
 *******************************************************************
 *   Kernel 13 -- 2-D PIC (Particle In Cell)
 *******************************************************************
 *    DO  13     L= 1,Loop
 *    DO  13    ip= 1,n
 *              i1= P(1,ip)
 *              j1= P(2,ip)
 *              i1=        1 + MOD2N(i1,64)
 *              j1=        1 + MOD2N(j1,64)
 *         P(3,ip)= P(3,ip)  + B(i1,j1)
 *         P(4,ip)= P(4,ip)  + C(i1,j1)
 *         P(1,ip)= P(1,ip)  + P(3,ip)
 *         P(2,ip)= P(2,ip)  + P(4,ip)
 *              i2= P(1,ip)
 *              j2= P(2,ip)
 *              i2=            MOD2N(i2,64)
 *              j2=            MOD2N(j2,64)
 *         P(1,ip)= P(1,ip)  + Y(i2+32)
 *         P(2,ip)= P(2,ip)  + Z(j2+32)
 *              i2= i2       + E(i2+32)
 *              j2= j2       + F(j2+32)
 *        H(i2,j2)= H(i2,j2) + 1.0
 * 13 CONTINUE
 */

//
// RDH modified kernel to fix zero-based indexing error.
//

class PIC_2D_Functor
{
public:
   PIC_2D_Functor(Real_ptr* p, Real_ptr* b, Real_ptr* c,
                  Real_ptr y, Real_ptr z,
                  Index_type* e, Index_type* f,
                  Real_ptr* h)
   : m_p(p), m_b(b), m_c(c),
     m_y(y), m_z(z),
     m_e(e), m_f(f),
     m_h(h) { ; }

   void operator() (Index_type ip)
   {
      Index_type i1, j1, i2, j2;
      i1 = (Index_type) m_p[ip][0];
      j1 = (Index_type) m_p[ip][1];
      i1 &= 64-1;
      j1 &= 64-1;
      m_p[ip][2] += m_b[j1][i1];
      m_p[ip][3] += m_c[j1][i1];
      m_p[ip][0] += m_p[ip][2];
      m_p[ip][1] += m_p[ip][3];
      i2 = (Index_type) m_p[ip][0];
      j2 = (Index_type) m_p[ip][1];
      i2 = ( i2 & 64-1 ) ;
      j2 = ( j2 & 64-1 ) ;
      m_p[ip][0] += m_y[i2+32];
      m_p[ip][1] += m_z[j2+32];
      i2 += m_e[i2+32];
      j2 += m_f[j2+32];
      m_h[j2][i2] += 1.0;
   }

   Real_ptr* m_p;
   Real_ptr* m_b;
   Real_ptr* m_c;
   Real_ptr m_y;
   Real_ptr m_z;
   Index_type* m_e;
   Index_type* m_f;
   Real_ptr* m_h;
};

#if defined(COMPILE_OMP_VARIANTS)
class PIC_2D_Functor_OMP
{
public:
   PIC_2D_Functor_OMP(Real_ptr* p, Real_ptr* b, Real_ptr* c,
                      Real_ptr y, Real_ptr z,
                      Index_type* e, Index_type* f,
                      Real_ptr* h)
   : m_p(p), m_b(b), m_c(c),
     m_y(y), m_z(z),
     m_e(e), m_f(f),
     m_h(h) { ; }

   void operator() (Index_type ip)
   {
      Index_type i1, j1, i2, j2;
      i1 = (Index_type) m_p[ip][0];
      j1 = (Index_type) m_p[ip][1];
      i1 &= 64-1;
      j1 &= 64-1;
      m_p[ip][2] += m_b[j1][i1];
      m_p[ip][3] += m_c[j1][i1];
      m_p[ip][0] += m_p[ip][2];
      m_p[ip][1] += m_p[ip][3];
      i2 = (Index_type) m_p[ip][0];
      j2 = (Index_type) m_p[ip][1];
      i2 = ( i2 & 64-1 ) ;
      j2 = ( j2 & 64-1 ) ;
      m_p[ip][0] += m_y[i2+32];
      m_p[ip][1] += m_z[j2+32];
      i2 += m_e[i2+32];
      j2 += m_f[j2+32];
      #pragma omp atomic
      m_h[j2][i2] += 1.0;
   }

   Real_ptr* m_p;
   Real_ptr* m_b;
   Real_ptr* m_c;
   Real_ptr m_y;
   Real_ptr m_z;
   Index_type* m_e;
   Index_type* m_f;
   Real_ptr* m_h;
};
#endif

#endif  // COMPILE_PIC_2D

#if defined(COMPILE_PIC_1D)
/*
 *******************************************************************
 *   Kernel 14 -- 1-D PIC (Particle In Cell)
 *******************************************************************
 *    DO   14   L= 1,Loop
 *    DO   141  k= 1,n
 *          VX(k)= 0.0
 *          XX(k)= 0.0
 *          IX(k)= INT(  GRD(k))
 *          XI(k)= REAL( IX(k))
 *         EX1(k)= EX   ( IX(k))
 *        DEX1(k)= DEX  ( IX(k))
 *41  CONTINUE
 *    DO   142  k= 1,n
 *          VX(k)= VX(k) + EX1(k) + (XX(k) - XI(k))*DEX1(k)
 *          XX(k)= XX(k) + VX(k)  + FLX
 *          IR(k)= XX(k)
 *          RX(k)= XX(k) - IR(k)
 *          IR(k)= MOD2N(  IR(k),2048) + 1
 *          XX(k)= RX(k) + IR(k)
 *42  CONTINUE
 *    DO  14    k= 1,n
 *    RH(IR(k)  )= RH(IR(k)  ) + 1.0 - RX(k)
 *    RH(IR(k)+1)= RH(IR(k)+1) + RX(k)
 *14  CONTINUE
 */

class PIC_1D_FunctorA
{
public:
   PIC_1D_FunctorA(Real_ptr vx, Real_ptr xx, Real_ptr xi,
                   Real_ptr ex, Real_ptr ex1, 
                   Real_ptr dex, Real_ptr dex1,
                   Index_type* ix,
                   Index_type* grd) 
   : m_vx(vx), m_xx(xx), m_xi(xi),
     m_ex(ex), m_ex1(ex1),
     m_dex(ex), m_dex1(dex1),
     m_ix(ix), 
     m_grd(grd) { ; }

   void operator() (Index_type k)
   {
      m_vx[k] = 0.0;
      m_xx[k] = 0.0;
      m_ix[k] = (Index_type) m_grd[k];
      m_xi[k] = (Real_type) m_ix[k];
      m_ex1[k] = m_ex[ m_ix[k] - 1 ];
      m_dex1[k] = m_dex[ m_ix[k] - 1 ];
   }

   Real_ptr m_vx;
   Real_ptr m_xx;
   Real_ptr m_xi;
   Real_ptr m_ex;
   Real_ptr m_ex1;
   Real_ptr m_dex;
   Real_ptr m_dex1;
   Index_type* m_ix;
   Index_type* m_grd;
};

class PIC_1D_FunctorB
{
public:
   PIC_1D_FunctorB(Real_ptr vx, Real_ptr xx, Real_ptr xi,
                   Real_ptr ex1,
                   Real_ptr dex1,
                   Real_ptr rx, 
                   Index_type* ir,
                   const Real_type flx)
   : m_vx(vx), m_xx(xx), m_xi(xi),
     m_ex1(ex1),
     m_dex1(dex1),
     m_rx(rx),
     m_ir(ir),
     m_flx(flx) { ; }

   void operator() (Index_type k)
   {
      m_vx[k] = m_vx[k] + m_ex1[k] + ( m_xx[k] - m_xi[k] )*m_dex1[k];
      m_xx[k] = m_xx[k] + m_vx[k]  + m_flx;
      m_ir[k] = (Index_type) m_xx[k];
      m_rx[k] = m_xx[k] - m_ir[k];
      m_ir[k] = ( m_ir[k] & (2048-1) ) + 1;
      m_xx[k] = m_rx[k] + m_ir[k];
   }

   Real_ptr m_vx;
   Real_ptr m_xx;
   Real_ptr m_xi;
   Real_ptr m_ex1;
   Real_ptr m_dex1;
   Real_ptr m_rx;
   Index_type* m_ir; 
   const Real_type m_flx;
};

class PIC_1D_FunctorC
{
public:
   PIC_1D_FunctorC(Real_ptr rh, Real_ptr rx,
                   Index_type* ir)
   : m_rh(rh), m_rx(rx),
     m_ir(ir) { ; }

   void operator() (Index_type k)
   {
      m_rh[ m_ir[k]-1 ] += 1.0 - m_rx[k];
      m_rh[ m_ir[k]   ] += m_rx[k];
   }

   Real_ptr m_rh;
   Real_ptr m_rx;
   Index_type* m_ir;
};

#endif
   
#if defined(COMPILE_HYDRO_2D)
/*
 *******************************************************************
 *   Kernel 18 - 2-D explicit hydrodynamics fragment
 *******************************************************************
 *       DO 75  L= 1,Loop
 *              T= 0.0037
 *              S= 0.0041
 *             KN= 6
 *             JN= n
 *       DO 70  k= 2,KN
 *       DO 70  j= 2,JN
 *        ZA(j,k)= (ZP(j-1,k+1)+ZQ(j-1,k+1)-ZP(j-1,k)-ZQ(j-1,k))
 *   .            *(ZR(j,k)+ZR(j-1,k))/(ZM(j-1,k)+ZM(j-1,k+1))
 *        ZB(j,k)= (ZP(j-1,k)+ZQ(j-1,k)-ZP(j,k)-ZQ(j,k))
 *   .            *(ZR(j,k)+ZR(j,k-1))/(ZM(j,k)+ZM(j-1,k))
 * 70    CONTINUE
 *       DO 72  k= 2,KN
 *       DO 72  j= 2,JN
 *        ZU(j,k)= ZU(j,k)+S*(ZA(j,k)*(ZZ(j,k)-ZZ(j+1,k))
 *   .                    -ZA(j-1,k) *(ZZ(j,k)-ZZ(j-1,k))
 *   .                    -ZB(j,k)   *(ZZ(j,k)-ZZ(j,k-1))
 *   .                    +ZB(j,k+1) *(ZZ(j,k)-ZZ(j,k+1)))
 *        ZV(j,k)= ZV(j,k)+S*(ZA(j,k)*(ZR(j,k)-ZR(j+1,k))
 *   .                    -ZA(j-1,k) *(ZR(j,k)-ZR(j-1,k))
 *   .                    -ZB(j,k)   *(ZR(j,k)-ZR(j,k-1))
 *   .                    +ZB(j,k+1) *(ZR(j,k)-ZR(j,k+1)))
 * 72    CONTINUE
 *       DO 75  k= 2,KN
 *       DO 75  j= 2,JN
 *        ZR(j,k)= ZR(j,k)+T*ZU(j,k)
 *        ZZ(j,k)= ZZ(j,k)+T*ZV(j,k)
 * 75    CONTINUE
 */
class HYDRO_2D_FunctorA
{
public:
   HYDRO_2D_FunctorA(Real_ptr* za, Real_ptr* zb,
                     Real_ptr* zm, Real_ptr* zp, Real_ptr* zq, Real_ptr* zr,
                     Index_type* k)
   : m_za(za), m_zb(zb),
     m_zm(zm), m_zp(zp), m_zq(zq), m_zr(zr),
     m_k(k) { ; }

   void operator() (Index_type j)
   {
      Index_type tk = *m_k;
      m_za[tk][j] = ( m_zp[tk+1][j-1] +m_zq[tk+1][j-1] -m_zp[tk][j-1] -m_zq[tk][j-1] )*
                 ( m_zr[tk][j] +m_zr[tk][j-1] ) / ( m_zm[tk][j-1] +m_zm[tk+1][j-1]);
      m_zb[tk][j] = ( m_zp[tk][j-1] +m_zq[tk][j-1] -m_zp[tk][j] -m_zq[tk][j] ) *
                 ( m_zr[tk][j] +m_zr[tk-1][j] ) / ( m_zm[tk][j] +m_zm[tk][j-1]);
   }

   Real_ptr* m_za;
   Real_ptr* m_zb;
   Real_ptr* m_zm;
   Real_ptr* m_zp;
   Real_ptr* m_zq;
   Real_ptr* m_zr;
   Index_type* m_k;
};

class HYDRO_2D_FunctorB
{
public:
   HYDRO_2D_FunctorB(Real_ptr* zu, Real_ptr* zv,
                     Real_ptr* za, Real_ptr* zb, Real_ptr* zr, Real_ptr* zz,
                     Real_type s, 
                     Index_type* k)
   : m_zu(zu), m_zv(zv),
     m_za(za), m_zb(zb), m_zr(zr), m_zz(zz),
     m_s(s),
     m_k(k) { ; }

   void operator() (Index_type j)
   {
      Index_type tk = *m_k;
      m_zu[tk][j] += m_s*( m_za[tk][j]   *( m_zz[tk][j] - m_zz[tk][j+1] ) -
                      m_za[tk][j-1] *( m_zz[tk][j] - m_zz[tk][j-1] ) -
                      m_zb[tk][j]   *( m_zz[tk][j] - m_zz[tk-1][j] ) +
                      m_zb[tk+1][j] *( m_zz[tk][j] - m_zz[tk+1][j] ) );
      m_zv[tk][j] += m_s*( m_za[tk][j]   *( m_zr[tk][j] - m_zr[tk][j+1] ) -
                      m_za[tk][j-1] *( m_zr[tk][j] - m_zr[tk][j-1] ) -
                      m_zb[tk][j]   *( m_zr[tk][j] - m_zr[tk-1][j] ) +
                      m_zb[tk+1][j] *( m_zr[tk][j] - m_zr[tk+1][j] ) );
   }

   Real_ptr* m_zu;
   Real_ptr* m_zv;
   Real_ptr* m_za;
   Real_ptr* m_zb;
   Real_ptr* m_zr;
   Real_ptr* m_zz;
   Real_type m_s;
   Index_type* m_k;
};

class HYDRO_2D_FunctorC
{
public:
   HYDRO_2D_FunctorC(Real_ptr* zrout, Real_ptr* zzout,
                     Real_ptr* zr, Real_ptr* zz,
                     Real_ptr* zu, Real_ptr* zv,
                     Real_type t,
                     Index_type* k)
   : m_zrout(zrout), m_zzout(zzout),
     m_zr(zr), m_zz(zz),
     m_zu(zu), m_zv(zv),
     m_t(t),
     m_k(k) { ; }
 
   void operator() (Index_type j)
   {       
      Index_type tk = *m_k;
      m_zrout[tk][j] = m_zr[tk][j] + m_t*m_zu[tk][j];
      m_zzout[tk][j] = m_zz[tk][j] + m_t*m_zv[tk][j];
   }

   Real_ptr* m_zrout;
   Real_ptr* m_zzout;
   Real_ptr* m_zr;
   Real_ptr* m_zz;
   Real_ptr* m_zu;
   Real_ptr* m_zv;
   Real_type m_t;
   Index_type* m_k;
};
#endif         

#if defined(COMPILE_GEN_LIN_RECUR)
/*
 *******************************************************************
 *   Kernel 19 -- general linear recurrence equations
 *******************************************************************
 *               KB5I= 0
 *           DO 194 L= 1,Loop
 *           DO 191 k= 1,n
 *         B5(k+KB5I)= SA(k) +STB5*SB(k)
 *               STB5= B5(k+KB5I) -STB5
 *191        CONTINUE
 *192        DO 193 i= 1,n
 *                  k= n-i+1
 *         B5(k+KB5I)= SA(k) +STB5*SB(k)
 *               STB5= B5(k+KB5I) -STB5
 *193        CONTINUE
 *194 CONTINUE
 */
class GEN_LIN_RECUR_FunctorA
{
public:
   GEN_LIN_RECUR_FunctorA(Real_ptr b5, Real_ptr sa, Real_ptr sb,
                          Real_type* stb5)
   : m_b5(b5), m_sa(sa), m_sb(sb),
     m_stb5(stb5),
     m_kb5i(0) { ; }

   void operator() (Index_type k)
   {
      m_b5[k+m_kb5i] = m_sa[k] + (*m_stb5)*m_sb[k];
      *m_stb5 = m_b5[k+m_kb5i] - (*m_stb5);
   }

   Real_ptr m_b5;
   Real_ptr m_sa;
   Real_ptr m_sb;
   Real_type* m_stb5;
   Index_type m_kb5i;
};

class GEN_LIN_RECUR_FunctorB
{
public:
   GEN_LIN_RECUR_FunctorB(Real_ptr b5, Real_ptr sa, Real_ptr sb,
                          Real_type* stb5,
                          Index_type len)
   : m_b5(b5), m_sa(sa), m_sb(sb),
     m_stb5(stb5), 
     m_len(len),
     m_kb5i(0) { ; }

   void operator() (Index_type i)
   {
      Index_type k = m_len - i ;
      m_b5[k+m_kb5i] = m_sa[k] + (*m_stb5)*m_sb[k];
      *m_stb5 = m_b5[k+m_kb5i] - (*m_stb5);
   }

   Real_ptr m_b5;
   Real_ptr m_sa;
   Real_ptr m_sb;
   Real_type* m_stb5;
   Index_type m_kb5i;
   Index_type m_len;
};
#endif

#if defined(COMPILE_DISC_ORD)
/*
 *******************************************************************
 *   Kernel 20 -- Discrete ordinates transport, cond recurrence on xx
 *******************************************************************
 *    DO 20 L= 1,Loop
 *    DO 20 k= 1,n
 *         DI= Y(k)-G(k)/( XX(k)+DK)
 *         DN= 0.2
 *         IF( DI.NE.0.0) DN= MAX( S,MIN( Z(k)/DI, T))
 *       X(k)= ((W(k)+V(k)*DN)* XX(k)+U(k))/(VX(k)+V(k)*DN)
 *    XX(k+1)= (X(k)- XX(k))*DN+ XX(k)
 * 20 CONTINUE
 */
class DISC_ORD_Functor 
{
public:
   DISC_ORD_Functor(Real_ptr x, Real_ptr y, Real_ptr z,
                    Real_ptr u, Real_ptr v, Real_ptr w, Real_ptr g,
                    Real_ptr xx, Real_ptr vx,
                    const Real_type s, const Real_type t, const Real_type dk)
   : m_x(x), m_y(y), m_z(z),
     m_u(u), m_v(v), m_w(w), m_g(g),
     m_xx(xx), m_vx(vx),
     m_s(s), m_t(t), m_dk(dk) { ; }

   void operator() (Index_type k)
   {
      Real_type di = m_y[k] - m_g[k] / ( m_xx[k] + m_dk );
      Real_type dn = 0.2;
      if ( di ) {
          dn = m_z[k]/di ;
          if ( m_t < dn ) dn = m_t;
          if ( m_s > dn ) dn = m_s;
      }
      m_x[k] = ( ( m_w[k] + m_v[k]*dn )* m_xx[k] + m_u[k] ) / 
                 ( m_vx[k] + m_v[k]*dn );
      m_xx[k+1] = ( m_x[k] - m_xx[k] )* dn + m_xx[k];
   }

   Real_ptr m_x;
   Real_ptr m_y;
   Real_ptr m_z;
   Real_ptr m_u;
   Real_ptr m_v;
   Real_ptr m_w;
   Real_ptr m_g;
   Real_ptr m_xx;
   Real_ptr m_vx;
   const Real_type m_s;
   const Real_type m_t;
   const Real_type m_dk;
};
#endif

#if defined(COMPILE_MAT_X_MAT)
/*
 *******************************************************************
 *   Kernel 21 -- matrix*matrix product
 *******************************************************************
 *    DO 21 L= 1,Loop
 *    DO 21 k= 1,25
 *    DO 21 i= 1,25
 *    DO 21 j= 1,n
 *    PX(i,j)= PX(i,j) +VY(i,k) * CX(k,j)
 * 21 CONTINUE
 */
class MAT_X_MAT_Functor
{
public:
   MAT_X_MAT_Functor(Real_ptr* px, Real_ptr* cx, Real_ptr* vy,
                     Index_type* k, Index_type* i)
   : m_px(px), m_cx(cx), m_vy(vy),
     m_k(k), m_i(i) { ; }

   void operator() (Index_type j)
   {
      m_px[j][*m_i] += m_vy[*m_k][*m_i] * m_cx[j][*m_k];
   }

   Real_ptr* m_px;
   Real_ptr* m_cx;
   Real_ptr* m_vy;
   Index_type* m_k;
   Index_type* m_i;
};
#endif

#if defined(COMPILE_PLANCKIAN)
/*
 *******************************************************************
 *   Kernel 22 -- Planckian distribution
 *******************************************************************
 *     EXPMAX= 20.0
 *       U(n)= 0.99*EXPMAX*V(n)
 *    DO 22 L= 1,Loop
 *    DO 22 k= 1,n
 *                                          Y(k)= U(k)/V(k)
 *       W(k)= X(k)/( EXP( Y(k)) -1.0)
 * 22 CONTINUE
 */
class PLANCKIAN_Functor
{
public:
   PLANCKIAN_Functor(Real_ptr x, Real_ptr y,
                     Real_ptr u, Real_ptr v, Real_ptr w)
   : m_x(x), m_y(y),
     m_u(u), m_v(v), m_w(w) { ; }

   void operator() (Index_type k)
   {
      m_y[k] = m_u[k] / m_v[k];
      m_w[k] = m_x[k] / ( exp( m_y[k] ) -1.0 );
   }

   Real_ptr m_x;
   Real_ptr m_y;
   Real_ptr m_u;
   Real_ptr m_v;
   Real_ptr m_w;
};
#endif

#if defined(COMPILE_IMP_HYDRO_2D)
/*
 *******************************************************************
 *   Kernel 23 -- 2-D implicit hydrodynamics fragment
 *******************************************************************
 *    DO 23  L= 1,Loop
 *    DO 23  j= 2,6
 *    DO 23  k= 2,n
 *          QA= ZA(k,j+1)*ZR(k,j) +ZA(k,j-1)*ZB(k,j) +
 *   .          ZA(k+1,j)*ZU(k,j) +ZA(k-1,j)*ZV(k,j) +ZZ(k,j)
 * 23  ZA(k,j)= ZA(k,j) +.175*(QA -ZA(k,j))
 */
class IMP_HYDRO_2D_Functor
{
public:
   IMP_HYDRO_2D_Functor(Real_ptr* za, Real_ptr* zb,
                        Real_ptr* zr, Real_ptr* zu,
                        Real_ptr* zv, Real_ptr* zz,
                        Index_type* j)
   : m_za(za), m_zb(zb),
     m_zr(zr), m_zu(zu),
     m_zv(zv), m_zz(zz),
     m_j(j) { ; }

   void operator() (Index_type k)
   {
      Index_type tj = *m_j;
      Real_type qa = m_za[tj+1][k]*m_zr[tj][k] + m_za[tj-1][k]*m_zb[tj][k] +
                     m_za[tj][k+1]*m_zu[tj][k] + m_za[tj][k-1]*m_zv[tj][k] +
                     m_zz[tj][k];
      m_za[tj][k] += 0.175*( qa - m_za[tj][k] );
   }

   Real_ptr* m_za;
   Real_ptr* m_zb;
   Real_ptr* m_zr;
   Real_ptr* m_zu;
   Real_ptr* m_zv;
   Real_ptr* m_zz;
   Index_type* m_j;
};
#endif

#if defined(COMPILE_FIND_FIRST_MIN)
/*
 *******************************************************************
 *   Kernel 24 -- find location of first minimum in array
 *******************************************************************
 *     X( n/2)= -1.0E+10
 *    DO 24  L= 1,Loop
 *           m= 1
 *    DO 24  k= 2,n
 *          IF( X(k).LT.X(m))  m= k
 * 24 CONTINUE
 */
class FIND_FIRST_MIN_Functor
{
public:
   FIND_FIRST_MIN_Functor(Real_ptr x, Index_type* m)
   : m_x(x),
     m_m(m) { ; }

   void operator() (Index_type k)
   {
      if ( m_x[k] < m_x[*m_m] ) *m_m = k;
   }

   Real_ptr m_x;
   Index_type* m_m;
};
#endif


#endif  // closing endif for header file include guard
