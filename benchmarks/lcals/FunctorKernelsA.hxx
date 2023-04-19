//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Header file with functors for LCALS "A" subset
//

#ifndef LCALSFunctorKernelsA_HXX
#define LCALSFunctorKernelsA_HXX

#include "LCALSParams.hxx"

#include<cmath>


#if defined(COMPILE_PRESSURE_CALC)
class PRESSURE_CALC_FunctorA
{
public:
   PRESSURE_CALC_FunctorA(Real_ptr bvc, Real_ptr compression,
                          Real_type cls)
   : m_bvc(bvc), m_compression(compression),
     m_cls(cls) { ; }

   void operator() (Index_type i)
   {
      m_bvc[i] = m_cls * (m_compression[i] + 1.0);
   }

   Real_ptr m_bvc;
   Real_ptr m_compression;
   const Real_type m_cls;
};

class PRESSURE_CALC_FunctorB
{
public:
   PRESSURE_CALC_FunctorB(Real_ptr p_new, 
                          Real_ptr bvc, Real_ptr e_old, Real_ptr vnewc,
                          Real_type p_cut, Real_type eosvmax, Real_type pmin)
   : m_p_new(p_new), m_bvc(bvc), m_e_old(e_old), m_vnewc(vnewc),
     m_p_cut(p_cut), m_eosvmax(eosvmax), m_pmin(pmin) { ; }

   void operator() (Index_type i)
   {
      m_p_new[i] = m_bvc[i] * m_e_old[i] ;

      if ( fabs(m_p_new[i]) <  m_p_cut )  m_p_new[i] = 0.0 ;

      if ( m_vnewc[i] >= m_eosvmax )  m_p_new[i] = 0.0 ;

      if ( m_p_new[i]  <  m_pmin )  m_p_new[i] = m_pmin ;
   }

   Real_ptr m_p_new;
   Real_ptr m_bvc;
   Real_ptr m_e_old;
   Real_ptr m_vnewc;
   const Real_type m_p_cut;
   const Real_type m_eosvmax;
   const Real_type m_pmin;
};
#endif

#if defined(COMPILE_ENERGY_CALC)
class ENERGY_CALC_FunctorA
{
public:
   ENERGY_CALC_FunctorA(Real_ptr e_new, 
                        Real_ptr e_old, Real_ptr delvc, Real_ptr p_old,
                        Real_ptr q_old, Real_ptr work) 
   : m_e_new(e_new), 
     m_e_old(e_old), m_delvc(delvc), m_p_old(p_old), 
     m_q_old(q_old), m_work(work) { ; }

   void operator() (Index_type i)
   {
      m_e_new[i] = m_e_old[i] - 0.5 * m_delvc[i] *
                            (m_p_old[i] + m_q_old[i]) + 0.5 * m_work[i];
   }

   Real_ptr m_e_new;
   Real_ptr m_e_old;
   Real_ptr m_delvc;
   Real_ptr m_p_old;
   Real_ptr m_q_old;
   Real_ptr m_work;
};

class ENERGY_CALC_FunctorB
{
public:
   ENERGY_CALC_FunctorB(Real_ptr q_new, 
                        Real_ptr delvc, 
                        Real_ptr compHalfStep, Real_ptr pHalfStep,
                        Real_ptr e_new, Real_ptr bvc, Real_ptr pbvc,
                        Real_ptr ql_old, Real_ptr qq_old,
                        Real_type rho0)
   : m_q_new(q_new), 
     m_delvc(delvc), m_compHalfStep(compHalfStep), m_pHalfStep(pHalfStep),
     m_e_new(e_new), m_bvc(bvc), m_pbvc(pbvc),
     m_ql_old(ql_old), m_qq_old(qq_old),
     m_rho0(rho0) { ; }

   void operator() (Index_type i)
   {
      if ( m_delvc[i] > 0.0 ) {
         m_q_new[i] = 0.0 ;
      }
      else {
         Real_type vhalf = 1.0 / (1.0 + m_compHalfStep[i]) ;
         Real_type ssc = ( m_pbvc[i] * m_e_new[i]
            + vhalf * vhalf * m_bvc[i] * m_pHalfStep[i] ) / m_rho0 ;

         if ( ssc <= 0.1111111e-36 ) {
            ssc = 0.3333333e-18 ;
         } else {
            ssc = sqrt(ssc) ;
         }

         m_q_new[i] = (ssc*m_ql_old[i] + m_qq_old[i]) ;
      }
   }

   Real_ptr m_q_new;
   Real_ptr m_delvc;
   Real_ptr m_compHalfStep;
   Real_ptr m_pHalfStep;
   Real_ptr m_e_new;
   Real_ptr m_bvc;
   Real_ptr m_pbvc;
   Real_ptr m_ql_old;
   Real_ptr m_qq_old;
   const Real_type m_rho0;
};

class ENERGY_CALC_FunctorC
{
public:
   ENERGY_CALC_FunctorC(Real_ptr e_new,
                        Real_ptr delvc,
                        Real_ptr pHalfStep,
                        Real_ptr p_old, Real_ptr q_old, Real_ptr q_new)
   : m_e_new(e_new),
     m_delvc(delvc), 
     m_pHalfStep(pHalfStep),
     m_p_old(p_old), m_q_old(q_old), m_q_new(q_new) { ; }

   void operator() (Index_type i)
   {
      m_e_new[i] = m_e_new[i] + 0.5 * m_delvc[i]
                    * ( 3.0*(m_p_old[i] + m_q_old[i])
                         - 4.0*(m_pHalfStep[i] + m_q_new[i])) ;
   }

   Real_ptr m_e_new;
   Real_ptr m_delvc;
   Real_ptr m_pHalfStep;
   Real_ptr m_p_old;
   Real_ptr m_q_old;
   Real_ptr m_q_new;
};

class ENERGY_CALC_FunctorD
{
public:
   ENERGY_CALC_FunctorD(Real_ptr e_new,
                        Real_ptr work,
                        Real_type e_cut, Real_type emin)
   : m_e_new(e_new), m_work(work),
     m_e_cut(e_cut), m_emin(emin) { ; }

   void operator() (Index_type i)
   {
      m_e_new[i] += 0.5 * m_work[i];

      if ( fabs(m_e_new[i]) < m_e_cut ) { m_e_new[i] = 0.0  ; }

      if ( m_e_new[i]  < m_emin ) { m_e_new[i] = m_emin ; }
   }

   Real_ptr m_e_new;
   Real_ptr m_work;
   const Real_type m_e_cut;
   const Real_type m_emin;
};

class ENERGY_CALC_FunctorE
{
public:
   ENERGY_CALC_FunctorE(Real_ptr q_new,
                        Real_ptr delvc,
                        Real_ptr e_new, Real_ptr p_new, Real_ptr vnewc, 
                        Real_ptr bvc, Real_ptr pbvc, 
                        Real_ptr ql_old, Real_ptr qq_old,
                        Real_ptr p_old, Real_ptr q_old, Real_ptr pHalfStep,
                        Real_type rho0, Real_type e_cut, Real_type emin)
   : m_q_new(q_new),
     m_delvc(delvc), 
     m_e_new(e_new), m_p_new(p_new), m_vnewc(vnewc),
     m_bvc(bvc), m_pbvc(pbvc), 
     m_ql_old(ql_old), m_qq_old(qq_old),
     m_p_old(p_old), m_q_old(q_old), m_pHalfStep(pHalfStep),
     m_rho0(rho0), m_e_cut(e_cut), m_emin(emin) { ; }

   void operator() (Index_type i)
   {
      Real_type q_tilde ;

      if (m_delvc[i] > 0.0) {
         q_tilde = 0. ;
      }
      else {
         Real_type ssc = ( m_pbvc[i] * m_e_new[i]
             + m_vnewc[i] * m_vnewc[i] * m_bvc[i] * m_p_new[i] ) / m_rho0 ;

         if ( ssc <= 0.1111111e-36 ) {
            ssc = 0.3333333e-18 ;
         } else {
            ssc = sqrt(ssc) ;
         }

         q_tilde = (ssc*m_ql_old[i] + m_qq_old[i]) ;
      }

      m_e_new[i] = m_e_new[i] - ( 7.0*(m_p_old[i] + m_q_old[i])
                             - 8.0*(m_pHalfStep[i] + m_q_new[i])
                             + (m_p_new[i] + q_tilde)) * m_delvc[i] / 6.0 ;


      if ( fabs(m_e_new[i]) < m_e_cut ) {
         m_e_new[i] = 0.0  ;
      }
      if ( m_e_new[i]  < m_emin ) {
         m_e_new[i] = m_emin ;
      }
   }

   Real_ptr m_q_new;
   Real_ptr m_delvc;
   Real_ptr m_e_new;
   Real_ptr m_p_new;
   Real_ptr m_vnewc;
   Real_ptr m_bvc;
   Real_ptr m_pbvc;
   Real_ptr m_ql_old;
   Real_ptr m_qq_old;
   Real_ptr m_p_old;
   Real_ptr m_q_old;
   Real_ptr m_pHalfStep;
   const Real_type m_rho0;
   const Real_type m_e_cut;
   const Real_type m_emin;
};

class ENERGY_CALC_FunctorF
{
public:
   ENERGY_CALC_FunctorF(Real_ptr q_new,
                        Real_ptr delvc,
                        Real_ptr e_new, Real_ptr p_new, Real_ptr vnewc,
                        Real_ptr bvc, Real_ptr pbvc,
                        Real_ptr ql_old, Real_ptr qq_old,
                        Real_type rho0, Real_type q_cut)
   : m_q_new(q_new),
     m_delvc(delvc),
     m_e_new(e_new), m_p_new(p_new), m_vnewc(vnewc),
     m_bvc(bvc), m_pbvc(pbvc),
     m_ql_old(ql_old), m_qq_old(qq_old),
     m_rho0(rho0), m_q_cut(q_cut) { ; }

   void operator() (Index_type i)
   {
      if ( m_delvc[i] <= 0.0 ) {
         Real_type ssc = ( m_pbvc[i] * m_e_new[i]
                 + m_vnewc[i] * m_vnewc[i] * m_bvc[i] * m_p_new[i] ) / m_rho0 ;

         if ( ssc <= 0.1111111e-36 ) {
            ssc = 0.3333333e-18 ;
         } else {
            ssc = sqrt(ssc) ;
         }

         m_q_new[i] = (ssc*m_ql_old[i] + m_qq_old[i]) ;

         if (fabs(m_q_new[i]) < m_q_cut) m_q_new[i] = 0.0 ;
      }
   }

   Real_ptr m_q_new;
   Real_ptr m_delvc;
   Real_ptr m_e_new;
   Real_ptr m_p_new;
   Real_ptr m_vnewc;
   Real_ptr m_bvc;
   Real_ptr m_pbvc;
   Real_ptr m_ql_old;
   Real_ptr m_qq_old;
   const Real_type m_rho0;
   const Real_type m_q_cut;
};
#endif

#if defined(COMPILE_VOL3D_CALC)
class VOL3D_CALC_Functor
{
public:
   VOL3D_CALC_Functor(Real_ptr vol,
                      Real_ptr x0, Real_ptr x1, Real_ptr x2, 
                      Real_ptr x3, Real_ptr x4, Real_ptr x5, 
                      Real_ptr x6, Real_ptr x7,
                      Real_ptr y0, Real_ptr y1, Real_ptr y2, 
                      Real_ptr y3, Real_ptr y4, Real_ptr y5, 
                      Real_ptr y6, Real_ptr y7,
                      Real_ptr z0, Real_ptr z1, Real_ptr z2, 
                      Real_ptr z3, Real_ptr z4, Real_ptr z5, 
                      Real_ptr z6, Real_ptr z7,
                      Real_type vnormq)
   : m_vol(vol),
     m_x0(x0), m_x1(x1), m_x2(x2),
     m_x3(x3), m_x4(x4), m_x5(x5),
     m_x6(x6), m_x7(x7),
     m_y0(y0), m_y1(y1), m_y2(y2),
     m_y3(y3), m_y4(y4), m_y5(y5),
     m_y6(y6), m_y7(y7),
     m_z0(z0), m_z1(z1), m_z2(z2),
     m_z3(z3), m_z4(z4), m_z5(z5),
     m_z6(z6), m_z7(z7),
     m_vnormq(vnormq) { ; }

   void operator() (Index_type i)
   {
      Real_type x71 = m_x7[i] - m_x1[i] ;
      Real_type x72 = m_x7[i] - m_x2[i] ;
      Real_type x74 = m_x7[i] - m_x4[i] ;
      Real_type x30 = m_x3[i] - m_x0[i] ;
      Real_type x50 = m_x5[i] - m_x0[i] ;
      Real_type x60 = m_x6[i] - m_x0[i] ;

      Real_type y71 = m_y7[i] - m_y1[i] ;
      Real_type y72 = m_y7[i] - m_y2[i] ;
      Real_type y74 = m_y7[i] - m_y4[i] ;
      Real_type y30 = m_y3[i] - m_y0[i] ;
      Real_type y50 = m_y5[i] - m_y0[i] ;
      Real_type y60 = m_y6[i] - m_y0[i] ;

      Real_type z71 = m_z7[i] - m_z1[i] ;
      Real_type z72 = m_z7[i] - m_z2[i] ;
      Real_type z74 = m_z7[i] - m_z4[i] ;
      Real_type z30 = m_z3[i] - m_z0[i] ;
      Real_type z50 = m_z5[i] - m_z0[i] ;
      Real_type z60 = m_z6[i] - m_z0[i] ;

      Real_type xps = x71 + x60 ;
      Real_type yps = y71 + y60 ;
      Real_type zps = z71 + z60 ;

      Real_type cyz = y72 * z30 - z72 * y30 ;
      Real_type czx = z72 * x30 - x72 * z30 ;
      Real_type cxy = x72 * y30 - y72 * x30 ;
      m_vol[i] = xps * cyz + yps * czx + zps * cxy ;

      xps = x72 + x50 ;
      yps = y72 + y50 ;
      zps = z72 + z50 ;

      cyz = y74 * z60 - z74 * y60 ;
      czx = z74 * x60 - x74 * z60 ;
      cxy = x74 * y60 - y74 * x60 ;
      m_vol[i] += xps * cyz + yps * czx + zps * cxy ;

      xps = x74 + x30 ;
      yps = y74 + y30 ;
      zps = z74 + z30 ;

      cyz = y71 * z50 - z71 * y50 ;
      czx = z71 * x50 - x71 * z50 ;
      cxy = x71 * y50 - y71 * x50 ;
      m_vol[i] += xps * cyz + yps * czx + zps * cxy ;

      m_vol[i] *= m_vnormq ;
   }

   Real_ptr m_vol;
   Real_ptr m_x0; Real_ptr m_x1; Real_ptr m_x2;
   Real_ptr m_x3; Real_ptr m_x4; Real_ptr m_x5;
   Real_ptr m_x6; Real_ptr m_x7;
   Real_ptr m_y0; Real_ptr m_y1; Real_ptr m_y2;
   Real_ptr m_y3; Real_ptr m_y4; Real_ptr m_y5;
   Real_ptr m_y6; Real_ptr m_y7;
   Real_ptr m_z0; Real_ptr m_z1; Real_ptr m_z2;
   Real_ptr m_z3; Real_ptr m_z4; Real_ptr m_z5;
   Real_ptr m_z6; Real_ptr m_z7;
   const Real_type m_vnormq; 
};
#endif


#if defined(COMPILE_DEL_DOT_VEC_2D)
class DEL_DOT_VEC_2D_Functor
{
public:
   DEL_DOT_VEC_2D_Functor(Index_type* zones,
                    Real_ptr div,
                    Real_ptr x1, Real_ptr x2, Real_ptr x3, Real_ptr x4,
                    Real_ptr y1, Real_ptr y2, Real_ptr y3, Real_ptr y4,
                    Real_ptr fx1, Real_ptr fx2, Real_ptr fx3, Real_ptr fx4,
                    Real_ptr fy1, Real_ptr fy2, Real_ptr fy3, Real_ptr fy4,
                    const Real_type ptiny)
   : m_zones(zones),
     m_div(div),
     m_x1(x1), m_x2(x2), m_x3(x3), m_x4(x4),
     m_y1(y1), m_y2(y2), m_y3(y3), m_y4(y4),
     m_fx1(fx1), m_fx2(fx2), m_fx3(fx3), m_fx4(fx4),
     m_fy1(fy1), m_fy2(fy2), m_fy3(fy3), m_fy4(fy4),
     m_ptiny(ptiny),
     m_half(0.5) { ; }

   void operator() (Index_type ii)
   {
      Index_type i  = m_zones[ii] ;

      Real_type xi  = m_half * ( m_x1[i]  + m_x2[i]  - m_x3[i]  - m_x4[i]  ) ;
      Real_type xj  = m_half * ( m_x2[i]  + m_x3[i]  - m_x4[i]  - m_x1[i]  ) ;

      Real_type yi  = m_half * ( m_y1[i]  + m_y2[i]  - m_y3[i]  - m_y4[i]  ) ;
      Real_type yj  = m_half * ( m_y2[i]  + m_y3[i]  - m_y4[i]  - m_y1[i]  ) ;

      Real_type fxi = m_half * ( m_fx1[i] + m_fx2[i] - m_fx3[i] - m_fx4[i] ) ;
      Real_type fxj = m_half * ( m_fx2[i] + m_fx3[i] - m_fx4[i] - m_fx1[i] ) ;

      Real_type fyi = m_half * ( m_fy1[i] + m_fy2[i] - m_fy3[i] - m_fy4[i] ) ;
      Real_type fyj = m_half * ( m_fy2[i] + m_fy3[i] - m_fy4[i] - m_fy1[i] ) ;

      Real_type rarea  = 1.0 / ( xi * yj - xj * yi + m_ptiny ) ;

      Real_type dfxdx  = rarea * ( fxi * yj - fxj * yi ) ;

      Real_type dfydy  = rarea * ( fyj * xi - fyi * xj ) ;

      Real_type affine = ( m_fy1[i] + m_fy2[i] + m_fy3[i] + m_fy4[i] ) /
                         ( m_y1[i]  + m_y2[i]  + m_y3[i]  + m_y4[i]  ) ;

      m_div[i] = dfxdx + dfydy + affine ;
   } 

   Index_type* m_zones;
   Real_ptr m_div;
   Real_ptr m_x1; Real_ptr m_x2; Real_ptr m_x3; Real_ptr m_x4;
   Real_ptr m_y1; Real_ptr m_y2; Real_ptr m_y3; Real_ptr m_y4;
   Real_ptr m_fx1; Real_ptr m_fx2; Real_ptr m_fx3; Real_ptr m_fx4;
   Real_ptr m_fy1; Real_ptr m_fy2; Real_ptr m_fy3; Real_ptr m_fy4;
   const Real_type m_ptiny; 
   const Real_type m_half; 
};
#endif


#if defined(COMPILE_COUPLE)
class COUPLE_Functor
{
public:
   COUPLE_Functor(Index_type imin, Index_type imax,
                  Index_type jmin, Index_type jmax,
                  Complex_ptr t0, Complex_ptr t1, Complex_ptr t2,
                  Complex_ptr denac, Complex_ptr denlw,
                  const Real_type dt,
                  const Real_type c10,
                  const Real_type fratio,
                  const Real_type r_fratio,
                  const Real_type c20,
                  const Complex_type ireal)
   : m_imin(imin), m_imax(imax),
     m_jmin(jmin), m_jmax(jmax),
     m_t0(t0), m_t1(t1), m_t2(t2),
     m_denac(denac), m_denlw(denlw),
     m_dt(dt),
     m_c10(c10),
     m_fratio(fratio),
     m_r_fratio(r_fratio),
     m_c20(c20),
     m_ireal(ireal) { ; }

   void operator() (Index_type k)
   {
      for (Index_type j = m_jmin; j < m_jmax; j++) {

         Index_type it0=    ((k)*(m_jmax+1) + (j))*(m_imax+1) ;
         Index_type idenac= ((k)*(m_jmax+2) + (j))*(m_imax+2) ;

         for (Index_type i = m_imin; i < m_imax; i++) {

            Complex_type c1 = m_c10 * m_denac[m_idenac+i];
            Complex_type c2 = m_c20 * m_denlw[m_it0+i];

            /* promote to doubles to avoid possible divide by zero
               errors later on. */
            Real_type c1re = real(c1);  Real_type c1im = imag(c1);
            Real_type c2re = real(c2);  Real_type c2im = imag(c2);

            /* compute lamda = sqrt(|c1|^2 + |c2|^2) using doubles
               to avoid underflow. */
            Real_type zlam = c1re*c1re + c1im*c1im +
                             c2re*c2re + c2im*c2im + 1.0e-34;
            zlam = sqrt(zlam);
            Real_type snlamt = sin(zlam * m_dt * 0.5);
            Real_type cslamt = cos(zlam * m_dt * 0.5);

            Complex_type a0t = m_t0[m_it0+i];
            Complex_type a1t = m_t1[m_it0+i];
            Complex_type a2t = m_t2[m_it0+i] * m_fratio;

            Real_type r_zlam= 1.0/zlam;
            c1 *= r_zlam;
            c2 *= r_zlam;
            Real_type zac1 = zabs2(c1);
            Real_type zac2 = zabs2(c2);

            /* compute new A0 */
            Complex_type z3 = ( c1 * a1t + c2 * a2t ) * snlamt ;
            m_t0[m_it0+i] = a0t * cslamt -  m_ireal * z3;

            /* compute new A1  */
            Real_type r = zac1 * cslamt + zac2;
            Complex_type z5 = c2 * a2t;
            Complex_type z4 = conj(c1) * z5 * (cslamt-1);
            z3 = conj(c1) * a0t * snlamt;
            m_t1[m_it0+i] = a1t * r + z4 - m_ireal * z3;

            /* compute new A2  */
            r = zac1 + zac2 * cslamt;
            z5 = c1 * a1t;  
            z4 = conj(c2) * z5 * (cslamt-1);
            z3 = conj(c2) * a0t * snlamt;
            m_t2[m_it0+i] = ( a2t * r + z4 - m_ireal * z3 ) * m_r_fratio;

         }  // i loop

      }  // j loop

   }

   Index_type m_imin; Index_type m_imax;
   Index_type m_jmin; Index_type m_jmax;

   Index_type m_it0;
   Index_type m_idenac;
   Complex_ptr m_t0; Complex_ptr m_t1; Complex_ptr m_t2;
   Complex_ptr m_denac; Complex_ptr m_denlw;

   const Real_type m_dt;
   const Real_type m_c10;
   const Real_type m_fratio;
   const Real_type m_r_fratio;
   const Real_type m_c20;
   const Complex_type m_ireal;
};
#endif


#if defined(COMPILE_FIR)
class FIR_Functor
{
public:
   FIR_Functor(Real_ptr in,
               Real_ptr out,
               Real_ptr coeff,
               const Index_type coefflen)  
   : m_in(in), m_out(out), m_coeff(coeff), m_coefflen(coefflen) { ; }

   void operator() (Index_type i)
   {
      Real_type sum = 0.0;
      for (Index_type j = 0; j < m_coefflen; ++j ) {
         sum += m_coeff[j]*m_in[i+j];
      }
      m_out[i] = sum; 
   }

   Real_ptr m_in;
   Real_ptr m_out;
   Real_ptr m_coeff;
   const Index_type m_coefflen;
};
#endif

#endif  // closing endif for header file include guard
