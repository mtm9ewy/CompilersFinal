//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Header file with functors for LCALS "B" subset
//

#ifndef LCALSFunctorKernelsB_HXX
#define LCALSFunctorKernelsB_HXX

#include "LCALSParams.hxx"


#if defined(COMPILE_INIT3)
class INIT3_Functor
{
public:
   INIT3_Functor(Real_ptr out1, Real_ptr out2, Real_ptr out3,
                 Real_ptr in1, Real_ptr in2) 
   : m_out1(out1), m_out2(out2), m_out3(out3),
     m_in1(in1), m_in2(in2) { ; }

   void operator() (Index_type i)
   {
      m_out1[i] = m_out2[i] = m_out3[i] = - m_in1[i] - m_in2[i];
   }

   Real_ptr m_out1;
   Real_ptr m_out2;
   Real_ptr m_out3;
   Real_ptr m_in1;
   Real_ptr m_in2;
};
#endif

#if defined(COMPILE_MULADDSUB)
class MULADDSUB_Functor
{
public:
   MULADDSUB_Functor(Real_ptr out1, Real_ptr out2, Real_ptr out3,
                     Real_ptr in1, Real_ptr in2) 
   : m_out1(out1), m_out2(out2), m_out3(out3),
     m_in1(in1), m_in2(in2) { ; }

   void operator() (Index_type i)
   {
      m_out1[i] = m_in1[i] * m_in2[i] ;
      m_out2[i] = m_in1[i] + m_in2[i] ;
      m_out3[i] = m_in1[i] - m_in2[i] ;
   }

   Real_ptr m_out1;
   Real_ptr m_out2;
   Real_ptr m_out3;
   Real_ptr m_in1;
   Real_ptr m_in2;
};
#endif

#if defined(COMPILE_IF_QUAD)
class IF_QUAD_Functor
{
public:
   IF_QUAD_Functor(Real_ptr a, Real_ptr b, Real_ptr c,
                   Real_ptr x1, Real_ptr x2)
   : m_a(a), m_b(b), m_c(c), m_x1(x1), m_x2(x2) { ; }

   void operator() (Index_type i)
   {
      Real_type s = m_b[i]*m_b[i] - 4.0*m_a[i]*m_c[i];
      if ( s >= 0 ) {
         s = sqrt(s);
         m_x2[i] = (-m_b[i]+s)/(2.0*m_a[i]);
         m_x1[i] = (-m_b[i]-s)/(2.0*m_a[i]);
      } else {
         m_x2[i] = 0.0;
         m_x1[i] = 0.0;
      }
   }

   Real_ptr m_a;
   Real_ptr m_b;
   Real_ptr m_c;
   Real_ptr m_x1;
   Real_ptr m_x2;
};
#endif

#if defined(COMPILE_TRAP_INT)
class TRAP_INT_Functor
{
public:
   TRAP_INT_Functor(Real_type* sumx, 
                    Real_type x0, Real_type h,
                    Real_type y, Real_type xp, Real_type yp) 
   : m_sumx(sumx), m_x0(x0), m_h(h), m_y(y), m_xp(xp), m_yp(yp) { ; }

   void operator() (Index_type i)
   {
      Real_type x = m_x0 + i*m_h;
      *m_sumx += trap_int_func(x, m_y, m_xp, m_yp);
   }

   Real_type* m_sumx;
   const Real_type m_x0;
   const Real_type m_h;
   const Real_type m_y;
   const Real_type m_xp;
   const Real_type m_yp;
};

#if defined(COMPILE_OMP_VARIANTS)
class TRAP_INT_Functor_OMP
{
public:
   TRAP_INT_Functor_OMP(Real_type* sumx,
                        Real_type x0, Real_type h,
                        Real_type y, Real_type xp, Real_type yp)
   : m_sumx(sumx), m_x0(x0), m_h(h), m_y(y), m_xp(xp), m_yp(yp) { ; }

   void operator() (Index_type i)
   {
      Real_type x = m_x0 + i*m_h;
      #pragma omp atomic
      *m_sumx += trap_int_func(x, m_y, m_xp, m_yp);
   }

   Real_type* m_sumx;
   const Real_type m_x0;
   const Real_type m_h;
   const Real_type m_y;
   const Real_type m_xp;
   const Real_type m_yp;
};
#endif

#endif // COMPILE_TRAP_INT


#endif  // closing endif for header file include guard
