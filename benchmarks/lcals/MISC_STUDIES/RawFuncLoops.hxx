//
//  THIS IS NOT OPEN SOURCE OR PUBLIC DOMAIN SOFTWARE
//
// See README-LCALS_license.txt for access and distribution restrictions
//

//
// Header file with prototypes of routines that implement "A" subset raw loops
//

#ifndef RawFuncLoops_HXX
#define RawFuncLoops_HXX

#include "LCALSSuite.hxx"

#include "SubsetDataA.hxx"

//
//  "A" subset
//
void doPRESSURE_CALC(LoopTimer& ltimer, int num_samples, Index_type len,
                     Real_ptr compression, Real_ptr bvc, Real_ptr p_new,
                     Real_ptr e_old, Real_ptr vnewc);

void doENERGY_CALC(LoopTimer& ltimer, int num_samples, Index_type len,
                   Real_ptr e_new, Real_ptr e_old, Real_ptr delvc,
                   Real_ptr p_new, Real_ptr p_old, Real_ptr q_new,
                   Real_ptr q_old, Real_ptr work, Real_ptr compHalfStep,
                   Real_ptr pHalfStep, Real_ptr bvc, Real_ptr pbvc,
                   Real_ptr ql_old, Real_ptr qq_old, Real_ptr vnewc);

void doVOL3D_CALC(LoopTimer& ltimer, int num_samples, ADomain& domain,
                  Real_ptr x, Real_ptr y, Real_ptr z, Real_ptr vol);

void doDEL_DOT_VEC_2D(LoopTimer& ltimer, int num_samples, ADomain& domain,
                      Real_ptr x, Real_ptr y, Real_ptr xdot, Real_ptr ydot,
                      Real_ptr div);

void doCOUPLE(LoopTimer& ltimer, int num_samples, ADomain& domain,
              Complex_ptr t0, Complex_ptr t1, Complex_ptr t2,
              Complex_ptr denac, Complex_ptr denlw);

Index_type doFIR(LoopTimer& ltimer, int num_samples, Index_type len,
                Real_ptr out, Real_ptr in);


//
//  "B" subset
//
void doINIT3(LoopTimer& ltimer, int num_samples, Index_type len,
             Real_ptr out1, Real_ptr out2, Real_ptr out3,
             Real_ptr in1, Real_ptr in2);

void doMULADDSUB(LoopTimer& ltimer, int num_samples, Index_type len,
                 Real_ptr out1, Real_ptr out2, Real_ptr out3,
                 Real_ptr in1, Real_ptr in2);

void doIF_QUAD(LoopTimer& ltimer, int num_samples, Index_type len,
               Real_ptr a, Real_ptr b, Real_ptr c,
               Real_ptr x1, Real_ptr x2);


//
//  "C" subset
//
void doHYDRO_1D(LoopTimer& ltimer, int num_samples, Index_type len,
                Real_ptr x, Real_ptr y, Real_ptr z);
         

#endif  // closing endif for header file include guard
