#ifndef HMM_CONFIG_HH
#define HMM_CONFIG_HH

// ! do we have a linear elliptic problem?
// if yes, #define LINEAR_PROBLEM
// #define LINEAR_PROBLEM

// ! TFR-HMM or simple HMM?
// #define TFR

// ! is an exact solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro EXACTSOLUTION_AVAILABLE

// ! is the homogenized solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro HOMOGENIZEDSOL_AVAILABLE
// (if HOMOGENIZEDSOL_AVAILABLE == true, it means that it can be computed. It still needs to be determined by using a
// homogenizer )
// #define HOMOGENIZEDSOL_AVAILABLE

// ! Do we solve the cell problems ad hoc or in a pre-process?
// the second possibility requires a data file where the solutions are saved (file becomes large)
// #define AD_HOC_COMPUTATION

// ! Do we have/want a fine-scale reference solution?
// #define FINE_SCALE_REFERENCE
#ifdef FINE_SCALE_REFERENCE
// load the precomputed fine scale reference from a file
 #define FSR_LOAD
 #ifndef FSR_LOAD
// compute the fine scale reference (on the fly)
  #define FSR_COMPUTE
  #ifdef FSR_COMPUTE
// Do we write the discrete fine-scale solution to a file? (for later usage)
   #define WRITE_FINESCALE_SOL_TO_FILE
  #endif       // FSR_COMPUTE
 #endif    // FSR_LOAD
#endif // FINE_SCALE_REFERENCE

// ! Do we have a HMM reference solution? (precomputed detailed HMM simulation)
// we might use a detailed HMM computation as a reference! (if it is available)
// #define HMM_REFERENCE

// ! Do we write the discrete HMM solution to a file? (for later usage)
#define WRITE_HMM_SOL_TO_FILE

// ! Do we want to use error estimation (a-posteriori estimate and adaptivity)?
// Not possible for ad-hoc computations! (in this case, error estimation is far too expensive)
#ifndef AD_HOC_COMPUTATION
//
// #define ERRORESTIMATION
// only possible if we use error estimation:
 #ifdef ERRORESTIMATION
// Do you want to allow adaptive mesh refinement?
// #define ADAPTIVE
 #endif // ifdef ERRORESTIMATION
#endif // ifndef AD_HOC_COMPUTATION

// ! Do we want to add a stochastic perturbation on the data?
// #define STOCHASTIC_PERTURBATION
#ifdef STOCHASTIC_PERTURBATION
// size of variance:
 #define VARIANCE 0.01
// ! Do we want to force the algorithm to come to an end?
// (was auch immer der Grund war, dass das Programm zuvor endlos lange weiter gelaufen ist. z.B. Tolerenzen nicht
// erreicht etc.)
 #define FORCE
#endif // ifdef STOCHASTIC_PERTURBATION

// ! if a computation was broken (after a certain HMM Newton step), we might want to resume to this computation,
// ! loading the solution of the last step that was succesfully carried out (it has to be saved somewhere!)
// (this only works for non-adaptive computations!)
// #define RESUME_TO_BROKEN_COMPUTATION

#ifdef RESUME_TO_BROKEN_COMPUTATION
// last HMM Newton step that was succesfully carried out, saving the iterate afterwards
 #define HMM_NEWTON_ITERATION_STEP 2
#else // ifdef RESUME_TO_BROKEN_COMPUTATION
      // default: we need a full computation. start with step 1:
 #define HMM_NEWTON_ITERATION_STEP 0
#endif // ifdef RESUME_TO_BROKEN_COMPUTATION

#if 1
// hmfemmain:
// !---------
 #include <dune/fem/solver/oemsolver/oemsolver.hh>

 #include <dune/fem/operator/discreteoperatorimp.hh>
 #include <dune/fem/misc/l2error.hh>
 #include <dune/fem/misc/l2norm.hh>
// #include <dune/fem/io/visual/grape/datadisp/errordisplay.hh>
// !-----------
#endif // if 1
#if 0
// poisson
// !--------

 #include <dune/fem/operator/matrix/spmatrix.hh>
 #include <dune/fem/operator/2order/lagrangematrixsetup.hh>

 #include <dune/fem/misc/l2norm.hh>
 #include <dune/fem/misc/h1norm.hh>
 #include <dune/fem/misc/mpimanager.hh>
 #include <dune/fem/io/parameter.hh>
 #include <dune/fem/io/visual/grape/datadisp/errordisplay.hh>
// !---------
#endif // if 0

#include <dune/multiscale/problems/elliptic_problems/model_problem_1/problem_specification.hh>

#endif // HMM_CONFIG_HH
