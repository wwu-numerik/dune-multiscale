#ifndef HMM_CONFIG_HH
#define HMM_CONFIG_HH

/****   all this should go into a default parameter file ********/
//! do we have a linear elliptic problem?
// if yes, #define LINEAR_PROBLEM
// #define LINEAR_PROBLEM

//! TFR-HMM or simple HMM?
// #define TFR


//! is the homogenized solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro HOMOGENIZEDSOL_AVAILABLE
// (if HOMOGENIZEDSOL_AVAILABLE == true, it means that it can be computed. It still needs to be determined by using a
// homogenizer )
// #define HOMOGENIZEDSOL_AVAILABLE

//! Do we solve the cell problems ad hoc or in a pre-process?
// the second possibility requires a data file where the solutions are saved (file becomes large)
// #define AD_HOC_COMPUTATION


//! Do we have a HMM reference solution? (precomputed detailed HMM simulation)
// we might use a detailed HMM computation as a reference! (if it is available)
// #define HMM_REFERENCE

//! Do we write the discrete HMM solution to a file? (for later usage)
#define WRITE_HMM_SOL_TO_FILE


//! if a computation was broken (after a certain HMM Newton step), we might want to resume to this computation,
//! loading the solution of the last step that was succesfully carried out (it has to be saved somewhere!)
// (this only works for non-adaptive computations!)
// #define RESUME_TO_BROKEN_COMPUTATION

/**** ***********************************************************/

#endif // HMM_CONFIG_HH
