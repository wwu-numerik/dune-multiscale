// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef ALGORITHM_STEP_HH
#define ALGORITHM_STEP_HH

#include <dune/multiscale/hmm/hmm_traits.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

struct HMMResult;

/**
 * \return true if a program should continue with a new newton step
 **/
bool process_hmm_newton_residual(typename HMMTraits::RangeType& relative_newton_error,
                                 typename HMMTraits::DiscreteFunctionType& hmm_solution,
                                 const typename HMMTraits::FEMMatrix& hmm_newton_matrix,
                                 const typename HMMTraits::DiscreteFunctionType& hmm_newton_rhs,
                                 const int hmm_iteration_step,
                                 const int loop_cycle,
                                 const double hmm_tolerance );


//! \TODO docme
HMMResult single_step( typename HMMTraits::GridPartType& gridPart,
        typename HMMTraits::GridPartType& gridPartFine,
        typename HMMTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
        typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
        const typename HMMTraits::DiffusionType& diffusion_op,
        const Dune::RightHandSideAssembler< typename HMMTraits::DiscreteFunctionType >& rhsassembler,
        typename HMMTraits::DiscreteFunctionType& hmm_solution,
        const typename HMMTraits::DiscreteFunctionType& reference_solution,
        const int loop_cycle );

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {

#endif // ALGORITHM_STEP_HH
