// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef ALGORITHM_STEP_HH
#define ALGORITHM_STEP_HH

#include <dune/multiscale/hmm/hmm_traits.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

struct HMMResult;

/**
 * \return true if a program should continue with a new newton step
 **/
bool process_hmm_newton_residual(typename CommonTraits::RangeType& relative_newton_error,
                                 typename CommonTraits::DiscreteFunctionType& hmm_solution,
                                 const typename HMMTraits::FEMMatrix& hmm_newton_matrix,
                                 const typename CommonTraits::DiscreteFunctionType& hmm_newton_rhs,
                                 const int hmm_iteration_step,
                                 const int loop_cycle,
                                 const double hmm_tolerance );


//! \TODO docme
HMMResult single_step( typename CommonTraits::GridPartType& gridPart,
        typename CommonTraits::GridPartType& gridPartFine,
        typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
        typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
        const typename CommonTraits::DiffusionType& diffusion_op,
        const Dune::RightHandSideAssembler< typename CommonTraits::DiscreteFunctionType >& rhsassembler,
        typename CommonTraits::DiscreteFunctionType& hmm_solution,
        const typename CommonTraits::DiscreteFunctionType& reference_solution,
        const int loop_cycle );

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {

#endif // ALGORITHM_STEP_HH
