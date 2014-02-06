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
class CellProblemNumberingManager;

/**
 * \return true if a program should continue with a new newton step
 **/
bool process_hmm_newton_residual(typename CommonTraits::RangeType& relative_newton_error,
                                 typename CommonTraits::DiscreteFunction_ptr& hmm_solution,
                                 CommonTraits::LinearOperatorType& hmm_newton_matrix,
                                 const typename CommonTraits::DiscreteFunctionType& hmm_newton_rhs,
                                 const int hmm_iteration_step, const int loop_cycle, const double hmm_tolerance);

//! \TODO docme
HMMResult single_step(typename CommonTraits::GridPartType& gridPart, typename CommonTraits::GridPartType& gridPartFine,
                      typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
                      typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                      const typename CommonTraits::DiffusionType& diffusion_op,
                      CommonTraits::DiscreteFunction_ptr& hmm_solution,
                      const typename CommonTraits::DiscreteFunctionType& reference_solution, const int loop_cycle);

//! The rhs-assemble()-methods for non-linear elliptic problems, solved with the heterogenous multiscale method
// ( requires reconstruction of old_u_H and local fine scale averages )

/**
 * old_u_H from the last iteration step to obtain some information about the
 * periodic discrete function space (space for the cell problems)
**/
static void assemble_for_HMM_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                                           const CommonTraits::DiscreteFunctionType& old_u_H,
                                           const CellProblemNumberingManager& cp_num_manager,
                                           const HMMTraits::PeriodicDiscreteFunctionType& dummy_func,
                                           CommonTraits::DiscreteFunctionType& rhsVector);
} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {

#endif // ALGORITHM_STEP_HH
