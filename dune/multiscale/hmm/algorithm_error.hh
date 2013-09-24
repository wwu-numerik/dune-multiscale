// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef ALGORITHM_ERROR_HH
#define ALGORITHM_ERROR_HH

#include <dune/multiscale/hmm/hmm_traits.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {

class CellProblemNumberingManager;

namespace Multiscale {
namespace HMM {

struct HMMResult;

//! Error Estimation
HMMResult estimate_error(const typename CommonTraits::GridPartType& gridPart,
                         const typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
                         const typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                         const typename CommonTraits::DiffusionType& diffusion_op,
                         const CellProblemNumberingManager& cp_num_manager,
                         const typename CommonTraits::DiscreteFunctionType& hmm_solution);

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {

#endif // ALGORITHM_ERROR_HH
