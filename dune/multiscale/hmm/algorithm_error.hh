// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef ALGORITHM_ERROR_HH
#define ALGORITHM_ERROR_HH

#include <dune/multiscale/hmm/hmm_traits.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

struct HMMResult;

//! Error Estimation
HMMResult estimate_error(
        const typename HMMTraits::GridPartType& gridPart,
        const typename HMMTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
        const typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
        const typename HMMTraits::DiffusionType& diffusion_op,
        const typename HMMTraits::CellProblemNumberingManagerType& cp_num_manager,
        const typename HMMTraits::DiscreteFunctionType& hmm_solution
         );

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {

#endif // ALGORITHM_ERROR_HH
