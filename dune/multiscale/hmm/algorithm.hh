// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_HMM_ALGORITHM_HH
#define DUNE_MS_HMM_ALGORITHM_HH

#include <dune/multiscale/hmm/hmm_traits.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

struct HMMResult;

//! \TODO docme
bool adapt(const HMMResult& result,
           const int loop_cycle,
           const double error_tolerance_,
           const typename HMMTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
           typename HMMTraits::AdaptationManagerType& adaptationManager
           );

//! the main hmm computation
void algorithm(typename HMMTraits::GridPointerType& macro_grid_pointer,   // grid pointer that belongs to the macro grid
               typename HMMTraits::GridPointerType& fine_macro_grid_pointer,   // grid pointer that belongs to the fine macro grid (for
                                                           // reference computations)
               typename HMMTraits::GridPointerType& periodic_grid_pointer,   // grid pointer that belongs to the periodic micro grid
               const std::string filename);

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {

#endif // DUNE_MS_HMM_ALGORITHM_HH
