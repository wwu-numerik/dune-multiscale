// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
#define DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {

namespace MsFEM {
class LocalGridList;
} // namespace MsFEM

/** Assembler for right rand side
 * We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
 * we call f the first Source and G the second Source
 **/
class RightHandSideAssembler {
private:
  typedef typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

public:
  /** assemble right hand side (if there is only one source - f):
   *  assemble-method for MsFEM in symmetric (non-Petrov-Galerkin) formulation
   *  rhsVector is the output parameter (kind of return value)
   **/
  static void assemble_msfem(const CommonTraits::DiscreteFunctionSpaceType& coarse_space,
                             const CommonTraits::SourceType& f, MsFEM::LocalGridList& subgrid_list,
                             CommonTraits::DiscreteFunctionType& rhsVector);
private:
  static std::vector<std::vector<RangeType>> coarseBaseEvals_;
  static std::vector<std::vector<JacobianRangeType>> coarseBaseJacs_;
  static bool cached_;

}; // end class
} // end namespace Multiscale
} // end namespace Dune

#endif
