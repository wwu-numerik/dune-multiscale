// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH
#define DUNE_RIGHT_HAND_SIDE_ASSEMBLER_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {

//! Assembler for right rand side
//! We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
//! we call f the first Source and G the second Source
namespace MsFEM {
class MacroMicroGridSpecifier;
class LocalGridList;
}  // namespace MsFEM

class RightHandSideAssembler {
private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef DomainFieldType TimeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  static const int dimension = GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;
  static const int quadratureOrder = 2*polynomialOrder + 2;

public:
  /** assemble standard right hand side:
   * if there is only one source (f) (there is no second source):
   * discreteFunction is an output parameter (kind of return value)
   **/
  static void assemble_fem(const CommonTraits::FirstSourceType& f, DiscreteFunctionType& rhsVector);

  /**
   * The rhs-assemble()-methods for linear elliptic problems
   * with non-homogeneous Dirichlet and Neumann boundary conditions:
   **/
  static void
  assemble_hmm_lod(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
           const DiscreteFunctionType& dirichlet_extension,
           const CommonTraits::NeumannBCType& neumann_bc, DiscreteFunctionType& rhsVector);

  /** assemble right hand side (if there is only one source - f):
   *  assemble-method for MsFEM in symmetric (non-Petrov-Galerkin) formulation
   *  rhsVector is the output parameter (kind of return value)
   **/
  static void assemble_msfem(const CommonTraits::FirstSourceType& f, MsFEM::MacroMicroGridSpecifier& specifier,
                                           MsFEM::LocalGridList& subgrid_list, DiscreteFunctionType& rhsVector);
};  // end class
} // end namespace Multiscale
} // end namespace Dune

#endif
