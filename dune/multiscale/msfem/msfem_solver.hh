// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/fem/functions/checks.hh>

#include "dune/multiscale/common/la_backend.hh"
#include "dune/multiscale/msfem/localproblems/subgrid-list.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class MacroMicroGridSpecifier;

class Elliptic_MsFEM_Solver {
private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpace;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;
  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename DiscreteFunctionSpace::GridType LocalGrid;
  typedef typename LocalGrid::Traits::LeafIndexSet LocalGridLeafIndexSet;
  typedef typename LocalGrid::Traits::LeafIndexSet CoarseGridLeafIndexSet;
  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpace::IteratorType HostgridIterator;
  typedef typename HostgridIterator::Entity LocalEntity;
  typedef typename LocalEntity::EntityPointer LocalEntityPointer;
  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  static const int faceCodim = 1;

  // --------------------------- subgrid typedefs ------------------------------------
  typedef MsFEMTraits::LocalGridListType LocalGridListType;
  typedef MsFEMTraits::LocalGridType LocalGridType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionSpaceType SubgridDiscreteFunctionSpaceType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionType SubgridDiscreteFunctionType;
  //!-----------------------------------------------------------------------------------------

  typedef typename BackendChooser<DiscreteFunctionSpace>::LinearOperatorType MsLinearOperatorTypeType;
  typedef typename BackendChooser<DiscreteFunctionSpace>::InverseOperatorType InverseOperatorType;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;

public:
  Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace);

private:
  // create a hostgrid function from a subgridfunction (projection for global continuity)
  // Note: the maximum gride levels for both underlying grids must be the same
  void subgrid_to_hostrid_projection(const SubgridDiscreteFunctionType& sub_func,
                                     DiscreteFunctionType& host_func) const;

  //! identify fine scale part of MsFEM solution (including the projection!)
  // ------------------------------------------------------------------------------------
  void identify_fine_scale_part(MacroMicroGridSpecifier& specifier, MsFEMTraits::LocalGridListType& subgrid_list,
                                const DiscreteFunctionType& coarse_msfem_solution,
                                DiscreteFunctionType& fine_scale_part) const;

public:

  //! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  //! then:
  //! A --> diffusion operator ('DiffusionOperatorType')
  //! b --> advective part ('AdvectionTermType')
  //! c --> reaction part ('ReactionTermType')
  //! f --> 'first' source term, scalar ('SourceTermType')
  //! G --> 'second' source term, vector valued ('SecondSourceTermType')
  //! homogenous Dirchilet boundary condition!:
  void solve_dirichlet_zero(const CommonTraits::DiffusionType& diffusion_op, const CommonTraits::FirstSourceType& f,
                            // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                            // n(T)-layers.
                            MacroMicroGridSpecifier& specifier, MsFEMTraits::LocalGridListType& subgrid_list,
                            DiscreteFunctionType& coarse_scale_part, DiscreteFunctionType& fine_scale_part,
                            DiscreteFunctionType& solution) const;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef Elliptic_MSEM_Solver_HH
