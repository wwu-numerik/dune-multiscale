// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef Elliptic_MSEM_Solver_HH
#define Elliptic_MSEM_Solver_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

#include <dune/common/fmatrix.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/solver/petscsolver.hh>
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>

#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class Elliptic_MsFEM_Solver {
private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpace;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;
  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename DiscreteFunctionSpace::GridType HostGrid;
  typedef typename HostGrid::Traits::LeafIndexSet HostGridLeafIndexSet;
  typedef typename HostGrid::Traits::LeafIndexSet CoarseGridLeafIndexSet;
  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpace::IteratorType HostgridIterator;
  typedef typename HostgridIterator::Entity HostEntity;
  typedef typename HostEntity::EntityPointer HostEntityPointer;
  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  static const int faceCodim = 1;

  // --------------------------- subgrid typedefs ------------------------------------
  typedef MsFEMTraits::SubGridListType SubGridListType;
  typedef MsFEMTraits::SubGridType SubGridType;
  typedef MsFEMTraits::SubGridDiscreteFunctionSpaceType SubgridDiscreteFunctionSpaceType;
  typedef MsFEMTraits::SubGridDiscreteFunctionType SubgridDiscreteFunctionType;
  //!-----------------------------------------------------------------------------------------

  typedef Dune::Fem::PetscLinearOperator<DiscreteFunctionType, DiscreteFunctionType> MsFEMMatrixType;
  typedef Dune::Fem::PetscInverseOperator<DiscreteFunctionType, MsFEMMatrixType> InverseOperatorType;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;

public:
  Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace);

private:
  // create a hostgrid function from a subgridfunction (projection for global continuity)
  // Note: the maximum gride levels for both underlying grids must be the same
  void subgrid_to_hostrid_projection(const SubgridDiscreteFunctionType& sub_func,
                                     DiscreteFunctionType& host_func) const;

  //! copy coarse scale part of MsFEM solution into a function defined on the fine grid
  // ------------------------------------------------------------------------------------
  void projectCoarseToFineScale(MacroMicroGridSpecifier& specifier, const DiscreteFunctionType& coarse_msfem_solution,
                                DiscreteFunctionType& coarse_scale_part) const;

  //! identify fine scale part of MsFEM solution (including the projection!)
  // ------------------------------------------------------------------------------------
  void identify_fine_scale_part(MacroMicroGridSpecifier& specifier, MsFEMTraits::SubGridListType& subgrid_list,
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
                            MacroMicroGridSpecifier& specifier, MsFEMTraits::SubGridListType& subgrid_list,
                            DiscreteFunctionType& coarse_scale_part, DiscreteFunctionType& fine_scale_part,
                            DiscreteFunctionType& solution) const;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef Elliptic_MSEM_Solver_HH
