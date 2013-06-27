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

#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {


class Elliptic_MsFEM_Solver
{
private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef DiscreteFunctionType DiscreteFunction;

  typedef typename DiscreteFunction::FunctionSpaceType FunctionSpace;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;

  typedef typename DiscreteFunctionSpace::GridType HostGrid;

  typedef typename HostGrid::Traits::LeafIndexSet HostGridLeafIndexSet;

  typedef typename HostGrid::Traits::LeafIndexSet CoarseGridLeafIndexSet;

  typedef typename DiscreteFunctionSpace::DomainType        DomainType;
  typedef typename DiscreteFunctionSpace::RangeType         RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;

  // typedef typename HostGrid ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // LevelEntityIteratorType;

  typedef typename DiscreteFunctionSpace::IteratorType HostgridIterator;

  typedef typename HostgridIterator::Entity HostEntity;

  typedef typename HostEntity::EntityPointer HostEntityPointer;

  // typedef typename HostGrid :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // HostGridLevelEntityIterator;

  enum { faceCodim = 1 };

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  // --------------------------- subgrid typedefs ------------------------------------

  typedef SubGrid< HostGrid::dimension, HostGrid > SubGridType;

  typedef Fem::LeafGridPart< SubGridType > SubGridPart;

  // typedef typename SubGridType ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // SubGridLevelEntityIteratorType;

  typedef Fem::LagrangeDiscreteFunctionSpace< FunctionSpace, SubGridPart, 1 >  // 1=POLORDER
    SubgridDiscreteFunctionSpace;

  typedef Fem::AdaptiveDiscreteFunction< SubgridDiscreteFunctionSpace > SubgridDiscreteFunction;

  typedef typename SubgridDiscreteFunctionSpace::IteratorType CoarseGridIterator;

  typedef typename CoarseGridIterator::Entity CoarseGridEntity;

  typedef typename CoarseGridEntity::EntityPointer CoarseGridEntityPointer;

  typedef typename SubgridDiscreteFunction::LocalFunctionType CoarseGridLocalFunction;

  typedef typename SubgridDiscreteFunctionSpace::LagrangePointSetType
  CoarseGridLagrangePointSet;


  //!-----------------------------------------------------------------------------------------

  //! --------------------- the standard matrix traits -------------------------------------

  struct MatrixTraits
  {
    typedef SubgridDiscreteFunctionSpace                          RowSpaceType;
    typedef SubgridDiscreteFunctionSpace                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                          StencilType;
    typedef Fem::ParallelScalarProduct< SubgridDiscreteFunctionSpace > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  //! --------------------------------------------------------------------------------------

  //! --------------------- type of fem stiffness matrix -----------------------------------

  typedef Fem::SparseRowMatrixOperator< DiscreteFunction, DiscreteFunction, MatrixTraits > MsFEMMatrix;

  //! --------------------------------------------------------------------------------------

  //! --------------- solver for the linear system of equations ----------------------------

  // use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric
  // ones. GMRES seems to be more stable, but is extremely slow!
  typedef /*OEMBICGSQOp*//*CGInverseOp*/ Fem::OEMBICGSTABOp< DiscreteFunction, MsFEMMatrix > InverseMsFEMMatrix;

  //! --------------------------------------------------------------------------------------

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;

public:
  Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace);

private:
  // create a hostgrid function from a subgridfunction (projection for global continuity)
  // Note: the maximum gride levels for both underlying grids must be the same
  void subgrid_to_hostrid_projection(const SubgridDiscreteFunction& sub_func, DiscreteFunction& host_func) const;


  //! copy coarse scale part of MsFEM solution into a function defined on the fine grid
  // ------------------------------------------------------------------------------------
  void identify_coarse_scale_part(MacroMicroGridSpecifier &specifier,
                                   const DiscreteFunction& coarse_msfem_solution,
                                   DiscreteFunction& coarse_scale_part ) const;


  //! identify fine scale part of MsFEM solution (including the projection!)
  // ------------------------------------------------------------------------------------
  void identify_fine_scale_part(MacroMicroGridSpecifier &specifier,
                                                          MsFEMTraits::SubGridListType& subgrid_list,
                                                          const DiscreteFunction& coarse_msfem_solution,
                                                          DiscreteFunction& fine_scale_part ) const;

public:

  //! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  //! then:
  //! A --> diffusion operator ('DiffusionOperatorType')
  //! b --> advective part ('AdvectionTermType')
  //! c --> reaction part ('ReactionTermType')
  //! f --> 'first' source term, scalar ('SourceTermType')
  //! G --> 'second' source term, vector valued ('SecondSourceTermType')
  //! homogenous Dirchilet boundary condition!:
  void solve_dirichlet_zero(const CommonTraits::DiffusionType& diffusion_op,
                            const CommonTraits::FirstSourceType& f,
                            // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                            // n(T)-layers.
                            MacroMicroGridSpecifier &specifier,
                            MsFEMTraits::SubGridListType& subgrid_list,
                            DiscreteFunction& coarse_scale_part,
                            DiscreteFunction& fine_scale_part,
                            DiscreteFunction& solution) const;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef Elliptic_MSEM_Solver_HH
