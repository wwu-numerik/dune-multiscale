// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MS_Elliptic_FEM_Solver_HH
#define MS_Elliptic_FEM_Solver_HH

#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

//! define a dummy mass term:
template< class FunctionSpaceImp >
class DummyMass;

//! \todo docme
class Elliptic_FEM_Solver
{
private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef DiscreteFunctionType DiscreteFunction;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;

  enum { faceCodim = 1 };

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  typedef typename LagrangePointSet::template Codim< faceCodim >
    ::SubEntityIteratorType
  FaceDofIterator;

  typedef DummyMass< DiscreteFunctionSpace > DummyMassType;

  //! --------------------- the standard matrix traits -------------------------------------

  struct MatrixTraits
  {
    typedef DiscreteFunctionSpace                          RowSpaceType;
    typedef DiscreteFunctionSpace                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                   StencilType;
    typedef ParallelScalarProduct< DiscreteFunctionSpace > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  //! --------------------------------------------------------------------------------------

  //! --------------------- type of fem stiffness matrix -----------------------------------

  typedef SparseRowMatrixOperator< DiscreteFunction, DiscreteFunction, MatrixTraits > FEMMatrix;

  //! --------------------------------------------------------------------------------------

  //! --------------- solver for the linear system of equations ----------------------------

  // use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric
  // ones.
  // GMRES seems to be more stable, but is extremely slow!
  // typedef OEMBICGSQOp/*OEMBICGSTABOp*/< DiscreteFunction, FEMMatrix > InverseFEMMatrix;
  typedef CGInverseOperator< DiscreteFunction, FEMMatrix > InverseFEMMatrix;

  //! --------------------------------------------------------------------------------------

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;

public:
  Elliptic_FEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace);

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
                            DiscreteFunction& solution) const;
};

} // namespace Multiscale
} // namespace Dune

#endif // #ifndef MS_Elliptic_FEM_Solver_HH
