#include "fem_solver.hh"

#include <dune/common/fmatrix.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#include <dune/multiscale/tools/righthandside_assembler.hh>
#include <dune/multiscale/fem/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

//! define a dummy mass term:
template< class FunctionSpaceImp >
class DummyMass
  : public Dune::Fem::Function< FunctionSpaceImp, DummyMass< FunctionSpaceImp > >
{
private:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef DummyMass< FunctionSpaceType >                     ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& y) const {
    DUNE_THROW(Dune::InvalidStateException, "Do not use the Dummy-class!!!");
    y = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType /*time*/,
                       RangeType& y) const {
    DSC_LOG_ERROR << "WARNING! Wrong call for 'evaluate' method of the MassTerm class (evaluate(x,t,y)). Return 0.0."
              << std::endl;
    return evaluate(x, y);
  }
};


Elliptic_FEM_Solver::Elliptic_FEM_Solver(const Elliptic_FEM_Solver::DiscreteFunctionSpace& discreteFunctionSpace)
    : discreteFunctionSpace_(discreteFunctionSpace)
{}


//! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
//! then:
//! A --> diffusion operator ('DiffusionOperatorType')
//! b --> advective part ('AdvectionTermType')
//! c --> reaction part ('ReactionTermType')
//! f --> 'first' source term, scalar ('SourceTermType')
//! G --> 'second' source term, vector valued ('SecondSourceTermType')
//! homogenous Dirchilet boundary condition!:
void Elliptic_FEM_Solver::solve_dirichlet_zero(const CommonTraits::DiffusionType& diffusion_op,
                          const CommonTraits::FirstSourceType& f,
                          Elliptic_FEM_Solver::DiscreteFunction& solution) const
{
    const GridPart& gridPart = discreteFunctionSpace_.gridPart();

    //! define the right hand side assembler tool
    // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
    RightHandSideAssembler< DiscreteFunctionType > rhsassembler;

    //! define the discrete (elliptic) operator that describes our problem
    // ( effect of the discretized differential operator on a certain discrete function )
    Multiscale::FEM::DiscreteEllipticOperator< DiscreteFunction, CommonTraits::DiffusionType, DummyMassType > discrete_elliptic_op(
                discreteFunctionSpace_,
                diffusion_op);
    // discrete elliptic operator (corresponds with FEM Matrix)

    //! (stiffness) matrix
    FEMMatrix fem_matrix("FEM stiffness matrix", discreteFunctionSpace_, discreteFunctionSpace_);

    //! right hand side vector
    // right hand side for the finite element method:
    DiscreteFunction fem_rhs("fem newton rhs", discreteFunctionSpace_);
    fem_rhs.clear();

    DSC_LOG_INFO << "Solving linear problem with standard FEM and resolution level "
                 << discreteFunctionSpace_.grid().maxLevel() << "." << std::endl;
    DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

    // to assemble the computational time
    Dune::Timer assembleTimer;
    // assemble the stiffness matrix
    discrete_elliptic_op.assemble_matrix(fem_matrix);

    DSC_LOG_INFO << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

    // assemble right hand side
    rhsassembler.assemble< 2* DiscreteFunctionSpace::polynomialOrder + 2 >(f, fem_rhs);

    // oneLinePrint( DSC_LOG_DEBUG , fem_rhs );

    // --- boundary treatment ---
    // set the dirichlet points to zero (in righ hand side of the fem problem)
    typedef typename DiscreteFunctionSpace::IteratorType EntityIterator;
    EntityIterator endit = discreteFunctionSpace_.end();
    for (EntityIterator it = discreteFunctionSpace_.begin(); it != endit; ++it)
    {
        IntersectionIterator iit = gridPart.ibegin(*it);
        const IntersectionIterator endiit = gridPart.iend(*it);
        for ( ; iit != endiit; ++iit)
        {
            if ( !(*iit).boundary() )
                continue;

            LocalFunction rhsLocal = fem_rhs.localFunction(*it);
            const LagrangePointSet& lagrangePointSet
                    = discreteFunctionSpace_.lagrangePointSet(*it);

            const int face = (*iit).indexInInside();

            auto faceIterator
                    = lagrangePointSet.beginSubEntity< faceCodim >(face);
            const auto faceEndIterator
                    = lagrangePointSet.endSubEntity< faceCodim >(face);
            for ( ; faceIterator != faceEndIterator; ++faceIterator)
                rhsLocal[*faceIterator] = 0;
        }
    }
    // --- end boundary treatment ---

    InverseFEMMatrix fem_biCGStab(fem_matrix, 1e-8, 1e-8, 20000, true);
    fem_biCGStab(fem_rhs, solution);

    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
    DSC_LOG_INFO << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl
                 << std::endl << std::endl;

    // oneLinePrint( DSC_LOG_DEBUG , solution );
} // solve_dirichlet_zero


} // namespace Multiscale
}


