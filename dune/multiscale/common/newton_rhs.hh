#ifndef DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH
#define DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH

#include <config.h>

#include <dune/multiscale/common/traits.hh>

#include <dune/fem/quadrature/quadrature.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/cell_problem_solver.hh>

#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/common/dirichletconstraints.hh>


#include <dune/stuff/common/logging.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace Multiscale {

//! Assembler for right rand side
//! We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
//! we call f the first Source and G the second Source

class NewtonRightHandSide {
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef DomainFieldType TimeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename GridType::Codim<0>::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef MsFEM::LocalSolutionManager LocalSolutionManagerType;

  static const int dimension = GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;
  static const int quadratureOrder = 2*polynomialOrder + 2;
public:
/**
 * The rhs-assemble()-methods for non-linear elliptic problems:
 * discreteFunction is an output parameter (kind of return value)
 **/
static void assemble_for_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                                       const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                       DiscreteFunctionType& rhsVector); // end method


/**
 * The rhs-assemble()-methods for non-linear elliptic problems
 * if there is a first source f and a lower order term F:
 * discreteFunction is an output parameter (kind of return value)
 **/
static void assemble_for_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                                       const CommonTraits::LowerOrderTermType& F,
                                       const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                                       DiscreteFunctionType& rhsVector); // end method

/**
 * The rhs-assemble()-methods for non-linear elliptic problems
 * if there is a first source f, a lower order term F
 * and Dirichlet and Neumann boundary conditions
 * discreteFunction is an output parameter (kind of return value)
 * \param old_u_H from the last iteration step
 * \param dirichlet_extension discrete function describing dirichlet extension
 **/
static void assemble_for_Newton_method(
    const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A, const CommonTraits::LowerOrderTermType& F,
    const DiscreteFunctionType& old_u_H,
    const DiscreteFunctionType& dirichlet_extension,
    const CommonTraits::NeumannBCType& neumann_bc, DiscreteFunctionType& rhsVector); // end method

//! The rhs-assemble()-methods for non-linear elliptic problems, solved with the heterogenous multiscale method
// ( requires reconstruction of old_u_H and local fine scale averages )

static void
assemble_for_HMM_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                               const DiscreteFunctionType& old_u_H, // old_u_H from the last iteration step
                               // to obtain some information about the periodic discrete function space (space
                               // for the cell problems)
                               const HMM::CellProblemNumberingManager& cp_num_manager,
                               const HMM::HMMTraits::PeriodicDiscreteFunctionType& dummy_func, DiscreteFunctionType& rhsVector); // end method
};

} // end namespace Multiscale
} // end namespace Dune
#endif // DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH
