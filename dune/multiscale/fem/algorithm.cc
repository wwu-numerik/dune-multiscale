#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/newton_rhs.hh>
#include <dune/multiscale/common/output_traits.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/fem/print_info.hh>
#include <dune/multiscale/problems/selector.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

#include <dune/fem/misc/l2norm.hh>

#include <string>
#include <fstream>

#include "fem_traits.hh"

namespace {
const std::string seperator_line =
    "---------------------------------------------------------------------------------\n";
}

namespace Dune {
namespace Multiscale {
namespace FEM {


//! \TODO docme
void solve(typename CommonTraits::DiscreteFunctionType& solution,
           const typename CommonTraits::DiscreteFunctionType& dirichlet_extension,
           const typename CommonTraits::DiscreteFunctionSpaceType& finerDiscreteFunctionSpace,
           const typename FEMTraits::EllipticOperatorType& discrete_elliptic_op) {
  typename CommonTraits::LinearOperatorType system_matrix("FEM Newton stiffness matrix", finerDiscreteFunctionSpace,
                                                          finerDiscreteFunctionSpace);
  typename CommonTraits::DiscreteFunctionType system_rhs("fem newton rhs", finerDiscreteFunctionSpace);
  system_rhs.clear();
  const auto f = Dune::Multiscale::Problem::getFirstSource();
  const auto diffusion_op = Problem::getDiffusion();
  const auto neumann_bc = Problem::getNeumannBC();

  DSC_LOG_INFO << "Solving linear problem." << std::endl;
  DSC_LOG_INFO << "Solving linear problem with standard FEM and resolution level "
               << DSC_CONFIG_GET("fem.grid_level", 4) << "." << std::endl;
  DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

  // to assemble the computational time
  Dune::Timer assembleTimer;

  // assemble the stiffness matrix
  discrete_elliptic_op.assemble_matrix(system_matrix);

  DSC_LOG_INFO << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

  const RightHandSideAssembler rhsassembler = {};
  rhsassembler.assemble(*f, *diffusion_op, dirichlet_extension, *neumann_bc, system_rhs);

  const auto boundaryInfo = Problem::getModelData()->boundaryInfo();
  DirichletConstraints<typename CommonTraits::DiscreteFunctionSpaceType> constraints(*boundaryInfo, finerDiscreteFunctionSpace);
  constraints.applyToOperator(system_matrix);
  constraints(dirichlet_extension, system_rhs);

  const typename FEMTraits::InverseOperatorType fem_biCGStab(
      system_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("global.cgsolver_verbose", false),
      DSC_CONFIG_GET("fem.algebraic_solver", "bcgs"), DSC_CONFIG_GET("preconditioner_type", std::string("sor")));
  fem_biCGStab(system_rhs, solution);

  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
  DSC_LOG_INFO << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
               << std::endl;
}


//! the main FEM computation
void algorithm(typename CommonTraits::GridPointerType& macro_grid_pointer, const std::string filename) {
  using namespace Dune;

  const auto problem_data = Problem::getModelData();
  print_info(*problem_data, DSC_LOG_INFO);
  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for the finite element problem
  typename CommonTraits::GridPartType gridPart(*macro_grid_pointer);
  //! --------------------------------------------------------------------------------------

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  typename CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  //! --------------------------------------------------------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const auto diffusion_op = Problem::getDiffusion();
  // lower order term F(x, u(x), grad u(x) )
  const auto lower_order_term = Problem::getLowerOrderTerm();
  // Dirichlet boundary condition
  const auto dirichlet_bc = Problem::getDirichletBC();

  // discrete function that takes the values of the Dirichlet BC on the Dirichlet Boundary nodes
  // and that is zero elsewhere
  typename CommonTraits::DiscreteFunctionType dirichlet_extension("Dirichlet extension", discreteFunctionSpace);
  dirichlet_extension.clear();
  setDirichletValues(*dirichlet_bc, dirichlet_extension);

  //! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
  const typename FEMTraits::EllipticOperatorType discrete_elliptic_op(discreteFunctionSpace, *diffusion_op,
                                                                      lower_order_term);

  //! solution vector
  // - By solution, we denote the "discrete solution" determined with FEM in the linear case or FEM-Newton (nonlinear
  // case)
  //    ( if the elliptic problem is linear, the 'solution' is determined without the Newton method )
  // - solution of the finite element method, where we use the Newton method to solve the non-linear system of equations
  auto discrete_solution = make_df_ptr<typename CommonTraits::DiscreteFunctionType>(filename + " FEM(-Newton) Solution",
                                                                discreteFunctionSpace);
  discrete_solution->clear();

  if (DSC_CONFIG_GET("problem.linear", true))
    solve(*discrete_solution, dirichlet_extension, discreteFunctionSpace, discrete_elliptic_op);
  else
    solve_nonlinear(*discrete_solution, dirichlet_extension, discreteFunctionSpace, discrete_elliptic_op,
                    lower_order_term, filename);
  *discrete_solution += dirichlet_extension;

  // write FEM solution to a file and produce a VTK output
  write_discrete_function(discrete_solution, "fem");

  ErrorCalculator(nullptr, discrete_solution.get()).print(DSC_LOG_INFO_0);
}

//! \TODO docme
void solve_nonlinear(typename CommonTraits::DiscreteFunctionType& solution,
                     const typename CommonTraits::DiscreteFunctionType& dirichlet_extension,
                     const typename CommonTraits::DiscreteFunctionSpaceType& finerDiscreteFunctionSpace,
                     const typename FEMTraits::EllipticOperatorType& discrete_elliptic_op,
                     const std::unique_ptr<const CommonTraits::LowerOrderTermType> &lower_order_term,
                     const std::string& filename)
{
  typename CommonTraits::LinearOperatorType system_matrix("FEM Newton stiffness matrix", finerDiscreteFunctionSpace,
                                                          finerDiscreteFunctionSpace);
  typename CommonTraits::DiscreteFunctionType system_rhs("fem newton rhs", finerDiscreteFunctionSpace);
  system_rhs.clear();
  const auto f = Dune::Multiscale::Problem::getFirstSource();
  const auto diffusion_op = Problem::getDiffusion();
  const auto neumann_bc = Problem::getNeumannBC();
  DSC_LOG_INFO << "Solving non-linear problem." << std::endl;
  DSC_LOG_INFO << "Solving nonlinear problem with FEM + Newton-Method. Resolution level of grid = "
               << DSC_CONFIG_GET("fem.grid_level", 4) << "." << std::endl;
  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;

  Dune::Timer assembleTimer;
  //! residual vector
  // current residual
  typename CommonTraits::DiscreteFunctionType residual(filename + "FEM Newton Residual", finerDiscreteFunctionSpace);
  residual.clear();

  typename CommonTraits::RangeType relative_newton_error_finescale = std::numeric_limits<typename CommonTraits::RangeType>::max();
  typename CommonTraits::RangeType rhs_L2_norm = std::numeric_limits<typename CommonTraits::RangeType>::max();

  int iteration_step = 1;
  // the Newton step for the FEM reference problem (solved with Newton Method):
  // L2-Norm of residual < tolerance ?
  double tolerance = 1e-06;
  const NewtonRightHandSide newton_rhs = {};
  while (relative_newton_error_finescale > tolerance) {
    // (here: solution = solution from the last iteration step)
    DSC_LOG_INFO << "Newton iteration " << iteration_step << ":" << std::endl;
    Dune::Timer stepAssembleTimer;
    // assemble the stiffness matrix
    discrete_elliptic_op.assemble_jacobian_matrix(solution, dirichlet_extension, system_matrix);

    DSC_LOG_INFO << "Time to assemble FEM Newton stiffness matrix for current iteration: "
                 << stepAssembleTimer.elapsed() << "s" << std::endl;

    newton_rhs.assemble_for_Newton_method(*f, *diffusion_op, *lower_order_term, solution,
                                                          dirichlet_extension, *neumann_bc, system_rhs);

    const Dune::Fem::L2Norm<typename CommonTraits::DiscreteFunctionType::GridPartType> l2norm(system_rhs.gridPart());
    rhs_L2_norm = l2norm.norm(system_rhs);
    if (rhs_L2_norm < 1e-10) {
      // residual solution almost identical to zero: break
      DSC_LOG_DEBUG << "Residual solution almost identical to zero. Therefore: break loop." << std::endl;
      DSC_LOG_DEBUG << "(L^2-Norm of current right hand side = " << rhs_L2_norm << " < 1e-10)" << std::endl;
      break;
    }
    // set Dirichlet Boundary to zero
    boundaryTreatment(system_rhs);

    const typename FEMTraits::InverseOperatorType fem_newton_biCGStab(system_matrix, 1e-8, 1e-8, 20000, true);
    fem_newton_biCGStab(system_rhs, residual);

    if (residual.dofsValid()) {
      solution += residual;
      relative_newton_error_finescale = l2norm.norm(residual);
      relative_newton_error_finescale /= l2norm.norm(solution);

      DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
      // residual solution almost identical to zero: break
      DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
      if (relative_newton_error_finescale <= tolerance) {
        DSC_LOG_INFO << "Since tolerance = " << tolerance << ": break loop." << std::endl;
      }
      residual.clear();
    } else {
      DSC_LOG_INFO << "WARNING! Invalid dofs in 'residual'." << std::endl;
      break;
    }
    iteration_step += 1;
  }

  if (DSC_CONFIG_GET("problem.linear", true)) {
    DSC_LOG_INFO << "Finite Element Problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
                 << std::endl;
  } else {
    DSC_LOG_INFO << "Nonlinear Finite Element Problem solved with Newton Method in " << assembleTimer.elapsed()
                 << "s." << std::endl << std::endl << std::endl;
  }
  DSC_LOG_INFO << seperator_line;

}

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
