#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/common/newton_rhs.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/fem/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/functions/norm.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <limits>
#include <sstream>
#include <string>

#include "dune/multiscale/common/dirichletconstraints.hh"
#include "fem_solver.hh"

namespace Dune {
namespace Multiscale {

Elliptic_FEM_Solver::Elliptic_FEM_Solver(const CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace)
  : discreteFunctionSpace_(discreteFunctionSpace) {}

void Elliptic_FEM_Solver::solve_linear(const CommonTraits::DiffusionType& diffusion_op,
                                       const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                                       const CommonTraits::FirstSourceType& f,
                                       CommonTraits::DiscreteFunctionType& solution, const bool use_smp) const {
  DSC_LOG_INFO << "Solving linear problem with standard FEM\n";

  // to assemble the computational time
  Dune::Timer assembleTimer;

  CommonTraits::LinearOperatorType fem_matrix("FEM stiffness matrix", discreteFunctionSpace_, discreteFunctionSpace_);
  if (use_smp) {
    FEM::SMPDiscreteEllipticOperator(discreteFunctionSpace_, diffusion_op).assemble_matrix(fem_matrix);
  } else {
    FEM::DiscreteEllipticOperator(discreteFunctionSpace_, diffusion_op, lower_order_term).assemble_matrix(fem_matrix);
  }

  DSC_LOG_INFO << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

  CommonTraits::DiscreteFunctionType fem_rhs("fem rhs", discreteFunctionSpace_);
  fem_rhs.clear();
  RightHandSideAssembler::assemble_fem(f, fem_rhs);

  const FEM::FEMTraits::InverseOperatorType inverse_op(
      fem_matrix, 1e-8, 1e-8, 5000, DSC_CONFIG_GET("global.cgsolver_verbose", false),
      DSC_CONFIG_GET("fem.algebraic_solver", "bcgs"), DSC_CONFIG_GET("fem.precond", "asm"), 1);
  fem_rhs.communicate();
  inverse_op.apply(fem_rhs, solution);

  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
  DSC_LOG_INFO << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
               << std::endl;
}

void
Elliptic_FEM_Solver::solve_nonlinear(const CommonTraits::DiffusionType& diffusion_op,
                                     const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                                     const CommonTraits::FirstSourceType& f,
                                     CommonTraits::DiscreteFunctionType& solution) const {
  DSC_LOG_INFO << "Solving non-linear problem." << std::endl;
  DSC_LOG_INFO << "Solving nonlinear problem with FEM + Newton-Method. Resolution level of grid = "
               << DSC_CONFIG_GET("fem.grid_level", 4) << "." << std::endl;
  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;

  Dune::Timer assembleTimer;
  //! residual vector
  // current residual
  CommonTraits::DiscreteFunctionType residual("FEM Newton Residual", discreteFunctionSpace_);
  residual.clear();

  const NewtonRightHandSide newton_rhs = {};
  CommonTraits::DiscreteFunctionType system_rhs("fem newton rhs", discreteFunctionSpace_);
  system_rhs.clear();

  typedef typename CommonTraits::DiscreteFunctionSpaceType::RangeType RangeType;

  RangeType relative_newton_error_finescale = std::numeric_limits<typename CommonTraits::RangeType>::max();

  //! (stiffness) matrix
  CommonTraits::LinearOperatorType fem_matrix("FEM stiffness matrix", discreteFunctionSpace_, discreteFunctionSpace_);

  FEM::DiscreteEllipticOperator discrete_elliptic_op(discreteFunctionSpace_, diffusion_op, lower_order_term);

  int iteration_step = 1;
  // the Newton step for the FEM reference problem (solved with Newton Method):
  // L2-Norm of residual < tolerance ?
  double tolerance = 1e-06;
  while (relative_newton_error_finescale > tolerance) {

    // (here: solution = solution from the last iteration step)
    DSC_LOG_INFO << "Newton iteration " << iteration_step << ":" << std::endl;
    Dune::Timer stepAssembleTimer;
    // assemble the stiffness matrix
    discrete_elliptic_op.assemble_jacobian_matrix(solution, fem_matrix);

    DSC_LOG_INFO << "Time to assemble FEM Newton stiffness matrix for current iteration: "
                 << stepAssembleTimer.elapsed() << "s" << std::endl;

    // assemble right hand side
    system_rhs.clear();
    newton_rhs.assemble_for_Newton_method(f, diffusion_op, *lower_order_term, solution, system_rhs);

    const auto rhs_L2_norm = DS::l2norm(system_rhs);
    if (rhs_L2_norm < 1e-10) {
      // residual solution almost identical to zero: break
      DSC_LOG_INFO << "Residual solution almost identical to zero. Therefore: break loop." << std::endl;
      DSC_LOG_INFO << "(L^2-Norm of current right hand side = " << rhs_L2_norm << " < 1e-10)" << std::endl;
      break;
    }
    // set Dirichlet Boundary to zero
    // --- boundary treatment ---
    // set the dirichlet points to zero (in righ hand side of the fem problem)
    Dune::Multiscale::getConstraintsFine(discreteFunctionSpace_).setValue(0.0, system_rhs);
    // --- end boundary treatment ---

    const FEM::FEMTraits::InverseOperatorType fem_newton_biCGStab(
        fem_matrix, 1e-8, 1e-8, 5000, true, "bcgs", DSC_CONFIG_GET("preconditioner_type", std::string("sor")));
    fem_newton_biCGStab.apply(system_rhs, residual);

    if (residual.dofsValid()) {
      solution += residual;
      relative_newton_error_finescale = DS::l2norm(residual);
      relative_newton_error_finescale /= DS::l2norm(solution);

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

  DSC_LOG_INFO << "Nonlinear Finite Element Problem solved with Newton Method in " << assembleTimer.elapsed() << "s."
               << std::endl << std::endl << std::endl;
}

void Elliptic_FEM_Solver::apply(const CommonTraits::DiffusionType& diffusion_op,
                                const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                                const CommonTraits::FirstSourceType& f, CommonTraits::DiscreteFunctionType& solution,
                                const bool use_smp) const {

  if (Problem::getModelData()->linear())
    solve_linear(diffusion_op, lower_order_term, f, solution, use_smp);
  else
    solve_nonlinear(diffusion_op, lower_order_term, f, solution);
}

} // namespace Multiscale
}
