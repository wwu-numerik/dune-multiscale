#ifndef DUNE_FEM_ALGORITHM_HH
#define DUNE_FEM_ALGORITHM_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>

#include <dune/stuff/common/logging.hh>

#include <string>
#include <fstream>

#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>

namespace {
  const std::string seperator_line = "---------------------------------------------------------------------------------\n";
}

//! set the dirichlet points to zero
template< class DiscreteFunctionType >
void boundaryTreatment(DiscreteFunctionType& rhs) {
  using namespace Dune::Stuff;
  const auto& discreteFunctionSpace = rhs.space();
  static const unsigned int faceCodim = 1;
  for (const auto& entity : discreteFunctionSpace)
  {
    for (const auto& intersection
         : Dune::Stuff::Common::intersectionRange(discreteFunctionSpace.gridPart(), entity))
    {
      if ( !intersection.boundary() )
        continue;
      auto rhsLocal = rhs.localFunction(entity);
      const auto face = intersection.indexInInside();
      for(auto point
          : Dune::Stuff::Common::lagrangePointSetRange<faceCodim>(rhs.space(), entity, face))
        rhsLocal[point] = 0;
    }
  }
} // boundaryTreatment


template <class FEM>
void solve(typename FEM::DiscreteFunctionType& solution,
                 const typename FEM::DiscreteFunctionSpaceType& finerDiscreteFunctionSpace,
                 const typename FEM::EllipticOperatorType& discrete_elliptic_op,
                 const std::string& filename,
                 const Dune::RightHandSideAssembler< typename FEM::DiscreteFunctionType >& rhsassembler)
{
  static const int fem_polorder = 2* FEM::DiscreteFunctionSpaceType::polynomialOrder + 2;

  //! *************************** Assembling the problem ****************************

  //! (stiffness) matrix
  typename FEM::FEMMatrix system_matrix("FEM Newton stiffness matrix", finerDiscreteFunctionSpace, finerDiscreteFunctionSpace);

  //! right hand side vector
  // right hand side for the finite element method with Newton solver:
  // ( also right hand side for the finer discrete function space )
  typename FEM::DiscreteFunctionType system_rhs("fem newton rhs", finerDiscreteFunctionSpace);
  system_rhs.clear();

  const typename FEM::FirstSourceType f;   // standard source f

  if (DSC_CONFIG_GET("problem.linear", true))
  {
    DSC_LOG_INFO << "Solving linear problem." << std::endl;
    DSC_LOG_INFO << "Solving linear problem with standard FEM and resolution level "
              << typename FEM::ModelProblemDataType().getRefinementLevelReferenceProblem() << "." << std::endl;
    DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

    // to assemble the computational time
    Dune::Timer assembleTimer;

    // assemble the stiffness matrix
    discrete_elliptic_op.assemble_matrix( system_matrix );

    DSC_LOG_INFO << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

    // assemble right hand side
    rhsassembler.template assemble< fem_polorder >(f, system_rhs);

    // set Dirichlet Boundary to zero
    boundaryTreatment(system_rhs);

    const typename FEM::InverseFEMMatrix fem_biCGStab(system_matrix, 1e-8, 1e-8, 20000, VERBOSE);
    fem_biCGStab(system_rhs, solution);

    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
    DSC_LOG_INFO << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
              << std::endl;
  } else {
    DSC_LOG_INFO << "Solving non-linear problem." << std::endl;
    DSC_LOG_INFO << "Solving nonlinear problem with FEM + Newton-Method. Resolution level of grid = "
              << typename FEM::ModelProblemDataType().getRefinementLevelReferenceProblem() << "." << std::endl;
    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;

    Dune::Timer assembleTimer;
    //! residual vector
    // current residual
    typename FEM::DiscreteFunctionType residual(filename + "FEM Newton Residual", finerDiscreteFunctionSpace);
    residual.clear();

    typename FEM::RangeType relative_newton_error_finescale = 10000.0;
    typename FEM::RangeType rhs_L2_norm = 10000.0;

    int iteration_step = 1;
    // the Newton step for the FEM reference problem (solved with Newton Method):
    // L2-Norm of residual < tolerance ?
    double tolerance = 1e-06;
    while (relative_newton_error_finescale > tolerance)
    {
      // (here: solution = solution from the last iteration step)
      DSC_LOG_INFO << "Newton iteration " << iteration_step << ":" << std::endl;
      Dune::Timer stepAssembleTimer;
      // assemble the stiffness matrix
      discrete_elliptic_op.assemble_jacobian_matrix(solution, system_matrix);

      DSC_LOG_INFO << "Time to assemble FEM Newton stiffness matrix for current iteration: "
                   << stepAssembleTimer.elapsed() << "s" << std::endl;

      // assemble right hand side
      const typename FEM::DiffusionType diffusion_op;
      rhsassembler.template assemble_for_Newton_method< fem_polorder >(f,
                                                                       diffusion_op,
                                                                       solution,
                                                                       system_rhs);

      const Dune::L2Norm< typename FEM::DiscreteFunctionType::GridPartType > l2norm(system_rhs.gridPart());
      rhs_L2_norm = l2norm.norm(system_rhs);
      if (rhs_L2_norm < 1e-10)
      {
        // residual solution almost identical to zero: break
        DSC_LOG_INFO << "Residual solution almost identical to zero. Therefore: break loop." << std::endl;
        DSC_LOG_INFO << "(L^2-Norm of current right hand side = " << rhs_L2_norm << " < 1e-10)" << std::endl;
        break;
      }
      // set Dirichlet Boundary to zero
      boundaryTreatment(system_rhs);

      const typename FEM::InverseFEMMatrix fem_newton_biCGStab(system_matrix, 1e-8, 1e-8, 20000, true);
      fem_newton_biCGStab(system_rhs, residual);

      if ( residual.dofsValid() )
      {
        solution += residual;
        relative_newton_error_finescale = l2norm.norm(residual);
        relative_newton_error_finescale /= l2norm.norm(solution);

        DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
        // residual solution almost identical to zero: break
        DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
        if (relative_newton_error_finescale <= tolerance)
        {
          DSC_LOG_INFO << "Since tolerance = " << tolerance << ": break loop." << std::endl;
        }
        residual.clear();
      } else {
        DSC_LOG_INFO << "WARNING! Invalid dofs in 'residual'." << std::endl;
        break;
      }
      iteration_step += 1;
    }
    DSC_LOG_INFO << "Problem with FEM + Newton-Method solved in " << assembleTimer.elapsed() << "s." << std::endl
                 << std::endl;
    DSC_LOG_INFO << seperator_line;
    DSC_LOG_INFO << "Problem with FEM + Newton-Method solved in " << assembleTimer.elapsed() << "s." << std::endl
              << std::endl << std::endl;
  }// end 'problem.linear <-> else'

  //! ********************** End of assembling the reference problem ***************************
}


template <class ProblemDataType>
void print_info(ProblemDataType info, std::ostream& out)
{
  // epsilon is specified in the parameter file
  // 'epsilon' in for instance A^{epsilon}(x) = A(x,x/epsilon)
  const double epsilon_ = DSC_CONFIG_GET("problem.epsilon", 1.0f);
  const int refinement_level_ = DSC_CONFIG_GET("fem.grid_level", 4);
  out << "Log-File for Elliptic Model Problem " << Problem::name << "." << std::endl << std::endl;
  if (DSC_CONFIG_GET("problem.linear", true))
    out << "Problem is declared as being LINEAR." << std::endl;
  else
    out << "Problem is declared as being NONLINEAR." << std::endl;

  if (info.has_exact_solution) {
    out << "Exact solution is available." << std::endl << std::endl;
  } else {
    out << "Exact solution is not available." << std::endl << std::endl;
  }
  out << "Computations were made for:" << std::endl << std::endl;
  out << "Refinement Level for Grid = " << refinement_level_ << std::endl << std::endl;

  out << "Epsilon = " << epsilon_ << std::endl << std::endl;
}


//! the main FEM computation
template < class FEMTraits >
void algorithm(typename FEMTraits::GridPointerType& macro_grid_pointer,   // grid pointer that belongs to the macro grid
               const std::string filename) {
  typedef FEMTraits FEM;
  using namespace Dune;

  const typename FEM::ModelProblemDataType problem_data;
  print_info(problem_data, DSC_LOG_INFO);
  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for the finite element problem
  typename FEM::GridPartType gridPart(*macro_grid_pointer);
  //! --------------------------------------------------------------------------------------

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  typename FEM::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  //! --------------------------------------------------------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const typename FEM::DiffusionType diffusion_op;

  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  Dune::RightHandSideAssembler< typename FEM::DiscreteFunctionType > rhsassembler;
  const typename FEM::FirstSourceType f;   // standard source f

  //! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
  const typename FEM::EllipticOperatorType discrete_elliptic_op(discreteFunctionSpace, diffusion_op);

  //! solution vector
  // - By solution, we denote the "discrete solution" determined with FEM in the linear case or FEM-Newton (nonlinear case)
  //    ( if the elliptic problem is linear, the 'solution' is determined without the Newton method )
  // - solution of the finite element method, where we use the Newton method to solve the non-linear system of equations
  typename FEM::DiscreteFunctionType discrete_solution(filename + " FEM(-Newton) Solution", discreteFunctionSpace);
  discrete_solution.clear();

  solve<FEM>(discrete_solution, discreteFunctionSpace, discrete_elliptic_op, filename, rhsassembler);

}

#endif // DUNE_FEM_ALGORITHM_HH
