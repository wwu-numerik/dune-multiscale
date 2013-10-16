#include "algorithm_step.hh"

#include <string>
#include <fstream>

#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/newton_rhs.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/hmm/result.hh>

#include <dune/multiscale/hmm/cell_problem_solver.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/elliptic_hmm_matrix_assembler.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/common/output_traits.hh>

#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "algorithm_error.hh"

namespace {
const std::string seperator_line =
    "---------------------------------------------------------------------------------\n";
}

namespace Dune {
namespace Multiscale {
namespace HMM {

//! set the dirichlet points to zero
struct BoundaryTreatment {
  template <class DiscreteFunctionType>
  static void apply(DiscreteFunctionType& rhs) {
    using namespace Dune::Stuff;
    const auto& discreteFunctionSpace = rhs.space();
    static const unsigned int faceCodim = 1;
    for (const auto& entity : discreteFunctionSpace) {
      for (const auto& intersection : DSC::intersectionRange(discreteFunctionSpace.gridPart(), entity)) {
        if (!intersection.boundary())
          continue;
        if (intersection.boundary() && (intersection.boundaryId() != 1))
          continue;

        auto rhsLocal = rhs.localFunction(entity);
        const auto face = intersection.indexInInside();
        for (auto point : DSC::lagrangePointSetRange<faceCodim>(rhs.space(), entity, face))
          rhsLocal[point] = 0;
      }
    }
  }
}; // boundaryTreatment

//! set the dirichlet points to the Dirichlet BC
template <class DirichletBC, class DiscreteFunctionType>
void setDirichletValues(DirichletBC& dirichlet_func, DiscreteFunctionType& func) {
  using namespace Dune::Stuff;
  const auto& discreteFunctionSpace = func.space();
  static const unsigned int faceCodim = 1;
  for (const auto& entity : discreteFunctionSpace) {
    for (const auto& intersection : DSC::intersectionRange(discreteFunctionSpace.gridPart(), entity)) {
      if (!intersection.boundary())
        continue;
      if (intersection.boundary() && (intersection.boundaryId() != 1))
        continue;

      auto funcLocal = func.localFunction(entity);
      const auto face = intersection.indexInInside();
      for (auto loc_point : DSC::lagrangePointSetRange<faceCodim>(func.space(), entity, face)) {
        const auto& global_point =
            entity.geometry().global(discreteFunctionSpace.lagrangePointSet(entity).point(loc_point));
        CommonTraits::RangeType dirichlet_value(0.0);
        dirichlet_func.evaluate(global_point, dirichlet_value);
        funcLocal[loc_point] = dirichlet_value;
      }
    }
  }
} // setDirichletValues

//! \TODO docme
void solve_hmm_problem_nonlinear(
    const typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
    const typename CommonTraits::DiffusionType& diffusion_op, typename CommonTraits::DiscreteFunctionType& hmm_solution,
    const CellProblemNumberingManager& cp_num_manager,
    const typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace, const int loop_cycle) {
  // the nonlinear case
  // solve cell problems in a preprocess
  //! -------------- solve and save the cell problems for the macroscopic base function set
  const CellProblemSolver cell_problem_solver(periodicDiscreteFunctionSpace, diffusion_op);
  const int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);
  DSC_LOG_INFO << "Start solving cell problems for " << number_of_grid_elements << " leaf entities..." << std::endl;
  // generate directory for cell problem data output
  // only for the case with test function reconstruction (=non-Petrov-Galerkin):
  if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
    // -------------- solve cell problems for the macro basefunction set ------------------------------
    // save the solutions of the cell problems for the set of macroscopic base functions

    cell_problem_solver.saveTheSolutions_baseSet(discreteFunctionSpace, cp_num_manager);

    DSC_LOG_INFO << "Solving the cell problems for the base function set succeeded." << std::endl;
    // end solving and saving cell problems
  }
  //! --------------- end solving and saving cell problems -----------------------------------------

  DSC_LOG_INFO << "Solving nonlinear HMM problem with Newton method." << std::endl;
  DSC_LOG_INFO << seperator_line << std::endl;

  DSC_PROFILER.startTiming("hmm.assemble");

  // just to provide some information
  typename HMMTraits::PeriodicDiscreteFunctionType dummy_periodic_func("a periodic dummy",
                                                                       periodicDiscreteFunctionSpace);
  dummy_periodic_func.clear();

  typename CommonTraits::RangeType relative_newton_error = 10000.0;
  typename CommonTraits::RangeType hmm_rhs_L2_norm = 10000.0;

  // number of HMM Newton step (1 = first step)
  // HMM_NEWTON_ITERATION_STEP' Netwon steps have been already performed,
  // the next one is 'HMM_NEWTON_ITERATION_STEP+1' = hmm_iteration_step
  int hmm_iteration_step = DSC_CONFIG_GET("HMM_NEWTON_ITERATION_STEP", 0) + 1;

  if (DSC_CONFIG_GET("RESUME_TO_BROKEN_COMPUTATION", false)) {
    const int refinement_level_macrogrid_ = DSC_CONFIG_GET("grid.refinement_level_macrogrid", 0);
    // std :: string location_hmm_newton_step_solution = "data/HMM/test/hmm_solution_discFunc_refLevel_5_NewtonStep_2";
    const std::string location_hmm_newton_step_solution =
        (boost::format("hmm_solution_discFunc_refLevel_%d_NewtonStep_%d") % refinement_level_macrogrid_ %
         DSC_CONFIG_GET("HMM_NEWTON_ITERATION_STEP", 0)).str();

    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader_hmm_newton_ref(location_hmm_newton_step_solution);
    discrete_function_reader_hmm_newton_ref.read(0, hmm_solution);
  }

  double old_error = 100.0;
  double error_decay = 0.0;

  // the Newton step for the nonlinear HMM problem:
  // L2-Norm of residual < tolerance ?
  const double hmm_tolerance = DSC_CONFIG_GET("problem.stochastic_pertubation", false)
                                   ? 1e-01 * DSC_CONFIG_GET("problem.stochastic_variance", 0.01)
                                   : 1e-05;

  //! matrix
  typename CommonTraits::LinearOperatorType hmm_newton_matrix("HMM Newton stiffness matrix", discreteFunctionSpace,
                                                              discreteFunctionSpace);
  //! define the elliptic hmm operator that describes our 'homogenized' macro problem
  // ( effect of the elliptic hmm operator on a certain discrete function )
  const DiscreteEllipticHMMOperator discrete_elliptic_hmm_op(discreteFunctionSpace, periodicDiscreteFunctionSpace,
                                                             diffusion_op, cp_num_manager);

  while (relative_newton_error > hmm_tolerance) {
    // (here: hmm_solution = solution from the last iteration step)
    DSC_PROFILER.startTiming("hmm.newton_step", hmm_iteration_step);
    DSC_LOG_INFO << "HMM Newton iteration " << hmm_iteration_step << ":" << std::endl;

    // solve cell problems for the solution of the last iteration step
    cell_problem_solver.saveTheSolutions_discFunc(hmm_solution);
    cell_problem_solver.saveTheJacCorSolutions_baseSet_discFunc(hmm_solution, cp_num_manager);

    // to assemble the computational time
    Dune::Timer stepHmmAssembleTimer;
    // assemble the stiffness matrix
    discrete_elliptic_hmm_op.assemble_jacobian_matrix(hmm_solution, hmm_newton_matrix);

    DSC_LOG_INFO << "Time to assemble HMM stiffness matrix for current Newton iteration: "
                 << stepHmmAssembleTimer.elapsed() << "s" << std::endl;

    DSC_LOG_INFO << "Assemble right hand side..." << std::endl;
    // assemble right hand side

    //! right hand side vector
    // right hand side for the hm finite element method with Newton solver:
    typename CommonTraits::DiscreteFunctionType hmm_newton_rhs("hmm rhs", discreteFunctionSpace);
    hmm_newton_rhs.clear();

    const auto f = Problem::getFirstSource();

    const NewtonRightHandSide newton_rhs = {};
    newton_rhs.assemble_for_HMM_Newton_method(
        *f, diffusion_op, hmm_solution, cp_num_manager, dummy_periodic_func, hmm_newton_rhs);
    DSC_LOG_INFO << "Right hand side assembled!" << std::endl;

    if (!(hmm_newton_rhs.dofsValid())) {
      DSC_LOG_INFO << "Right hand side invalid!" << std::endl;
      DSC_LOG_INFO << "Right hand side invalid!" << std::endl;
      DUNE_THROW(Dune::InvalidStateException, "Right hand side invalid!");
    } else {
      DSC_LOG_INFO << "Right hand side valid ";
    }

    // starting value for the Newton method
    typename CommonTraits::DiscreteFunctionType zero_func_coarse(" constant zero function coarse ",
                                                                 discreteFunctionSpace);
    zero_func_coarse.clear();

    const Dune::Fem::LPNorm<typename CommonTraits::GridPartType> l2norm(hmm_newton_rhs.space().gridPart(), 2);
    hmm_rhs_L2_norm = l2norm.distance(zero_func_coarse, hmm_newton_rhs);

    DSC_LOG_INFO << "with L^2-Norm = " << hmm_rhs_L2_norm << "." << std::endl;
    DSC_LOG_INFO << "Assembled right hand side, with L^2-Norm of RHS = " << hmm_rhs_L2_norm << "." << std::endl;

    if (hmm_rhs_L2_norm < 1e-10) {
      // residual solution almost identical to zero: break
      DSC_LOG_INFO << "HMM residual solution almost identical to zero. Therefore: break loop." << std::endl;
      DSC_LOG_INFO << "(L^2-Norm of current right hand side = " << hmm_rhs_L2_norm << " < 1e-10)" << std::endl;
      break;
    }

    // set Dirichlet Boundary to zero
    BoundaryTreatment::apply(hmm_newton_rhs);

    if (!process_hmm_newton_residual(relative_newton_error, hmm_solution, hmm_newton_matrix, hmm_newton_rhs,
                                     hmm_iteration_step, loop_cycle, hmm_tolerance)) {
      break; // invalid dofs in residual or
    }

    if (relative_newton_error > hmm_tolerance) {
      auto newton_step_time = DSC_PROFILER.stopTiming("hmm.newton_step", hmm_iteration_step) / 1000.f;
      DSC_LOG_INFO << std::endl << "Total time for current HMM Newton step = " << newton_step_time << "s." << std::endl
                   << std::endl;

      error_decay = relative_newton_error / old_error;
      old_error = relative_newton_error;
      // maximum number of Newton iterations
      if ((hmm_iteration_step >= 20) || (error_decay >= 0.95)) {
        DSC_LOG_INFO << std::endl
                     << "Reached constant or inceasing error decay or maximum number of Newton iterations:  break loop."
                     << std::endl;
        DSC_LOG_INFO << "....................................................." << std::endl << std::endl;
        break;
      }
      DSC_LOG_INFO << "....................................................." << std::endl << std::endl;
    }

    hmm_iteration_step += 1;

  } // while( relative_newton_error > hmm_tolerance )

  const auto elapsed = DSC_PROFILER.stopTiming("hmm.assemble");
  DSC_LOG_INFO << seperator_line << "HMM problem with Newton method solved in " << elapsed / 1000.f << "s." << std::endl
               << std::endl;

  if (DSC_CONFIG_GET("hmm.adaptivity", false)) {
    //! TODO which section the local time needs to be added to
    // or if it's necessary at all
    // total_hmm_time += DSC_PROFILER.stopTiming("hmmAssemble");
  }
} // solve_cell_problems_nonlinear

//! \TODO docme
void
solve_hmm_problem_linear(const typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                         const typename CommonTraits::DiffusionType& diffusion_op,
                         typename CommonTraits::DiscreteFunctionType& hmm_solution,
                         const CellProblemNumberingManager& cp_num_manager,
                         const typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace) {
  {
    //! -------------- solve and save the cell problems for the base function set --------------------------------------
    const CellProblemSolver cell_problem_solver(periodicDiscreteFunctionSpace, diffusion_op);
    const int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);
    DSC_LOG_INFO << "Solving cell problems for " << number_of_grid_elements << " leaf entities." << std::endl;
    // -------------- solve cell problems for the macro basefunction set ------------------------------
    // save the solutions of the cell problems for the set of macroscopic base functions
    cell_problem_solver.saveTheSolutions_baseSet(discreteFunctionSpace, cp_num_manager);
    // ------------- end solving and saving cell problems for the macro basefunction set --------------
    //! --------------- end solving and saving cell problems -----------------------------------------
  }

  DSC_LOG_INFO << "Solving linear HMM problem." << std::endl;
  DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

  // to assemble the computational time
  Dune::Timer hmmAssembleTimer;

  //! matrix
  typename CommonTraits::LinearOperatorType hmm_matrix("HMM Newton stiffness matrix", discreteFunctionSpace,
                                                       discreteFunctionSpace);
  const DiscreteEllipticHMMOperator discrete_elliptic_hmm_op(discreteFunctionSpace, periodicDiscreteFunctionSpace,
                                                             diffusion_op, cp_num_manager);
  // assemble the hmm stiffness matrix
  discrete_elliptic_hmm_op.assemble_matrix(hmm_matrix);
  // to print the matrix, use:   hmm_matrix.print();

  DSC_LOG_INFO << "Time to assemble HMM macro stiffness matrix: " << hmmAssembleTimer.elapsed() << "s" << std::endl;

  // assemble right hand side
  //! right hand side vector
  // if we have some additional source term (-div G), define:
  const auto G = Problem::getSecondSource();
  const auto f = Problem::getFirstSource(); // standard source f
                                            // lower order term F(x, u(x), grad u(x) ):

  // Dirichlet boundary condition
  const auto dirichlet_bc = Problem::getDirichletBC();
  // Neumann boundary condition
  const auto neumann_bc = Problem::getNeumannBC();

  // - div ( A^{\epsilon} \nabla u^{\epsilon} ) + F(x, u^{\epsilon}, \nabla u^{\epsilon}) = f - div G

  // discrete function that takes the values of the Dirichlet BC on the Dirichlet Boundary nodes
  // and that is zero elsewhere
  typename CommonTraits::DiscreteFunctionType dirichlet_extension("Dirichlet extension", discreteFunctionSpace);
  dirichlet_extension.clear();
  setDirichletValues(*dirichlet_bc, dirichlet_extension);

  // right hand side for the hm finite element method with Newton solver:
  typename CommonTraits::DiscreteFunctionType hmm_rhs("hmm rhs", discreteFunctionSpace);
  hmm_rhs.clear();
  const RightHandSideAssembler<typename CommonTraits::DiscreteFunctionType> rhsassembler = {};
  rhsassembler.assemble<2 * CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2>(
      *f, diffusion_op, dirichlet_extension, *neumann_bc, hmm_rhs);

  // set Dirichlet Boundary to zero
  BoundaryTreatment::apply(hmm_rhs);

  typename HMMTraits::InverseOperatorType hmm_biCGStab(hmm_matrix, 1e-8, 1e-8, 20000,
                                                       DSC_CONFIG_GET("global.cgsolver_verbose", false));
  hmm_biCGStab(hmm_rhs, hmm_solution);
  DSC_LOG_INFO << seperator_line << "Linear HMM problem solved in " << hmmAssembleTimer.elapsed() << "s." << std::endl
               << std::endl;

  hmm_solution += dirichlet_extension;
}

//! \TODO docme
bool process_hmm_newton_residual(typename CommonTraits::RangeType& relative_newton_error,
                                 typename CommonTraits::DiscreteFunctionType& hmm_solution,
                                 typename CommonTraits::LinearOperatorType& hmm_newton_matrix,
                                 const typename CommonTraits::DiscreteFunctionType& hmm_newton_rhs,
                                 const int hmm_iteration_step, const int loop_cycle, const double hmm_tolerance) {
  //! residual vector
  // current residual
  typename CommonTraits::DiscreteFunctionType hmm_newton_residual("HMM Newton Residual", hmm_solution.space());
  hmm_newton_residual.clear();
  const int refinement_level_macrogrid_ = DSC_CONFIG_GET("grid.refinement_level_macrogrid", 0);

  double hmm_biCG_tolerance = 1e-8;
  bool hmm_solution_convenient = false;
  while (!hmm_solution_convenient) {
    hmm_newton_residual.clear();
    const typename HMMTraits::InverseOperatorType hmm_newton_biCGStab(hmm_newton_matrix, 1e-8, hmm_biCG_tolerance,
                                                                      20000, false);

    hmm_newton_biCGStab(hmm_newton_rhs, hmm_newton_residual);
    hmm_solution_convenient = hmm_newton_residual.dofsValid();

    if (hmm_biCG_tolerance > 1e-4) {
      DSC_LOG_INFO << "WARNING! Iteration step " << hmm_iteration_step << ". Invalid dofs in 'hmm_newton_residual'."
                   << std::endl;
      DUNE_THROW(Dune::InvalidStateException, "Right hand side invalid!");
    }
    hmm_biCG_tolerance *= 10.0;
  }

  if (!hmm_newton_residual.dofsValid()) {
    DSC_LOG_ERROR << "WARNING! Invalid dofs in 'hmm_newton_residual'." << std::endl;
    return false;
  }
  hmm_solution += hmm_newton_residual;

  // write the solution after the current HMM Newton step to a file
  // for adaptive computations, the saved solution is not suitable for a later usage
  if (DSC_CONFIG_GET("hmm.adaptivity", false) && DSC_CONFIG_GET("WRITE_HMM_SOL_TO_FILE", true)) {
    std::string fname = (boost::format("hmm_solution_discFunc_refLevel_%d_NewtonStep_%d") %
                         refinement_level_macrogrid_ % hmm_iteration_step).str();
    DiscreteFunctionWriter(fname).append(hmm_solution);

    // writing paraview data output
    // general output parameters
    Dune::Multiscale::OutputParameters outputparam;

    // create and initialize output class
    typename OutputTraits::IOTupleType hmm_solution_newton_step_series(&hmm_solution);
    outputparam.set_prefix((boost::format("hmm_solution_%d_NewtonStep_%d") % loop_cycle % hmm_iteration_step).str());
    typename OutputTraits::DataOutputType hmmsol_dataoutput(hmm_solution.space().gridPart().grid(),
                                                            hmm_solution_newton_step_series, outputparam);
    // write data
    hmmsol_dataoutput.writeData(1.0 /*dummy*/, "hmm-solution-NewtonStep");
  }

  // || u^(n+1) - u^(n) ||_L2
  const Dune::Fem::L2Norm<typename CommonTraits::DiscreteFunctionType::GridPartType> l2norm(
      hmm_newton_residual.gridPart());
  relative_newton_error = l2norm.norm(hmm_newton_residual);
  // || u^(n+1) - u^(n) ||_L2 / || u^(n+1) ||_L2
  relative_newton_error = relative_newton_error / l2norm.norm(hmm_solution);

  DSC_LOG_INFO << "Relative L2 HMM Newton iteration error = " << relative_newton_error << std::endl;

  // residual solution almost identical to zero: break
  if (relative_newton_error <= hmm_tolerance) {
    const auto newton_step_time = DSC_PROFILER.stopTiming("hmm.newton_step", hmm_iteration_step);
    DSC_LOG_INFO << std::endl << "Total time for current HMM Newton step = " << newton_step_time << "ms." << std::endl
                 << std::endl;
    DSC_LOG_INFO << "Since HMM-tolerance = " << hmm_tolerance << ": break loop." << std::endl;
    DSC_LOG_INFO << "....................................................." << std::endl << std::endl;
    return false;
  }
  return true;
}

//! \TODO docme
void step_data_output(const typename CommonTraits::GridPartType& gridPart,
                      const typename CommonTraits::GridPartType& gridPartFine,
                      typename CommonTraits::DiscreteFunctionType& hmm_solution, const int loop_cycle) {
  //! --------------- writing data output ---------------------
  // general output parameters
  Dune::Multiscale::OutputParameters outputparam;

  // --------- data output hmm solution --------------

  // create and initialize output class
  typename OutputTraits::IOTupleType hmm_solution_series(&hmm_solution);
  if (DSC_CONFIG_GET("hmm.adaptivity", false)) {
    outputparam.set_prefix((boost::format("hmm_solution_%d") % loop_cycle).str());
  } else {
    outputparam.set_prefix("hmm_solution");
  }
  typename OutputTraits::DataOutputType hmmsol_dataoutput(gridPart.grid(), hmm_solution_series, outputparam);

  // write data
  hmmsol_dataoutput.writeData(1.0 /*dummy*/, "hmm-solution");
  // -------------------------------------------------------

  if (Problem::getModelData()->hasExactSolution()) {
    // --------- data output discrete exact solution --------------

    // create and initialize output class
    const auto u = Problem::getExactSolution();
    const OutputTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", *u, gridPartFine);
    typename OutputTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("exact_solution");
    typename OutputTraits::ExSolDataOutputType exactsol_dataoutput(gridPartFine.grid(), exact_solution_series,
                                                                   outputparam);
    std::stringstream outstring;
    outstring << "exact-solution";
    exactsol_dataoutput.writeData(1.0 /*dummy*/, outstring.str());
    // clear the std::stringstream:
    outstring.str(std::string());
    // -------------------------------------------------------
  }
  //!-------------------------------------------------------------
}

//! \TODO docme
HMMResult single_step(typename CommonTraits::GridPartType& gridPart, typename CommonTraits::GridPartType& gridPartFine,
                      typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
                      typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                      const typename CommonTraits::DiffusionType& diffusion_op,
                      typename CommonTraits::DiscreteFunctionType& hmm_solution,
                      const typename CommonTraits::DiscreteFunctionType& reference_solution, const int loop_cycle) {
  DSC_LOG_INFO << std::endl << "Solving HMM-macro-problem for " << discreteFunctionSpace.size()
               << " unkowns and polynomial order " << CommonTraits::DiscreteFunctionSpaceType::polynomialOrder << "."
               << std::endl << std::endl;

  const Dune::Fem::LPNorm<typename CommonTraits::DiscreteFunctionType::GridPartType> l2norm(gridPart, 2);

  // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain cell problem
  CellProblemNumberingManager cp_num_manager(discreteFunctionSpace);

  if (DSC_CONFIG_GET("problem.linear", true))
    solve_hmm_problem_linear(periodicDiscreteFunctionSpace, diffusion_op, hmm_solution, cp_num_manager,
                             discreteFunctionSpace);
  else // for a given loop cycle of the Newton scheme:
    solve_hmm_problem_nonlinear(periodicDiscreteFunctionSpace, diffusion_op, hmm_solution, cp_num_manager,
                                discreteFunctionSpace, loop_cycle);

  const auto errors = estimate_error(gridPart, discreteFunctionSpace, periodicDiscreteFunctionSpace, diffusion_op,
                                     cp_num_manager, hmm_solution);

  const int refinement_level_macrogrid_ = DSC_CONFIG_GET("hmm.coarse_grid_level", 4);
  // for adaptive computations, the saved solution is not suitable for a later usage
  if (!DSC_CONFIG_GET("hmm.adaptivity", false) && DSC_CONFIG_GET("hmm.write_to_file", false)) {
    DiscreteFunctionWriter((boost::format("hmm_solution_discFunc_refLevel_%d") % refinement_level_macrogrid_).str())
        .append(hmm_solution);
  }

  DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;

  //! ----------------- compute L2-errors -------------------

  // L2 error with reference solution
  if (DSC_CONFIG_GET("problem.reference_solution", false)) {
    DSC_PROFILER.startTiming("hmm.timeadapt");

    // project HMM solution in finer discrete function space
    typename CommonTraits::DiscreteFunctionType projected_hmm_solution("HMM Solution Projection",
                                                                       reference_solution.space());
    projected_hmm_solution.clear();
    Dune::Stuff::HeterogenousProjection<DSG::EntityInlevelSearch>::project(hmm_solution /*source*/,
                                                                           projected_hmm_solution /*target*/);

    const auto hmm_error = l2norm.distance(projected_hmm_solution, reference_solution);

    // old (expensive) hack to deal with discrete functions, defined on different grids
    // (should do the same as the heterogenous projection above - could therefore be used for comparison)
    /*
    static const int hmm_polorder = 2* CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2;
    Dune::ImprovedL2Error< typename CommonTraits::DiscreteFunctionType > impL2error;
    typename CommonTraits::RangeType hmm_error = impL2error.template norm_adaptive_grids_2< hmm_polorder >(
      hmm_solution,
      reference_solution);

    const auto timeadapt = DSC_PROFILER.stopTiming("hmm.timeadapt") / 1000.f;
    // if it took longer then 1 minute to compute the error:
    if (timeadapt > 60)
    {
      DSC_LOG_INFO << "WARNING! EXPENSIVE! Error assembled in " << timeadapt << "s." << std::endl << std::endl;
    }
    */

    DSC_LOG_INFO << "|| u_hmm - u_reference ||_L2 =  " << hmm_error << std::endl << std::endl;
  }

  // L2 errors with exact solution
  if (Problem::getModelData()->hasExactSolution()) {
    const auto u = Problem::getExactSolution();
    OutputTraits::DiscreteExactSolutionType gf("Exact Solution", *u, hmm_solution.space().gridPart());

    const auto exact_hmm_error = l2norm.distance(gf, hmm_solution);

    DSC_LOG_INFO << "|| u_hmm - u_exact ||_L2 =  " << exact_hmm_error << std::endl << std::endl;
    if (DSC_CONFIG_GET("problem.reference_solution", false)) {
      const auto reference_sol_error = l2norm.distance(gf, reference_solution);

      DSC_LOG_INFO << "|| u_reference - u_exact ||_L2 =  " << reference_sol_error << std::endl << std::endl;
    }
  }

  if (DSC_CONFIG_GET("hmm.error_estimation", false)) {
    DSC_LOG_INFO << "Estimated error = " << errors.estimated_error << "." << std::endl;
    DSC_LOG_INFO << "In detail:" << std::endl;
    DSC_LOG_INFO << "   Estimated source error = " << errors.estimated_source_error << "." << std::endl;
    DSC_LOG_INFO << "   Estimated approximation error = " << errors.estimated_approximation_error << "." << std::endl;
    DSC_LOG_INFO << "   Estimated residual error = " << errors.estimated_residual_error << ", where:" << std::endl;
    DSC_LOG_INFO << "        contribution of macro jumps = " << errors.estimated_residual_error_macro_jumps << " and "
                 << std::endl;
    DSC_LOG_INFO << "        contribution of micro jumps = " << errors.estimated_residual_error_micro_jumps;
    if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
      DSC_LOG_INFO << " and " << std::endl;
      DSC_LOG_INFO << "   Estimated tfr error = " << errors.estimated_tfr_error << "." << std::endl;
    }
    DSC_LOG_INFO << std::endl;
  }
  //! -------------------------------------------------------

  step_data_output(gridPart, gridPartFine, hmm_solution, loop_cycle);
  return errors;
}

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {
