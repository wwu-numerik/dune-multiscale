#ifndef ALGORITHM_STEP_HH
#define ALGORITHM_STEP_HH

#include "algorithm_error.hh"


#include <string>
#include <fstream>

#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/tools/improved_l2error.hh>
#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/solver.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>

namespace {
  const std::string seperator_line = "---------------------------------------------------------------------------------\n";
}

/**
 * \return true if a program should continue with a new newton step
 **/
template < class HMM >
bool process_hmm_newton_residual(typename HMM::RangeType& relative_newton_error,
                                 typename HMM::DiscreteFunctionType& hmm_solution,
                                 const typename HMM::FEMMatrix& hmm_newton_matrix,
                                 const typename HMM::DiscreteFunctionType& hmm_newton_rhs,
                                 const int hmm_iteration_step,
                                 const int loop_cycle,
                                 const double hmm_tolerance );

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

template <class HMM>
void solve_cell_problems_nonlinear(const typename HMM::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                                   const typename HMM::DiffusionType& diffusion_op,
                                   typename HMM::DiscreteFunctionType& hmm_solution,
                                   const typename HMM::CellProblemNumberingManagerType& cp_num_manager,
                                   const typename HMM::DiscreteFunctionSpaceType& discreteFunctionSpace,
                                   const Dune::RightHandSideAssembler< typename HMM::DiscreteFunctionType >& rhsassembler,
                                   const int loop_cycle) {
    // the nonlinear case
    // solve cell problems in a preprocess

    //! -------------- solve and save the cell problems for the macroscopic base function set
    const Dune::CellProblemSolver< typename HMM::PeriodicDiscreteFunctionType,
        typename HMM::DiffusionType> cell_problem_solver(periodicDiscreteFunctionSpace, diffusion_op );
    const int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);
    DSC_LOG_INFO << "Start solving cell problems for " << number_of_grid_elements << " leaf entities..." << std::endl;
      // generate directory for cell problem data output
    // only for the case with test function reconstruction (=non-Petrov-Galerkin):
    if ( !DSC_CONFIG_GET("hmm.petrov_galerkin", true ) ) {
      // -------------- solve cell problems for the macro basefunction set ------------------------------
      // save the solutions of the cell problems for the set of macroscopic base functions

      cell_problem_solver.template saveTheSolutions_baseSet< typename HMM::DiscreteFunctionType >(discreteFunctionSpace,
                                                                           cp_num_manager);

      DSC_LOG_INFO << "Solving the cell problems for the base function set succeeded." << std::endl;
      // end solving and saving cell problems
    }
    //! --------------- end solving and saving cell problems -----------------------------------------

    DSC_LOG_INFO << "Solving nonlinear HMM problem with Newton method." << std::endl;
    DSC_LOG_INFO << seperator_line << std::endl;

    DSC_PROFILER.startTiming("hmmAssemble");

    // just to provide some information
    typename HMM::PeriodicDiscreteFunctionType dummy_periodic_func("a periodic dummy", periodicDiscreteFunctionSpace);
    dummy_periodic_func.clear();


    typename HMM::RangeType relative_newton_error = 10000.0;
    typename HMM::RangeType hmm_rhs_L2_norm = 10000.0;

    // number of HMM Newton step (1 = first step)
    // HMM_NEWTON_ITERATION_STEP' Netwon steps have been already performed,
    // the next one is 'HMM_NEWTON_ITERATION_STEP+1' = hmm_iteration_step
    int hmm_iteration_step = DSC_CONFIG_GET("HMM_NEWTON_ITERATION_STEP", 0) + 1;

    if (DSC_CONFIG_GET("RESUME_TO_BROKEN_COMPUTATION", false)) {
      const int refinement_level_macrogrid_ = DSC_CONFIG_GET("grid.refinement_level_macrogrid", 0);
      // std :: string location_hmm_newton_step_solution = "data/HMM/test/hmm_solution_discFunc_refLevel_5_NewtonStep_2";
      const std::string location_hmm_newton_step_solution =
          (boost::format("hmm_solution_discFunc_refLevel_%d_NewtonStep_%d")
           % refinement_level_macrogrid_ % DSC_CONFIG_GET("HMM_NEWTON_ITERATION_STEP", 0)).str();

      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader_hmm_newton_ref(location_hmm_newton_step_solution);
      discrete_function_reader_hmm_newton_ref.read(0, hmm_solution);
    }

    double old_error = 100.0;
    double error_decay = 0.0;

    // the Newton step for the nonlinear HMM problem:
    // L2-Norm of residual < tolerance ?
    const double hmm_tolerance = DSC_CONFIG_GET("problem.stochastic_pertubation", false)
                                  ? 1e-01 * DSC_CONFIG_GET("problem.stochastic_variance",  0.01)
                                  : 1e-05;

    //! matrix
    typename HMM::FEMMatrix hmm_newton_matrix("HMM Newton stiffness matrix", discreteFunctionSpace, discreteFunctionSpace);
    //! define the elliptic hmm operator that describes our 'homogenized' macro problem
    // ( effect of the elliptic hmm operator on a certain discrete function )
    const typename HMM::EllipticHMMOperatorType discrete_elliptic_hmm_op(discreteFunctionSpace,
                                                     periodicDiscreteFunctionSpace,
                                                     diffusion_op,
                                                     cp_num_manager);

    while (relative_newton_error > hmm_tolerance)
    {
      // (here: hmm_solution = solution from the last iteration step)
      DSC_PROFILER.startTiming("newton_step_", hmm_iteration_step);
      DSC_LOG_INFO << "HMM Newton iteration " << hmm_iteration_step << ":" << std::endl;

      // solve cell problems for the solution of the last iteration step
      cell_problem_solver.template saveTheSolutions_discFunc< typename HMM::DiscreteFunctionType >(hmm_solution);
      cell_problem_solver.template saveTheJacCorSolutions_baseSet_discFunc< typename HMM::DiscreteFunctionType >(hmm_solution,
                                                                                          cp_num_manager);

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
      typename HMM::DiscreteFunctionType hmm_newton_rhs("hmm rhs", discreteFunctionSpace);
      hmm_newton_rhs.clear();

      const typename HMM::FirstSourceType f;   // standard source f
      rhsassembler.template assemble_for_HMM_Newton_method< HMM::assembler_order >(
        f,
        diffusion_op,
        hmm_solution,
        cp_num_manager,
        dummy_periodic_func,
        hmm_newton_rhs);
      DSC_LOG_INFO << "Right hand side assembled!" << std::endl;

      if ( !( hmm_newton_rhs.dofsValid() ) )
      {
        DSC_LOG_INFO << "Right hand side invalid!" << std::endl;
        DSC_LOG_INFO << "Right hand side invalid!" << std::endl;
        DUNE_THROW(Dune::InvalidStateException, "Right hand side invalid!");
      } else {
        DSC_LOG_INFO << "Right hand side valid ";
      }

      // starting value for the Newton method
      typename HMM::DiscreteFunctionType zero_func_coarse( " constant zero function coarse ", discreteFunctionSpace);
      zero_func_coarse.clear();

      const Dune::L2Error< typename HMM::DiscreteFunctionType > l2error;
      hmm_rhs_L2_norm = l2error.template norm2< HMM::assembler_order >(zero_func_coarse, hmm_newton_rhs);

      DSC_LOG_INFO << "with L^2-Norm = " << hmm_rhs_L2_norm << "." << std::endl;
      DSC_LOG_INFO << "Assembled right hand side, with L^2-Norm of RHS = " << hmm_rhs_L2_norm << "." << std::endl;

      if (hmm_rhs_L2_norm < 1e-10)
      {
        // residual solution almost identical to zero: break
        DSC_LOG_INFO << "HMM residual solution almost identical to zero. Therefore: break loop." << std::endl;
        DSC_LOG_INFO << "(L^2-Norm of current right hand side = " << hmm_rhs_L2_norm << " < 1e-10)" << std::endl;
        break;
      }

      // set Dirichlet Boundary to zero
      boundaryTreatment(hmm_newton_rhs);

      if (!process_hmm_newton_residual<HMM>(relative_newton_error, hmm_solution, hmm_newton_matrix, hmm_newton_rhs,
                                  hmm_iteration_step, loop_cycle, hmm_tolerance)) {
        break;//invalid dofs in residual or
      }

      if (relative_newton_error > hmm_tolerance)
      {
        auto newton_step_time = DSC_PROFILER.stopTiming("newton_step_", hmm_iteration_step) / 1000.f;
        DSC_LOG_INFO << std::endl << "Total time for current HMM Newton step = " << newton_step_time << "s."
                    << std::endl << std::endl;

        error_decay = relative_newton_error / old_error;
        old_error = relative_newton_error;
        // maximum number of Newton iterations
        if ( (hmm_iteration_step >= 20) || (error_decay >= 0.95) )
        {
          DSC_LOG_INFO << std::endl
                    << "Reached constant or inceasing error decay or maximum number of Newton iterations:  break loop."
                    << std::endl;
          DSC_LOG_INFO << "....................................................." << std::endl << std::endl;
          break;
        }
        DSC_LOG_INFO << "....................................................." << std::endl << std::endl;
      }

      hmm_iteration_step += 1;

    }       // while( relative_newton_error > hmm_tolerance )

    const auto elapsed = DSC_PROFILER.stopTiming("hmmAssemble");
    DSC_LOG_INFO << seperator_line << "HMM problem with Newton method solved in " << elapsed / 1000.f << "s." << std::endl
              << std::endl;

    if (DSC_CONFIG_GET("hmm.adaptivity", false)) {
      //! TODO which section the local time needs to be added to
      // or if it's necessary at all
      //total_hmm_time += DSC_PROFILER.stopTiming("hmmAssemble");
    }
} //solve_cell_problems_nonlinear

template <class HMM>
void solve_cell_problems_linear(const typename HMM::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                                const typename HMM::DiffusionType& diffusion_op,
                                typename HMM::DiscreteFunctionType& hmm_solution,
                                const typename HMM::CellProblemNumberingManagerType& cp_num_manager,
                                const typename HMM::DiscreteFunctionSpaceType& discreteFunctionSpace,
                                const Dune::RightHandSideAssembler< typename HMM::DiscreteFunctionType >& rhsassembler) {

  //! -------------- solve and save the cell problems for the base function set --------------------------------------
  const Dune::CellProblemSolver< typename HMM::PeriodicDiscreteFunctionType,
      typename HMM::DiffusionType> cell_problem_solver(periodicDiscreteFunctionSpace, diffusion_op );
  const int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);
  DSC_LOG_INFO << "Solving cell problems for " << number_of_grid_elements << " leaf entities." << std::endl;
  // -------------- solve cell problems for the macro basefunction set ------------------------------
  // save the solutions of the cell problems for the set of macroscopic base functions
  cell_problem_solver.template saveTheSolutions_baseSet< typename HMM::DiscreteFunctionType >(discreteFunctionSpace,
                                                                       cp_num_manager);
  // ------------- end solving and saving cell problems for the macro basefunction set --------------
  //! --------------- end solving and saving cell problems -----------------------------------------

  DSC_LOG_INFO << "Solving linear HMM problem." << std::endl;
  DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

  // to assemble the computational time
  Dune::Timer hmmAssembleTimer;

  //! matrix
  typename HMM::FEMMatrix hmm_newton_matrix("HMM Newton stiffness matrix", discreteFunctionSpace, discreteFunctionSpace);
  const typename HMM::EllipticHMMOperatorType discrete_elliptic_hmm_op(discreteFunctionSpace,
                                                   periodicDiscreteFunctionSpace,
                                                   diffusion_op,
                                                   cp_num_manager);
  // assemble the hmm stiffness matrix
  discrete_elliptic_hmm_op.assemble_matrix(hmm_newton_matrix);
  // to print the matrix, use:   hmm_newton_matrix.print();

  DSC_LOG_INFO << "Time to assemble HMM macro stiffness matrix: " << hmmAssembleTimer.elapsed() << "s" << std::endl;

  // assemble right hand side
  //! right hand side vector
  const typename HMM::FirstSourceType f;   // standard source f
  // right hand side for the hm finite element method with Newton solver:
  typename HMM::DiscreteFunctionType hmm_newton_rhs("hmm rhs", discreteFunctionSpace);
  hmm_newton_rhs.clear();
  rhsassembler.template assemble< 2* HMM::DiscreteFunctionSpaceType::polynomialOrder + 2 >(f, hmm_newton_rhs);

  // set Dirichlet Boundary to zero
  boundaryTreatment(hmm_newton_rhs);

  typename HMM::InverseFEMMatrix hmm_biCGStab(hmm_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE);
  hmm_biCGStab(hmm_newton_rhs, hmm_solution);
  DSC_LOG_INFO << seperator_line << "Linear HMM problem solved in " << hmmAssembleTimer.elapsed() << "s." << std::endl << std::endl;
}

template < class HMM >
bool process_hmm_newton_residual(typename HMM::RangeType& relative_newton_error,
                                 typename HMM::DiscreteFunctionType& hmm_solution,
                                 const typename HMM::FEMMatrix& hmm_newton_matrix,
                                 const typename HMM::DiscreteFunctionType& hmm_newton_rhs,
                                 const int hmm_iteration_step,
                                 const int loop_cycle,
                                 const double hmm_tolerance ) {
  //! residual vector
  // current residual
  typename HMM::DiscreteFunctionType hmm_newton_residual("HMM Newton Residual", hmm_solution.space());
  hmm_newton_residual.clear();
  const int refinement_level_macrogrid_ = DSC_CONFIG_GET("grid.refinement_level_macrogrid", 0);

  double hmm_biCG_tolerance = 1e-8;
  bool hmm_solution_convenient = false;
  while (!hmm_solution_convenient)
  {
    hmm_newton_residual.clear();
    const typename HMM::InverseFEMMatrix hmm_newton_biCGStab(hmm_newton_matrix,
                                           1e-8, hmm_biCG_tolerance, 20000, VERBOSE);

    hmm_newton_biCGStab(hmm_newton_rhs, hmm_newton_residual);

    if ( hmm_newton_residual.dofsValid() )
    {
      hmm_solution_convenient = true;
    }

    if (hmm_biCG_tolerance > 1e-4)
    {
      DSC_LOG_INFO << "WARNING! Iteration step " << hmm_iteration_step << ". Invalid dofs in 'hmm_newton_residual'."
                << std::endl;
      DUNE_THROW(Dune::InvalidStateException, "Right hand side invalid!");
    }
    hmm_biCG_tolerance *= 10.0;
  }



  if ( !hmm_newton_residual.dofsValid() ) {
    DSC_LOG_ERROR << "WARNING! Invalid dofs in 'hmm_newton_residual'." << std::endl;
    return false;
  }
  hmm_solution += hmm_newton_residual;

  // write the solution after the current HMM Newton step to a file
  // for adaptive computations, the saved solution is not suitable for a later usage
  if (DSC_CONFIG_GET("hmm.adaptivity", false) && DSC_CONFIG_GET("WRITE_HMM_SOL_TO_FILE", true)) {
    std::string fname = (boost::format("/hmm_solution_discFunc_refLevel_%d_NewtonStep_%d")
                         % refinement_level_macrogrid_ % hmm_iteration_step).str();
    DiscreteFunctionWriter(fname).append(hmm_solution);

    // writing paraview data output
    // general output parameters
    Dune::myDataOutputParameters outputparam;

    // create and initialize output class
    typename HMM::IOTupleType hmm_solution_newton_step_series(&hmm_solution);
    outputparam.set_prefix((boost::format("hmm_solution_%d_NewtonStep_%d") % loop_cycle % hmm_iteration_step).str());
    typename HMM::DataOutputType hmmsol_dataoutput(hmm_solution.space().gridPart().grid(),
                                                   hmm_solution_newton_step_series, outputparam);
    // write data
    hmmsol_dataoutput.writeData( 1.0 /*dummy*/, "hmm-solution-NewtonStep" );
  }

  // || u^(n+1) - u^(n) ||_L2
  const Dune::L2Norm< typename HMM::DiscreteFunctionType::GridPartType > l2norm(hmm_newton_residual.gridPart());
  relative_newton_error = l2norm.norm(hmm_newton_residual);
  // || u^(n+1) - u^(n) ||_L2 / || u^(n+1) ||_L2
  relative_newton_error = relative_newton_error / l2norm.norm(hmm_solution);

  DSC_LOG_INFO << "Relative L2 HMM Newton iteration error = " << relative_newton_error << std::endl;

  // residual solution almost identical to zero: break
  if (relative_newton_error <= hmm_tolerance)
  {
      const auto newton_step_time = DSC_PROFILER.stopTiming("newton_step_", hmm_iteration_step);
      DSC_LOG_INFO << std::endl << "Total time for current HMM Newton step = " << newton_step_time << "ms."
                << std::endl << std::endl;
      DSC_LOG_INFO << "Since HMM-tolerance = " << hmm_tolerance << ": break loop." << std::endl;
      DSC_LOG_INFO << "....................................................." << std::endl << std::endl;
      return false;
  }
  return true;
}

template <class HMM>
void step_data_output(const typename HMM::GridPartType& gridPart,
                      const typename HMM::GridPartType& gridPartFine,
                      typename HMM::DiscreteFunctionType& hmm_solution,
                      const int loop_cycle) {
  //! --------------- writing data output ---------------------
  // general output parameters
  Dune::myDataOutputParameters outputparam;

  // sequence stamp
  std::stringstream outstring;

  // --------- data output hmm solution --------------

  // create and initialize output class
  typename HMM::IOTupleType hmm_solution_series(&hmm_solution);
  if (DSC_CONFIG_GET("hmm.adaptivity", false)) {
    outputparam.set_prefix((boost::format("hmm_solution_%d") % loop_cycle).str());
  }
  else {
    outputparam.set_prefix("hmm_solution");
  }
  typename HMM::DataOutputType hmmsol_dataoutput(gridPart.grid(), hmm_solution_series, outputparam);

  // write data
  hmmsol_dataoutput.writeData( 1.0 /*dummy*/, "hmm-solution" );
  // -------------------------------------------------------

  if (HMM::ModelProblemDataType::has_exact_solution) {
    // --------- data output discrete exact solution --------------

    // create and initialize output class
    typename HMM::ExactSolutionType u;
    typename HMM::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPartFine);
    typename HMM::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("exact_solution");
    typename HMM::ExSolDataOutputType exactsol_dataoutput(gridPartFine.grid(), exact_solution_series, outputparam);

    // write data
    outstring << "exact-solution";
    exactsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
    // clear the std::stringstream:
    outstring.str( std::string() );
    // -------------------------------------------------------
  }

  #ifdef WRITE_FINESCALE_SOL_TO_FILE
  // --------- data output reference solution (fine fem newton computation) --------------

  // create and initialize output class
  typename HMM::IOTupleType fem_newton_solution_series(&fem_newton_solution);
  outputparam.set_prefix("reference_solution");
  typename HMM::DataOutputType fem_newton_dataoutput(gridPartFine.grid(), fem_newton_solution_series, outputparam);

  // write data
  outstring << "reference_solution";
  fem_newton_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str( std::string() );

  // -------------------------------------------------------
  #endif // ifdef WRITE_FINESCALE_SOL_TO_FILE

  #ifdef HOMOGENIZEDSOL_AVAILABLE
  // --------- data output homogenized solution --------------

  // create and initialize output class
  IOTupleType homogenized_solution_series(&homogenized_solution);
  outputparam.set_prefix("homogenized_solution");
  DataOutputType homogenized_solution_dataoutput(gridPartFine.grid(), homogenized_solution_series, outputparam);

  // write data
  outstring << "homogenized_solution";
  homogenized_solution_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str( std::string() );

  // -------------------------------------------------------
  #endif // ifdef HOMOGENIZEDSOL_AVAILABLE
//!-------------------------------------------------------------
}

template < class HMMTraits >
HMMResult<HMMTraits> single_step( typename HMMTraits::GridPartType& gridPart,
        typename HMMTraits::GridPartType& gridPartFine,
        typename HMMTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
        typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
        const typename HMMTraits::DiffusionType& diffusion_op,
        const Dune::RightHandSideAssembler< typename HMMTraits::DiscreteFunctionType >& rhsassembler,
        typename HMMTraits::DiscreteFunctionType& hmm_solution,
        const typename HMMTraits::DiscreteFunctionType& fem_newton_solution,
        const int loop_cycle ) {
    typedef HMMTraits HMM;
    DSC_LOG_INFO << std::endl << "Solving HMM-macro-problem for " << discreteFunctionSpace.size()
              << " unkowns and polynomial order "
              << HMM::DiscreteFunctionSpaceType::polynomialOrder << "."
              << std::endl << std::endl;

    // if we have some additional source term (-div G), define:
    const typename HMM::SecondSourceType G;
    // - div ( A^{\epsilon} \nabla u^{\epsilon} ) = f - div G

    const Dune::L2Error< typename HMM::DiscreteFunctionType > l2error;
    // ----------------------------------------------------------------------------------------------//
    // ----------------------- THE DISCRETE HMM OPERATOR -----------------------------------//
    // ----------------------------------------------------------------------------------------------//

    // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain cell problem
    typename HMM::CellProblemNumberingManagerType cp_num_manager(discreteFunctionSpace);


    // ----------------------------------------------------------------------------------------------//
    // ----------------------------------------------------------------------------------------------//
    // ----------------------------------------------------------------------------------------------//

    if (DSC_CONFIG_GET("problem.linear", true))
      solve_cell_problems_linear<HMM>(periodicDiscreteFunctionSpace, diffusion_op, hmm_solution,
                                    cp_num_manager, discreteFunctionSpace, rhsassembler);
    else
      solve_cell_problems_nonlinear<HMM>(periodicDiscreteFunctionSpace, diffusion_op, hmm_solution,
                                    cp_num_manager, discreteFunctionSpace, rhsassembler, loop_cycle);

    const auto errors = estimate_error<HMM>(gridPart, discreteFunctionSpace, periodicDiscreteFunctionSpace,
                   diffusion_op, cp_num_manager, hmm_solution);

    const int refinement_level_macrogrid_ = DSC_CONFIG_GET("grid.refinement_level_macrogrid", 0);
    // for adaptive computations, the saved solution is not suitable for a later usage
    if (!DSC_CONFIG_GET("hmm.adaptivity", false) && DSC_CONFIG_GET("WRITE_HMM_SOL_TO_FILE", false)) {
      DiscreteFunctionWriter((boost::format("/hmm_solution_discFunc_refLevel_%d") % refinement_level_macrogrid_).str()
                             ).append(hmm_solution);
    }

    if (DSC_CONFIG_GET("fsr", true) && DSC_CONFIG_GET("WRITE_FINESCALE_SOL_TO_FILE", true))
    {
      const int refinement_level_referenceprob_ = typename HMM::ModelProblemDataType().getRefinementLevelReferenceProblem();
      DiscreteFunctionWriter((boost::format("/finescale_solution_discFunc_refLevel_%d") % refinement_level_referenceprob_).str()
          ).append(fem_newton_solution);
    }

    //! ******************** End of assembling and solving the HMM problem ***************************
    DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;
    //! ----------------- compute L2-errors -------------------
    if (DSC_CONFIG_GET("fsr", true))
    {
      DSC_PROFILER.startTiming("timeadapt");

      static const int hmm_polorder = 2* HMM::DiscreteFunctionSpaceType::polynomialOrder + 2;
      // expensive hack to deal with discrete functions, defined on different grids
      Dune::ImprovedL2Error< typename HMM::DiscreteFunctionType > impL2error;
      typename HMM::RangeType hmm_error = impL2error.template norm_adaptive_grids_2< hmm_polorder >(
        hmm_solution,
        fem_newton_solution);

      DSC_LOG_INFO << "|| u_hmm - u_fine_scale ||_L2 =  " << hmm_error << std::endl << std::endl;

      const auto timeadapt = DSC_PROFILER.stopTiming("timeadapt") / 1000.f;
      // if it took longer then 1 minute to compute the error:
      if (timeadapt > 60)
      {
        DSC_LOG_INFO << "WARNING! EXPENSIVE! Error assembled in " << timeadapt << "s." << std::endl << std::endl;
      }
    }

    #ifdef HMM_REFERENCE
    {
      DSC_PROFILER.startTiming("timeadapthmmref", hmm_iteration_step)

      RangeType hmm_vs_hmm_ref_error = impL2error.norm_adaptive_grids_2< hmm_polorder >(
        hmm_solution,
        hmm_reference_solution);

      DSC_LOG_INFO << "|| u_hmm - u_hmm_ref ||_L2 =  " << hmm_vs_hmm_ref_error << std::endl << std::endl;

      auto timeadapthmmref = DSC_PROFILER.stopTiming("timeadapthmmref") / 1000.f;
      // if it took longer then 1 minute to compute the error:
      if (timeadapthmmref > 60)
      {
        DSC_LOG_INFO << "WARNING! EXPENSIVE! Error assembled in " << timeadapthmmref << "s." << std::endl << std::endl;
      }
    }
    #endif // ifdef HMM_REFERENCE

    #ifdef HOMOGENIZEDSOL_AVAILABLE

    if (DSC_CONFIG_GET("fsr", true))
    {
      // not yet modified according to a generalized L2-error, here, homogenized_solution and fem_newton_solution still need
      // to be defined on the same grid!
      RangeType hom_error = l2error.norm2< hmm_polorder >(homogenized_solution, fem_newton_solution);
      DSC_LOG_INFO << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std::endl << std::endl;
      if ( DSC_LOG_INFO.is_open() )
      {
        DSC_LOG_INFO << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std::endl;
      }
    }

    RangeType hom_hmm_error = l2error.norm2< hmm_polorder >(hmm_solution,
                                                                                                 homogenized_solution);

    DSC_LOG_INFO << "|| u_hom - u_hmm ||_L2 =  " << hom_hmm_error << std::endl << std::endl;

    #endif // ifdef HOMOGENIZEDSOL_AVAILABLE

    if (HMM::ModelProblemDataType::has_exact_solution)
    {
      const typename HMM::ExactSolutionType u;
      const typename HMM::RangeType exact_hmm_error = l2error.template norm< typename HMM::ExactSolutionType >(u,
                                                                    hmm_solution,
                                                                    2 * HMM::DiscreteFunctionSpaceType::polynomialOrder + 2);

      DSC_LOG_INFO << "|| u_hmm - u_exact ||_L2 =  " << exact_hmm_error << std::endl << std::endl;
      if (DSC_CONFIG_GET("fsr", true))
      {
        typename HMM::RangeType fem_newton_error = l2error.template norm< typename HMM::ExactSolutionType >(u,
                                                                       fem_newton_solution,
                                                                       2 * HMM::DiscreteFunctionSpaceType::polynomialOrder + 2);

        DSC_LOG_INFO << "|| u_fem_newton - u_exact ||_L2 =  " << fem_newton_error << std::endl << std::endl;
      }
    }

    if (DSC_CONFIG_GET("ERRORESTIMATION", false)) {
      DSC_LOG_INFO << "Estimated error = " << errors.estimated_error << "." << std::endl;
      DSC_LOG_INFO << "In detail:" << std::endl;
      DSC_LOG_INFO << "   Estimated source error = " << errors.estimated_source_error << "." << std::endl;
      DSC_LOG_INFO << "   Estimated approximation error = " << errors.estimated_approximation_error << "." << std::endl;
      DSC_LOG_INFO << "   Estimated residual error = " << errors.estimated_residual_error << ", where:" << std::endl;
      DSC_LOG_INFO << "        contribution of macro jumps = " << errors.estimated_residual_error_macro_jumps << " and " << std::endl;
      DSC_LOG_INFO << "        contribution of micro jumps = " << errors.estimated_residual_error_micro_jumps << " and " << std::endl;
      if ( !DSC_CONFIG_GET("hmm.petrov_galerkin", true ) )
        DSC_LOG_INFO << "   Estimated tfr error = " << errors.estimated_tfr_error << "." << std::endl;
    }
    //! -------------------------------------------------------

    step_data_output<HMM>(gridPart, gridPartFine, hmm_solution, loop_cycle);
    return errors;
}

#endif // ALGORITHM_STEP_HH
