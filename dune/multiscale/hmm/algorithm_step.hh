#ifndef ALGORITHM_STEP_HH
#define ALGORITHM_STEP_HH

#include "algorithm_error.hh"


#include <string>
#include <fstream>

#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>


//! set the dirichlet points to zero
template< class EntityType, class DiscreteFunctionType >
void boundaryTreatment(const EntityType& entity, DiscreteFunctionType& rhs) {
  static const int faceCodim = 1;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType                          DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType                                  LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::LagrangePointSetType                          LagrangePointSetType;
  typedef typename DiscreteFunctionSpaceType::GridPartType                                  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType                                   IntersectionIteratorType;
  typedef typename LagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType FaceDofIteratorType;

  const DiscreteFunctionSpaceType& discreteFunctionSpace = rhs.space();
  const GridPartType& gridPart = discreteFunctionSpace.gridPart();
  IntersectionIteratorType it = gridPart.ibegin(entity);
  const IntersectionIteratorType endit = gridPart.iend(entity);
  for ( ; it != endit; ++it)
  {
    if ( !(*it).boundary() )
      continue;

    LocalFunctionType rhsLocal = rhs.localFunction(entity);
    const LagrangePointSetType& lagrangePointSet
      = discreteFunctionSpace.lagrangePointSet(entity);

    const int face = (*it).indexInInside();
    FaceDofIteratorType faceIterator
      = lagrangePointSet.template beginSubEntity< faceCodim >(face);
    const FaceDofIteratorType faceEndIterator
      = lagrangePointSet.template endSubEntity< faceCodim >(face);
    for ( ; faceIterator != faceEndIterator; ++faceIterator)
      rhsLocal[*faceIterator] = 0;
  }
} // boundaryTreatment

template < class HMMTraits >
HMMResult<HMMTraits>
    single_step(
        typename HMMTraits::GridPartType& gridPart,
        typename HMMTraits::GridPartType& gridPartFine,
        typename HMMTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
        typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
        const typename HMMTraits::DiffusionType& diffusion_op,
        const Dune::RightHandSideAssembler< typename HMMTraits::DiscreteFunctionType >& rhsassembler,
        std::ofstream& data_file,
        const std::string filename,
        typename HMMTraits::DiscreteFunctionType& hmm_solution
        )
{
    typedef HMMTraits HMM;
    DSC_LOG_INFO << std::endl << "Solving HMM-macro-problem for " << discreteFunctionSpace.size()
              << " unkowns and polynomial order "
              << HMM::DiscreteFunctionSpaceType::polynomialOrder << "."
              << std::endl << std::endl;

    // if we have some additional source term (-div G), define:
    const typename HMM::SecondSourceType G;
    // - div ( A^{\epsilon} \nabla u^{\epsilon} ) = f - div G
    // ! Ueberdenken, ob wir das nicht rausschmeisen und nur im Hintergrund fuer die Zellprobleme verwenden:
    // define mass (just for cell problems \lambda w - \div A \nabla w = rhs)
    const typename HMM::MassTermType mass;
    // dummy coefficient (mass, advection, etc.)
    const typename HMM::DefaultDummyFunctionType dummy_coeff;
    // define (first) source term:
    const typename HMM::FirstSourceType f;   // standard source f

    static const int hmm_polorder = 2* HMM::DiscreteFunctionSpaceType::polynomialOrder + 2;
    const Dune::L2Error< typename HMM::DiscreteFunctionType > l2error;
    // ----------------------------------------------------------------------------------------------//
    // ----------------------- THE DISCRETE HMM OPERATOR -----------------------------------//
    // ----------------------------------------------------------------------------------------------//

    // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain cell problem
    typename HMM::CellProblemNumberingManagerType cp_num_manager(discreteFunctionSpace);

    // ! define the elliptic hmm operator that describes our 'homogenized' macro problem
    // ( effect of the elliptic hmm operator on a certain discrete function )
    const typename HMM::EllipticHMMOperatorType discrete_elliptic_hmm_op(discreteFunctionSpace,
                                                     periodicDiscreteFunctionSpace,
                                                     diffusion_op,
                                                     cp_num_manager,
                                                     filename);

    // ----------------------------------------------------------------------------------------------//
    // ----------------------------------------------------------------------------------------------//
    // ----------------------------------------------------------------------------------------------//

    // ! matrix
    typename HMM::FEMMatrix hmm_newton_matrix("HMM Newton stiffness matrix", discreteFunctionSpace, discreteFunctionSpace);

    // ! right hand side vector
    // right hand side for the hm finite element method with Newton solver:
    typename HMM::DiscreteFunctionType hmm_newton_rhs("hmm rhs", discreteFunctionSpace);
    hmm_newton_rhs.clear();


    // starting value for the Newton method
    typename HMM::DiscreteFunctionType zero_func_coarse(filename + " constant zero function coarse ", discreteFunctionSpace);
    zero_func_coarse.clear();

    if (DSC_CONFIG.get("problem.linear", true)) {
      // solve cell problems in a preprocess, if AD_HOC_COMPUTATION is not defined
      #ifndef AD_HOC_COMPUTATION
      // ! -------------- solve and save the cell problems for the base function set --------------------------------------
      Dune::CellProblemSolver< HMM > cell_problem_solver(periodicDiscreteFunctionSpace,
                                                                                           diffusion_op,
                                                                                           data_file /*optinal*/);
      int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);
      DSC_LOG_INFO << "Solving cell problems for " << number_of_grid_elements << " leaf entities." << std::endl;
      // generate directory for cell problem data output
      std::string cell_path = "path";// + "/cell_problems/";assert(false);//path was not in scope
      Dune::Stuff::Common::Filesystem::testCreateDirectory(cell_path);
      // -------------- solve cell problems for the macro basefunction set ------------------------------
      // save the solutions of the cell problems for the set of macroscopic base functions
      cell_problem_solver.template saveTheSolutions_baseSet< typename HMM::DiscreteFunctionType >(discreteFunctionSpace,
                                                                           cp_num_manager,
                                                                           cell_path);
      // ------------- end solving and saving cell problems for the macro basefunction set --------------
      // ! --------------- end solving and saving cell problems -----------------------------------------
      #endif       // AD_HOC_COMPUTATION

      DSC_LOG_INFO << "Solving linear HMM problem." << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "Solving linear HMM problem." << std::endl;
        data_file << "------------------------------------------------------------------------------" << std::endl;
      }

      // to assemble the computational time
      Dune::Timer hmmAssembleTimer;

      // assemble the hmm stiffness matrix
      discrete_elliptic_hmm_op.assemble_matrix(hmm_newton_matrix);
      // to print the matrix, use:   hmm_newton_matrix.print();

      DSC_LOG_INFO << "Time to assemble HMM macro stiffness matrix: " << hmmAssembleTimer.elapsed() << "s" << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "Time to assemble HMM macro stiffness matrix: " << hmmAssembleTimer.elapsed() << "s" << std::endl;
      }

      // assemble right hand side
      rhsassembler.template assemble< 2* HMM::DiscreteFunctionSpaceType::polynomialOrder + 2 >(f, hmm_newton_rhs);

      // set Dirichlet Boundary to zero
      typedef typename HMM::DiscreteFunctionSpaceType::IteratorType IteratorType;
      IteratorType endit = discreteFunctionSpace.end();
      for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
        boundaryTreatment(*it, hmm_newton_rhs);

      typename HMM::InverseFEMMatrix hmm_biCGStab(hmm_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE);
      hmm_biCGStab(hmm_newton_rhs, hmm_solution);

      DSC_LOG_INFO << "Linear HMM problem solved in " << hmmAssembleTimer.elapsed() << "s." << std::endl << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "---------------------------------------------------------------------------------" << std::endl;
        data_file << "Linear HMM problem solved in " << hmmAssembleTimer.elapsed() << "s." << std::endl << std::endl
                  << std::endl;
      }

      } else {
      // the nonlinear case
      // solve cell problems in a preprocess, if AD_HOC_COMPUTATION is not defined
      #ifndef AD_HOC_COMPUTATION
      // ! -------------- solve and save the cell problems for the macroscopic base function set
      // --------------------------------------
      Dune::CellProblemSolver< HMMTraits > cell_problem_solver(periodicDiscreteFunctionSpace,
                                                                                           diffusion_op,
                                                                                           data_file /*optinal*/);

      int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);

      DSC_LOG_INFO << "Start solving cell problems for " << number_of_grid_elements << " leaf entities..." << std::endl;

      // generate directory for cell problem data output
      Dune::Stuff::Common::Filesystem::testCreateDirectory("data/HMM/" + filename + "/cell_problems/");

      // only for the case with test function reconstruction:
      #ifdef TFR
      // -------------- solve cell problems for the macro basefunction set ------------------------------
      // save the solutions of the cell problems for the set of macroscopic base functions

      cell_problem_solver.saveTheSolutions_baseSet< DiscreteFunctionType >(discreteFunctionSpace,
                                                                           cp_num_manager,
                                                                           filename + "/cell_problems/");

      DSC_LOG_INFO << "Solving the cell problems for the base function set succeeded." << std::endl;
      // end solving and saving cell problems
      #endif // ifdef TFR
             // ! --------------- end solving and saving cell problems -----------------------------------------
      #endif       // AD_HOC_COMPUTATION

      DSC_LOG_INFO << "Solving nonlinear HMM problem." << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "Solving nonlinear HMM problem with Newton method." << std::endl;
        data_file << "---------------------------------------------------------------------------------" << std::endl;
      }

      Dune::Timer hmmAssembleTimer;

      // just to provide some information
      typename HMM::PeriodicDiscreteFunctionType dummy_periodic_func("a periodic dummy", periodicDiscreteFunctionSpace);
      dummy_periodic_func.clear();

      // ! residual vector
      // current residual
      typename HMM::DiscreteFunctionType hmm_newton_residual(filename + "HMM Newton Residual", discreteFunctionSpace);
      hmm_newton_residual.clear();

      typename HMM::RangeType relative_newton_error = 10000.0;
      typename HMM::RangeType hmm_rhs_L2_norm = 10000.0;

      // number of HMM Newton step (1 = first step)
      // HMM_NEWTON_ITERATION_STEP' Netwon steps have been already performed,
      // the next one is 'HMM_NEWTON_ITERATION_STEP+1' = hmm_iteration_step
      int hmm_iteration_step = HMM_NEWTON_ITERATION_STEP + 1;

      #ifdef RESUME_TO_BROKEN_COMPUTATION
      // std :: string location_hmm_newton_step_solution = "data/HMM/test/hmm_solution_discFunc_refLevel_5_NewtonStep_2";
      char fnewtonname[50];
      sprintf(fnewtonname,
              "/hmm_solution_discFunc_refLevel_%d_NewtonStep_%d",
              refinement_level_macrogrid_,
              HMM_NEWTON_ITERATION_STEP);
      std::string fnewtonname_s(fnewtonname);
      std::string location_hmm_newton_step_solution = "data/HMM/" + filename + fnewtonname_s;

      bool reader_open = false;

      // reader for the cell problem data file:
      typename HMM::DiscreteFunctionReader discrete_function_reader_hmm_newton_ref( (location_hmm_newton_step_solution).c_str() );
      discrete_function_reader_hmm_newton_ref.open();

      discrete_function_reader_hmm_newton_ref.read(0, hmm_solution);
      #endif       // RESUME_TO_BROKEN_COMPUTATION

      double old_error = 100.0;
      double error_decay = 0.0;

      // the Newton step for the nonlinear HMM problem:
      // L2-Norm of residual < tolerance ?
      const double hmm_tolerance = DSC_CONFIG.get("problem.stochastic_pertubation", false)
                                    ? 1e-01 * DSC_CONFIG.get("problem.stochastic_variance",  0.01)
                                    : 1e-05;
      const int refinement_level_macrogrid_ = DSC_CONFIG.get("grid.refinement_level_macrogrid", 0);

      while (relative_newton_error > hmm_tolerance)
      {
        // (here: hmm_solution = solution from the last iteration step)
        long double newton_step_time = clock();
        DSC_LOG_INFO << "HMM Newton iteration " << hmm_iteration_step << ":" << std::endl;
        if ( data_file.is_open() )
        {
          data_file << "HMM Newton iteration " << hmm_iteration_step << ":" << std::endl;
        }

        #ifndef AD_HOC_COMPUTATION
        // solve cell problems for the solution of the last iteration step
        cell_problem_solver.template saveTheSolutions_discFunc< typename HMM::DiscreteFunctionType >(hmm_solution, filename + "/cell_problems/");
        cell_problem_solver.template saveTheJacCorSolutions_baseSet_discFunc< typename HMM::DiscreteFunctionType >(hmm_solution,
                                                                                            cp_num_manager,
                                                                                            filename + "/cell_problems/");
        #endif // ifndef AD_HOC_COMPUTATION

        // to assemble the computational time
        Dune::Timer stepHmmAssembleTimer;
        // assemble the stiffness matrix
        discrete_elliptic_hmm_op.assemble_jacobian_matrix(hmm_solution, hmm_newton_matrix);

        DSC_LOG_INFO << "Time to assemble HMM stiffness matrix for current Newton iteration: "
                  << stepHmmAssembleTimer.elapsed() << "s" << std::endl;
        if ( data_file.is_open() )
        {
          data_file << "Time to assemble HMM stiffness matrix for current Newton iteration: "
                    << stepHmmAssembleTimer.elapsed() << "s" << std::endl;
        }
        DSC_LOG_INFO << "Assemble right hand side..." << std::endl;
        // assemble right hand side
        const int assembler_order = 2* HMM::DiscreteFunctionSpaceType::polynomialOrder + 2;
        rhsassembler.template assemble_for_HMM_Newton_method< assembler_order >(
          f,
          diffusion_op,
          hmm_solution,
          cp_num_manager,
          dummy_periodic_func,
          hmm_newton_rhs,
          filename);
        DSC_LOG_INFO << "Right hand side assembled!" << std::endl;

        if ( !( hmm_newton_rhs.dofsValid() ) )
        {
          DSC_LOG_INFO << "Right hand side invalid!" << std::endl;
          data_file << "Right hand side invalid!" << std::endl;
          DUNE_THROW(Dune::InvalidStateException, "Right hand side invalid!");
        } else {
          DSC_LOG_INFO << "Right hand side valid ";
        }

        hmm_rhs_L2_norm = l2error.template norm2< assembler_order >(zero_func_coarse, hmm_newton_rhs);

        DSC_LOG_INFO << "with L^2-Norm = " << hmm_rhs_L2_norm << "." << std::endl;
        data_file << "Assembled right hand side, with L^2-Norm of RHS = " << hmm_rhs_L2_norm << "." << std::endl;

        if (hmm_rhs_L2_norm < 1e-10)
        {
          // residual solution almost identical to zero: break
          if ( data_file.is_open() )
          {
            data_file << "HMM residual solution almost identical to zero. Therefore: break loop." << std::endl;
            data_file << "(L^2-Norm of current right hand side = " << hmm_rhs_L2_norm << " < 1e-10)" << std::endl;
          }
          break;
        }

        // set Dirichlet Boundary to zero
        typedef typename HMM::DiscreteFunctionSpaceType::IteratorType IteratorType;
        IteratorType endit = discreteFunctionSpace.end();
        for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
          boundaryTreatment(*it, hmm_newton_rhs);

        #ifndef AD_HOC_COMPUTATION
        double hmm_biCG_tolerance = 1e-8;
        bool hmm_solution_convenient = false;
        while (hmm_solution_convenient == false)
        {
          hmm_newton_residual.clear();
          typename HMM::InverseFEMMatrix hmm_newton_biCGStab(hmm_newton_matrix,
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
        #else         // AD_HOC_COMPUTATION
        InverseFEMMatrix hmm_newton_biCGStab(hmm_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE);
        hmm_newton_biCGStab(hmm_newton_rhs, hmm_newton_residual);
        #endif         // AD_HOC_COMPUTATION


        if ( hmm_newton_residual.dofsValid() )
        {
          hmm_solution += hmm_newton_residual;

          // write the solution after the current HMM Newton step to a file
          #ifdef WRITE_HMM_SOL_TO_FILE
          // for adaptive computations, the saved solution is not suitable for a later usage
          #ifndef ADAPTIVE
          char fname[50];
          sprintf(fname,
                  "/hmm_solution_discFunc_refLevel_%d_NewtonStep_%d",
                  refinement_level_macrogrid_,
                  hmm_iteration_step);
          std::string fname_s(fname);

          std::string location = "data/HMM/" + filename + fname_s;
          DiscreteFunctionWriter dfw( (location).c_str() );
          if (dfw.is_open())
            dfw.append(hmm_solution);

          // if you want an utput for all newton steps, even for an adaptive computation, use:
          // #endif

          // writing paraview data output

          // general output parameters
          myDataOutputParameters outputparam;
          outputparam.set_path("data/HMM/" + filename);

          // sequence stamp
          std::stringstream outstring;

          // create and initialize output class
          typename HMM::IOTupleType hmm_solution_newton_step_series(&hmm_solution);
          #ifdef ADAPTIVE
          char hmm_prefix[50];
          sprintf(hmm_prefix, "hmm_solution_%d_NewtonStep_%d", loop_cycle, hmm_iteration_step);
          #else // ifdef ADAPTIVE
          char hmm_prefix[50];
          sprintf(hmm_prefix, "hmm_solution_NewtonStep_%d", hmm_iteration_step);
          #endif // ifdef ADAPTIVE
          outputparam.set_prefix(hmm_prefix);
          typename HMM::DataOutputType hmmsol_dataoutput(gridPart.grid(), hmm_solution_newton_step_series, outputparam);

          // write data
          outstring << "hmm-solution-NewtonStep";
          hmmsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
          // clear the std::stringstream:
          outstring.str( std::string() );
          #endif               // ADAPTIVE
          #endif           // WRITE_HMM_SOL_TO_FILE

          // || u^(n+1) - u^(n) ||_L2
          relative_newton_error = l2error.template norm2< assembler_order >(hmm_newton_residual,
                                                                                                     zero_func_coarse);
          // || u^(n+1) - u^(n) ||_L2 / || u^(n+1) ||_L2
          relative_newton_error = relative_newton_error / l2error.template norm2< assembler_order >(
            hmm_solution,
            zero_func_coarse);

          DSC_LOG_INFO << "Relative L2 HMM Newton iteration error = " << relative_newton_error << std::endl;

          // residual solution almost identical to zero: break
          if ( data_file.is_open() )
          {
            data_file << "Relative L2 HMM Newton iteration error = " << relative_newton_error << std::endl;
            if (relative_newton_error <= hmm_tolerance)
            {
              newton_step_time = clock() - newton_step_time;
              newton_step_time = newton_step_time / CLOCKS_PER_SEC;
              if ( data_file.is_open() )
              {
                data_file << std::endl << "Total time for current HMM Newton step = " << newton_step_time << "s."
                          << std::endl << std::endl;
              }
              data_file << "Since HMM-tolerance = " << hmm_tolerance << ": break loop." << std::endl;
              data_file << "....................................................." << std::endl << std::endl;
            }
          }

          hmm_newton_residual.clear();
        } else {
          DSC_LOG_INFO << "WARNING! Invalid dofs in 'hmm_newton_residual'." << std::endl;
          break;
        }

        hmm_iteration_step += 1;

        if (relative_newton_error > hmm_tolerance)
        {
          newton_step_time = clock() - newton_step_time;
          newton_step_time = newton_step_time / CLOCKS_PER_SEC;
          if ( data_file.is_open() )
          {
            data_file << std::endl << "Total time for current HMM Newton step = " << newton_step_time << "s."
                      << std::endl << std::endl;

            error_decay = relative_newton_error / old_error;
            old_error = relative_newton_error;
            // maximum number of Newton iterations
            if ( (hmm_iteration_step >= 20) || (error_decay >= 0.95) )
            {
              data_file << std::endl
                        << "Reached constant or inceasing error decay or maximum number of Newton iterations:  break loop."
                        << std::endl;
              data_file << "....................................................." << std::endl << std::endl;
              break;
            }
            data_file << "....................................................." << std::endl << std::endl;
          }
        }
      }       // while( relative_newton_error > hmm_tolerance )

      DSC_LOG_INFO << "HMM problem with Newton method solved in " << hmmAssembleTimer.elapsed() << "s." << std::endl
                << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "---------------------------------------------------------------------------------" << std::endl;
        data_file << "HMM problem with Newton method solved in " << hmmAssembleTimer.elapsed() << "s." << std::endl
                  << std::endl << std::endl;
      }
      #ifdef ADAPTIVE
      total_hmm_time += hmmAssembleTimer.elapsed();
      #endif
    } //if not problem.linear

    const auto retval = estimate_error<HMM>(gridPart, gridPartFine, discreteFunctionSpace, periodicDiscreteFunctionSpace,
                   diffusion_op, rhsassembler, data_file, filename, cp_num_manager, hmm_solution);

    const int refinement_level_macrogrid_ = DSC_CONFIG.get("grid.refinement_level_macrogrid", 0);
    #ifdef WRITE_HMM_SOL_TO_FILE
    // for adaptive computations, the saved solution is not suitable for a later usage
    #ifndef ADAPTIVE
    char fname[40];
    sprintf(fname, "/hmm_solution_discFunc_refLevel_%d", refinement_level_macrogrid_);
    std::string fname_s(fname);
    std::string location = "data/HMM/" + filename + fname_s;
    DiscreteFunctionWriter dfw( (location).c_str() );
    if (dfw.is_open())
      dfw.append(hmm_solution);
    #endif // ifndef ADAPTIVE
    #endif   // #ifdef WRITE_HMM_SOL_TO_FILE

    #ifdef FINE_SCALE_REFERENCE
    #ifdef WRITE_FINESCALE_SOL_TO_FILE
    bool fine_writer_is_open = false;
    char fine_fname[50];
    sprintf(fine_fname, "/finescale_solution_discFunc_refLevel_%d", refinement_level_referenceprob_);
    std::string fine_fname_s(fine_fname);
    std::string fine_location = "data/HMM/" + filename + fine_fname_s;
    DiscreteFunctionWriter fine_dfw( (fine_location).c_str() );
    fine_writer_is_open = fine_dfw.open();
    if (fine_writer_is_open)
      fine_dfw.append(fem_newton_solution);
    #endif // ifdef WRITE_FINESCALE_SOL_TO_FILE
    #endif // ifdef FINE_SCALE_REFERENCE

    // ! ******************** End of assembling and solving the HMM problem ***************************

    DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "The L2 errors:" << std::endl << std::endl;
    }

    // ! ----------------- compute L2-errors -------------------

    #ifdef FINE_SCALE_REFERENCE
    long double timeadapt = clock();

    RangeType hmm_error = impL2error.norm_adaptive_grids_2< hmm_polorder >(
      hmm_solution,
      fem_newton_solution);

    DSC_LOG_INFO << "|| u_hmm - u_fine_scale ||_L2 =  " << hmm_error << std::endl << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "|| u_hmm - u_fine_scale ||_L2 =  " << hmm_error << std::endl;
    }

    timeadapt = clock() - timeadapt;
    timeadapt = timeadapt / CLOCKS_PER_SEC;

    // if it took longer then 1 minute to compute the error:
    if (timeadapt > 60)
    {
      DSC_LOG_INFO << "WARNING! EXPENSIVE! Error assembled in " << timeadapt << "s." << std::endl << std::endl;

      if ( data_file.is_open() )
      {
        DSC_LOG_INFO << "WARNING! EXPENSIVE! Error assembled in " << timeadapt << "s." << std::endl << std::endl;
      }
    }
    #endif // ifdef FINE_SCALE_REFERENCE

    #ifdef HMM_REFERENCE
    long double timeadapthmmref = clock();

    RangeType hmm_vs_hmm_ref_error = impL2error.norm_adaptive_grids_2< hmm_polorder >(
      hmm_solution,
      hmm_reference_solution);

    DSC_LOG_INFO << "|| u_hmm - u_hmm_ref ||_L2 =  " << hmm_vs_hmm_ref_error << std::endl << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "|| u_hmm - u_hmm_ref ||_L2 =  " << hmm_vs_hmm_ref_error << std::endl;
    }

    timeadapthmmref = clock() - timeadapthmmref;
    timeadapthmmref = timeadapthmmref / CLOCKS_PER_SEC;

    // if it took longer then 1 minute to compute the error:
    if (timeadapthmmref > 60)
    {
      DSC_LOG_INFO << "WARNING! EXPENSIVE! Error assembled in " << timeadapthmmref << "s." << std::endl << std::endl;

      if ( data_file.is_open() )
      {
        DSC_LOG_INFO << "WARNING! EXPENSIVE! Error assembled in " << timeadapthmmref << "s." << std::endl << std::endl;
      }
    }
    #endif // ifdef HMM_REFERENCE

    #ifdef HOMOGENIZEDSOL_AVAILABLE

    #ifdef FINE_SCALE_REFERENCE
    // not yet modified according to a generalized L2-error, here, homogenized_solution and fem_newton_solution still need
    // to be defined on the same grid!
    RangeType hom_error = l2error.norm2< hmm_polorder >(homogenized_solution,
                                                                                             fem_newton_solution);

    DSC_LOG_INFO << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std::endl << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std::endl;
    }
    #endif // ifdef FINE_SCALE_REFERENCE

    RangeType hom_hmm_error = l2error.norm2< hmm_polorder >(hmm_solution,
                                                                                                 homogenized_solution);

    DSC_LOG_INFO << "|| u_hom - u_hmm ||_L2 =  " << hom_hmm_error << std::endl << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "|| u_hom - u_hmm ||_L2 =  " << hom_hmm_error << std::endl;
    }
    #endif // ifdef HOMOGENIZEDSOL_AVAILABLE

    if (HMM::ModelProblemDataType::has_exact_solution)
    {
      const typename HMM::ExactSolutionType u;
      const typename HMM::RangeType exact_hmm_error = l2error.template norm< typename HMM::ExactSolutionType >(u,
                                                                    hmm_solution,
                                                                    2 * HMM::DiscreteFunctionSpaceType::polynomialOrder + 2);

      DSC_LOG_INFO << "|| u_hmm - u_exact ||_L2 =  " << exact_hmm_error << std::endl << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "|| u_hmm - u_exact ||_L2 =  " << exact_hmm_error << std::endl;
      }

      #ifdef FINE_SCALE_REFERENCE
      typename HMM::RangeType fem_newton_error = l2error.norm< ExactSolutionType >(u,
                                                                     fem_newton_solution,
                                                                     2 * HMM::DiscreteFunctionSpaceType::polynomialOrder + 2);

      DSC_LOG_INFO << "|| u_fem_newton - u_exact ||_L2 =  " << fem_newton_error << std::endl << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "|| u_fem_newton - u_exact ||_L2 =  " << fem_newton_error << std::endl;
      }
      #endif // ifdef FINE_SCALE_REFERENCE
    }


    #ifdef ERRORESTIMATION
    DSC_LOG_INFO << "Estimated error = " << estimated_error << "." << std::endl;
    DSC_LOG_INFO << "In detail:" << std::endl;
    DSC_LOG_INFO << "   Estimated source error = " << estimated_source_error << "." << std::endl;
    DSC_LOG_INFO << "   Estimated approximation error = " << estimated_approximation_error << "." << std::endl;
    DSC_LOG_INFO << "   Estimated residual error = " << estimated_residual_error << ", where:" << std::endl;
    DSC_LOG_INFO << "        contribution of macro jumps = " << estimated_residual_error_macro_jumps << " and " << std::endl;
    DSC_LOG_INFO << "        contribution of micro jumps = " << estimated_residual_error_micro_jumps << " and " << std::endl;
    #ifdef TFR
    DSC_LOG_INFO << "   Estimated tfr error = " << estimated_tfr_error << "." << std::endl;
    #endif
    if ( data_file.is_open() )
    {
      data_file << "Estimated error = " << estimated_error << "." << std::endl;
      data_file << "In detail:" << std::endl;
      data_file << "   Estimated source error = " << estimated_source_error << "." << std::endl;
      data_file << "   Estimated approximation error = " << estimated_approximation_error << "." << std::endl;
      data_file << "   Estimated residual error = " << estimated_residual_error << ", where:" << std::endl;
      data_file << "        contribution of macro jumps = " << estimated_residual_error_macro_jumps << " and "
                << std::endl;
      data_file << "        contribution of micro jumps = " << estimated_residual_error_micro_jumps << " and "
                << std::endl;
      #ifdef TFR
      data_file << "   Estimated tfr error = " << estimated_tfr_error << "." << std::endl;
      #endif
    }
    #endif // ifdef ERRORESTIMATION
    // ! -------------------------------------------------------

    // ! --------------- writing data output ---------------------
    // general output parameters
    myDataOutputParameters outputparam;
    outputparam.set_path("data/HMM/" + filename);

    // sequence stamp
    std::stringstream outstring;

    // --------- data output hmm solution --------------

    // create and initialize output class
    typename HMM::IOTupleType hmm_solution_series(&hmm_solution);
    #ifdef ADAPTIVE
    char hmm_prefix[30];
    sprintf(hmm_prefix, "hmm_solution_%d", loop_cycle);
    outputparam.set_prefix(hmm_prefix);
    #else // ifdef ADAPTIVE
    outputparam.set_prefix("hmm_solution");
    #endif // ifdef ADAPTIVE
    typename HMM::DataOutputType hmmsol_dataoutput(gridPart.grid(), hmm_solution_series, outputparam);

    // write data
    outstring << "hmm-solution";
    hmmsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
    // clear the std::stringstream:
    outstring.str( std::string() );
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
    // !-------------------------------------------------------------
    return retval;
}

#endif // ALGORITHM_STEP_HH
