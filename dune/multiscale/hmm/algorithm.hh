#ifndef DUNE_MS_HMM_ALGORITHM_HH
#define DUNE_MS_HMM_ALGORITHM_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/solver.hh>
#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>
#include <dune/stuff/common/logging.hh>

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

//! \todo DOCME
template< class Stream, class DiscFunc >
void oneLinePrint(Stream& stream, const DiscFunc& func) {
  typedef typename DiscFunc::ConstDofIteratorType
  DofIteratorType;
  DofIteratorType it = func.dbegin();
  stream << "\n" << func.name() << ": [ ";
  for ( ; it != func.dend(); ++it)
    stream << std::setw(5) << *it << "  ";

  stream << " ] " << std::endl;
} // oneLinePrint

//! \todo replace me with Stuff::something
template < class DiscreteFunctionSpaceType >
typename DiscreteFunctionSpaceType::RangeType get_size_of_domain(DiscreteFunctionSpaceType& discreteFunctionSpace) {
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  IteratorType endit = discreteFunctionSpace.end();
  typename DiscreteFunctionSpaceType::RangeType size_of_domain = 0.0;
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    Dune::CachingQuadrature< GridPartType, 0 > entityQuadrature(*it, 0);
    // get geoemetry of entity
    const typename GridPartType::GridType::template Codim< 0 >::Geometry& geometry = it->geometry();
    const double volumeEntity = entityQuadrature.weight(0)
                                * geometry.integrationElement( entityQuadrature.point(0) );
    size_of_domain += volumeEntity;
  }
  return size_of_domain;
} // get_size_of_domain

template < class ProblemDataType >
void print_info(const ProblemDataType& info, std::ofstream& data_file)
{
  // epsilon is specified in ModelProblemData, which is specified in problem_specification.hh
  // 'epsilon' in for instance A^{epsilon}(t,x) = A(t,x/epsilon)
  const double epsilon_ = info.getEpsilon();
  // estimated epsilon (specified in ModelProblemData)
  // estimated epsilon in case epsilon is unknown
  const double epsilon_est_ = info.getEpsilonEstimated();
  // edge length of the cells in the cells, belonging to the cell problems
  // note that (delta/epsilon_est) needs to be a positive integer!
  // edge length of the cells in the cell proplems,
  const double delta_ = info.getDelta();
  const int refinement_level_macrogrid_ = DSC_CONFIG.get("grid.refinement_level_macrogrid", 0);
  if ( data_file.is_open() )
  {
    data_file << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << "." << std::endl
              << std::endl;
    if (DSC_CONFIG.get("problem.linear", true))
      data_file << "Problem is declared as being LINEAR." << std::endl;
    else
      data_file << "Problem is declared as being NONLINEAR." << std::endl;

    if (ProblemDataType::has_exact_solution) {
      data_file << "Exact solution is available." << std::endl << std::endl;
    } else {
      data_file << "Exact solution is not available." << std::endl << std::endl;
    }
    data_file << "Computations were made for:" << std::endl << std::endl;
    data_file << "Refinement Level for (uniform) Macro Grid = " << refinement_level_macrogrid_ << std::endl;
    const int refinement_level_cellgrid = DSC_CONFIG.get("grid.refinement_level_cellgrid", 1);
    data_file << "Refinement Level for Periodic Micro Grid = " << refinement_level_cellgrid << std::endl << std::endl;
    #ifdef TFR
    data_file << "We use TFR-HMM (HMM with test function reconstruction)." << std::endl;
    #else
    data_file << "We use HMM without test function reconstruction (NO TFR)." << std::endl;
    #endif // ifdef TFR
    #ifdef AD_HOC_COMPUTATION
    data_file << "Cell problems are solved ad hoc (where required)." << std::endl << std::endl;
    #else
    data_file << "Cell problems are solved and saved (in a pre-process)." << std::endl << std::endl;
    #ifdef ERRORESTIMATION
    data_file << "Error estimation activated!" << std::endl << std::endl;
    #endif
    #endif // ifdef AD_HOC_COMPUTATION
    data_file << "Epsilon = " << epsilon_ << std::endl;
    data_file << "Estimated Epsilon = " << epsilon_est_ << std::endl;
    data_file << "Delta (edge length of cell-cube) = " << delta_ << std::endl;
    if (DSC_CONFIG.get("problem.stochastic_pertubation", false))
      data_file << std::endl << "Stochastic perturbation added. Variance = " << DSC_CONFIG.get("problem.stochastic_variance", 0.01) << std::endl;
    #ifdef ADAPTIVE
    //only used in adaptive config
    const double error_tolerance_ = DSC_CONFIG.get("problem.error_tolerance", 1e-6);
    data_file << std::endl << "Adaptive computation. Global error tolerance for program abort = "
              << error_tolerance_ << std::endl;
    #endif // ifdef ADAPTIVE
    data_file << std::endl << std::endl;
  }
}

//! the main hmm computation
template < class ProblemDataType, class HMMTraits >
void algorithm(const ProblemDataType& problem_data,
               const std::string& /*UnitCubeName*/,
               typename HMMTraits::GridPointerType& macro_grid_pointer,   // grid pointer that belongs to the macro grid
               typename HMMTraits::GridPointerType& fine_macro_grid_pointer,   // grid pointer that belongs to the fine macro grid (for
                                                           // reference computations)
               typename HMMTraits::GridPointerType& periodic_grid_pointer,   // grid pointer that belongs to the periodic micro grid
               int /*refinement_difference*/,   // refinement difference for the macro grid (problem-to-solve vs. reference
                                            // problem)
               std::ofstream& data_file,
               const std::string filename) {
  typedef HMMTraits HMM;
  print_info(problem_data, data_file);
  // ! ---- tools ----
  // model problem data
// UNUSED  Problem::ModelProblemData problem_info;
// set of hmm parameters/information
// UNUSED  Multiscale::HMMParameters method_info;

  Dune::L2Error< typename HMM::DiscreteFunctionType > l2error;
  // expensive hack to deal with discrete functions, defined on different grids
// UNUSED  ImprovedL2Error< DiscreteFunctionType > impL2error;

  // ! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for HMM-macro-problem
  typename HMM::GridPartType gridPart(*macro_grid_pointer);
  // grid part for the periodic function space, required for HMM-cell-problems
  typename HMM::PeriodicGridPartType periodicGridPart(*periodic_grid_pointer);
  // auxiliary grid part for the periodic function space, required for HMM-cell-problems
  typename HMM::GridPartType auxiliaryGridPart(*periodic_grid_pointer);
  // auxiliaryGridPart for the error estimator (the auxiliaryGridPart yields an intersection iterator, which can not be
  // done by the periodicGridPart)
  // grid part for the global function space, required for the detailed fine-scale computation (very high resolution)
  typename HMM::GridPartType gridPartFine(*fine_macro_grid_pointer);

// UNUSED  GridType& grid = gridPart.grid();
// UNUSED  GridType& gridFine = gridPartFine.grid();
// ! --------------------------------------------------------------------------------------

  // ! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  typename HMM::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  // the global-problem function space for the reference computation:
  typename HMM::DiscreteFunctionSpaceType finerDiscreteFunctionSpace(gridPartFine);
  // the local-problem function space (containing periodic functions):
  typename HMM::PeriodicDiscreteFunctionSpaceType periodicDiscreteFunctionSpace(periodicGridPart);
  // and the corresponding auxiliary one:
// UNUSED  DiscreteFunctionSpaceType auxiliaryDiscreteFunctionSpace(auxiliaryGridPart);
// ! --------------------------------------------------------------------------------------

  // ! --------------------------- coefficient functions ------------------------------------
  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  typename HMM::DiffusionType diffusion_op;
  // define (first) source term:
  typename HMM::FirstSourceType f;   // standard source f
  // if we have some additional source term (-div G), define:
  typename HMM::SecondSourceType G;
  // - div ( A^{\epsilon} \nabla u^{\epsilon} ) = f - div G
  // ! Ueberdenken, ob wir das nicht rausschmeisen und nur im Hintergrund fuer die Zellprobleme verwenden:
  // define mass (just for cell problems \lambda w - \div A \nabla w = rhs)
  typename HMM::MassTermType mass;
  // dummy coefficient (mass, advection, etc.)
  typename HMM::DefaultDummyFunctionType dummy_coeff;
  // exact solution unknown?

  // ! --------------------------------------------------------------------------------------

  // ! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  Dune::RightHandSideAssembler< typename HMM::DiscreteFunctionType > rhsassembler;

  // ----------------------------------------------------------------------------------------------//
  // ----------------------- THE DISCRETE FEM OPERATOR -----------------------------------//
  // ----------------------------------------------------------------------------------------------//
  // ! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
// UNUSED EllipticOperatorType discrete_elliptic_op(finerDiscreteFunctionSpace, diffusion_op);
// ----------------------------------------------------------------------------------------------//
// ----------------------------------------------------------------------------------------------//
// ----------------------------------------------------------------------------------------------//
// UNUSED  RangeType size_of_domain = get_size_of_domain(discreteFunctionSpace);

  #ifdef HOMOGENIZEDSOL_AVAILABLE
  if (DSC_CONFIG.get("problem.linear", true)) {
    std::string unit_cell_location = "../dune/multiscale/grids/cell_grids/unit_cube.dgf";
    Dune::FieldMatrix< RangeType, dimension, dimension > A_hom;

    // descretized homogenizer:
    Homogenizer< GridType, DiffusionType > disc_homogenizer(unit_cell_location);
    A_hom = disc_homogenizer.getHomTensor(diffusion_op);

    typedef FieldMatrix< RangeType, dimension, dimension > FieldMatrixType;
    Problem::HomDiffusion< FunctionSpaceType, FieldMatrixType > hom_diffusion_op(A_hom);

    typedef DiscreteEllipticOperator< DiscreteFunctionType,
                                      Problem::HomDiffusion< FunctionSpaceType,
                                                             FieldMatrixType >, MassTermType > HomEllipticOperatorType;

    HomEllipticOperatorType hom_discrete_elliptic_op(finerDiscreteFunctionSpace, hom_diffusion_op);

    FEMMatrix hom_stiff_matrix("homogenized stiffness matrix", finerDiscreteFunctionSpace, finerDiscreteFunctionSpace);

    DiscreteFunctionType hom_rhs("homogenized rhs", finerDiscreteFunctionSpace);
    hom_rhs.clear();

    DiscreteFunctionType homogenized_solution(filename + " Homogenized Solution", finerDiscreteFunctionSpace);
    homogenized_solution.clear();
    hom_discrete_elliptic_op.assemble_matrix(hom_stiff_matrix);

    rhsassembler.assemble< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(f, hom_rhs);

    // set Dirichlet Boundary to zero
    typedef DiscreteFunctionSpaceType::IteratorType IteratorType;
    IteratorType hom_endit = finerDiscreteFunctionSpace.end();
    for (IteratorType fine_it = finerDiscreteFunctionSpace.begin(); fine_it != hom_endit; ++fine_it)
      boundaryTreatment(*fine_it, hom_rhs);

    InverseFEMMatrix hom_biCGStab(hom_stiff_matrix, 1e-8, 1e-8, 20000, VERBOSE);
    hom_biCGStab(hom_rhs, homogenized_solution);
  } else {

  }
  #endif   // HOMOGENIZEDSOL_AVAILABLE

  #ifdef FINE_SCALE_REFERENCE
  #ifdef FSR_COMPUTE
  // ! *******************************************************************

  // starting value for the Newton method
  DiscreteFunctionType zero_func(filename + " constant zero function ", finerDiscreteFunctionSpace);
  zero_func.clear();

  // ! *************************** Assembling the reference problem ****************************
  // ( fine scale reference solution = fem_newton_solution )

  // ! (stiffness) matrix
  FEMMatrix fem_newton_matrix("FEM Newton stiffness matrix", finerDiscreteFunctionSpace, finerDiscreteFunctionSpace);

  // ! right hand side vector
  // right hand side for the finite element method with Newton solver:
  // ( also right hand side for the finer discrete function space )
  DiscreteFunctionType fem_newton_rhs("fem newton rhs", finerDiscreteFunctionSpace);
  fem_newton_rhs.clear();

  // ! solution vector
  // solution of the finite element method, where we used the Newton method to solve the non-linear system of equations
  // in general this will be an accurate approximation of the exact solution, that is why we it also called reference
  // solution
  DiscreteFunctionType fem_newton_solution(filename + " Reference (FEM Newton) Solution", finerDiscreteFunctionSpace);
  fem_newton_solution.clear();
  // By fem_newton_solution, we denote the "fine scale reference solution" (used for comparison)
  // ( if the elliptic problem is linear, the 'fem_newton_solution' is determined without the Newton method )

  if (DSC_CONFIG.get("problem.linear", true)) {
    DSC_LOG_INFO << "Solving linear problem." << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "Solving linear problem with standard FEM and resolution level "
                << problem_info.getRefinementLevelReferenceProblem() << "." << std::endl;
      data_file << "------------------------------------------------------------------------------" << std::endl;
    }

    // to assemble the computational time
    Dune::Timer assembleTimer;

    // assemble the stiffness matrix
    discrete_elliptic_op.assemble_matrix(fem_newton_matrix);

    DSC_LOG_INFO << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
    }

    // assemble right hand side
    rhsassembler.assemble< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(f, fem_newton_rhs);

    // set Dirichlet Boundary to zero
    typedef DiscreteFunctionSpaceType::IteratorType IteratorType;
    IteratorType fine_endit = finerDiscreteFunctionSpace.end();
    for (IteratorType fine_it = finerDiscreteFunctionSpace.begin(); fine_it != fine_endit; ++fine_it)
      boundaryTreatment(*fine_it, fem_newton_rhs);

    InverseFEMMatrix fem_biCGStab(fem_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE);
    fem_biCGStab(fem_newton_rhs, fem_newton_solution);

    if ( data_file.is_open() )
    {
      data_file << "---------------------------------------------------------------------------------" << std::endl;
      data_file << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
                << std::endl;
    }
  } else {
    DSC_LOG_INFO << "Solving non-linear problem." << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "Solving nonlinear problem with FEM + Newton-Method. Resolution level of grid = "
                << problem_info.getRefinementLevelReferenceProblem() << "." << std::endl;
      data_file << "---------------------------------------------------------------------------------" << std::endl;
    }

    Dune::Timer assembleTimer;

    // ! residual vector
    // current residual
    DiscreteFunctionType fem_newton_residual(filename + "FEM Newton Residual", finerDiscreteFunctionSpace);
    fem_newton_residual.clear();

    RangeType relative_newton_error_finescale = 10000.0;
    RangeType rhs_L2_norm = 10000.0;

    int iteration_step = 1;
    // the Newton step for the FEM reference problem (solved with Newton Method):
    // L2-Norm of residual < tolerance ?
    double tolerance = 1e-06;
    while (relative_newton_error_finescale > tolerance)
    {
      // (here: fem_newton_solution = solution from the last iteration step)
      DSC_LOG_INFO << "Newton iteration " << iteration_step << ":" << std::endl;
      if ( data_file.is_open() )
      {
        data_file << "Newton iteration " << iteration_step << ":" << std::endl;
      }
      Dune::Timer stepAssembleTimer;
      // assemble the stiffness matrix
      discrete_elliptic_op.assemble_jacobian_matrix(fem_newton_solution, fem_newton_matrix);

      DSC_LOG_INFO << "Time to assemble FEM Newton stiffness matrix for current iteration: "
                << stepAssembleTimer.elapsed() << "s" << std::endl;

      // assemble right hand side
      rhsassembler.assemble_for_Newton_method< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(f,
                                                                                                   diffusion_op,
                                                                                                   fem_newton_solution,
                                                                                                   fem_newton_rhs);

      rhs_L2_norm = l2error.norm2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(fem_newton_rhs, zero_func);
      if (rhs_L2_norm < 1e-10)
      {
        // residual solution almost identical to zero: break
        if ( data_file.is_open() )
        {
          data_file << "Residual solution almost identical to zero. Therefore: break loop." << std::endl;
          data_file << "(L^2-Norm of current right hand side = " << rhs_L2_norm << " < 1e-10)" << std::endl;
        }
        break;
      }
      // set Dirichlet Boundary to zero
      typedef DiscreteFunctionSpaceType::IteratorType IteratorType;
      IteratorType fine_endit = finerDiscreteFunctionSpace.end();
      for (IteratorType fine_it = finerDiscreteFunctionSpace.begin(); fine_it != fine_endit; ++fine_it)
        boundaryTreatment(*fine_it, fem_newton_rhs);

      InverseFEMMatrix fem_newton_biCGStab(fem_newton_matrix, 1e-8, 1e-8, 20000, true);
      fem_newton_biCGStab(fem_newton_rhs, fem_newton_residual);

      if ( fem_newton_residual.dofsValid() )
      {
        fem_newton_solution += fem_newton_residual;
        relative_newton_error_finescale = l2error.norm2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(
          fem_newton_residual,
          zero_func);
        relative_newton_error_finescale /= l2error.norm2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(
          fem_newton_solution,
          zero_func);

        DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
        // residual solution almost identical to zero: break
        if ( data_file.is_open() )
        {
          data_file << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
          if (relative_newton_error_finescale <= tolerance)
          {
            data_file << "Since tolerance = " << tolerance << ": break loop." << std::endl;
          }
        }
        fem_newton_residual.clear();
      } else {
        DSC_LOG_INFO << "WARNING! Invalid dofs in 'fem_newton_residual'." << std::endl;
        break;
      }
      iteration_step += 1;
    }
    DSC_LOG_INFO << "Problem with FEM + Newton-Method solved in " << assembleTimer.elapsed() << "s." << std::endl
              << std::endl;
    if ( data_file.is_open() )
    {
      data_file << "---------------------------------------------------------------------------------" << std::endl;
      data_file << "Problem with FEM + Newton-Method solved in " << assembleTimer.elapsed() << "s." << std::endl
                << std::endl << std::endl;
    }
  }// end 'problem.linear <-> else'

  // ! ********************** End of assembling the reference problem ***************************
  #endif       // FSR_COMPUTE

  #ifdef FSR_LOAD

  DiscreteFunctionType fem_newton_solution(filename + " Reference (FEM Newton) Solution", finerDiscreteFunctionSpace);
  fem_newton_solution.clear();

  char modeprob[50];
  sprintf( modeprob, "/Model_Problem_%d", problem_info.get_Number_of_Model_Problem() );
  std::string modeprob_s(modeprob);

  char reference_solution_directory[50];
  sprintf(reference_solution_directory, "/reference_solution_ref_%d", refinement_level_referenceprob_);
  std::string reference_solution_directory_s(reference_solution_directory);

  char reference_solution_name[50];
  sprintf(reference_solution_name, "/finescale_solution_discFunc_refLevel_%d", refinement_level_referenceprob_);
  std::string reference_solution_name_s(reference_solution_name);

  std::string location_fine_scale_ref = "data/HMM/" + modeprob_s + reference_solution_directory_s
                                        + reference_solution_name_s;

  bool reader_is_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_ref( (location_fine_scale_ref).c_str() );
  discrete_function_reader_ref.open();

  discrete_function_reader_ref.read(0, fem_newton_solution);
  DSC_LOG_INFO << "fine scale reference read." << std::endl;

  #endif       // FSR_LOADd
  #endif   // FINE_SCALE_REFERENCE

  // noch per Hand die Daten eingetragen:
  #ifdef HMM_REFERENCE
  const int gridLevel_refHMM = 10;       // Macro_'gridLevel'

  std::string macroGridName_refHMM;
  problem_info.getMacroGridFile(macroGridName_refHMM);
  DSC_LOG_INFO << "loading dgf: " << macroGridName_refHMM << std::endl;

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer_refHMM(macroGridName_refHMM);
  // refine the grid 'gridLevel_refHMM' times:
  macro_grid_pointer_refHMM->globalRefine(gridLevel_refHMM);

  GridPartType gridPart_refHMM(*macro_grid_pointer_refHMM);
  GridType& grid_refHMM = gridPart_refHMM.grid();

  DiscreteFunctionSpaceType discreteFunctionSpace_refHMM(gridPart_refHMM);

  DiscreteFunctionType hmm_reference_solution(filename + " Reference (HMM) Solution", discreteFunctionSpace_refHMM);
  hmm_reference_solution.clear();

  std::string location_hmm_ref = "data/HMM/Model_Problem_1/Macro_10_Micro_8/hmm_solution_discFunc_refLevel_10";
  bool hmm_ref_reader_is_open = false;
  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_hmm_ref( (location_hmm_ref).c_str() );
  discrete_function_reader_hmm_ref.open();
  discrete_function_reader_hmm_ref.read(0, hmm_reference_solution);
  DSC_LOG_INFO << "HMM reference read." << std::endl;
  #endif   // HMM_REFERENCE

  // ! ************************* Assembling and solving the HMM problem ****************************

  #ifdef ADAPTIVE
  // number of the loop cycle of the while-loop
  int loop_cycle = 1;
  double total_hmm_time = 0.0;
  bool repeat = true;
  while (repeat == true)
  {
    if ( data_file.is_open() )
    {
      data_file << "########################### LOOP CYCLE " << loop_cycle << " ###########################"
                << std::endl << std::endl << std::endl;
    }
  #endif   // ADAPTIVE

  DSC_LOG_INFO << std::endl << "Solving HMM-macro-problem for " << discreteFunctionSpace.size()
            << " unkowns and polynomial order "
            << HMM::DiscreteFunctionSpaceType::polynomialOrder << "."
            << std::endl << std::endl;

  // ----------------------------------------------------------------------------------------------//
  // ----------------------- THE DISCRETE HMM OPERATOR -----------------------------------//
  // ----------------------------------------------------------------------------------------------//

  // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain cell problem
  typename HMM::CellProblemNumberingManagerType cp_num_manager(discreteFunctionSpace);

  // ! define the elliptic hmm operator that describes our 'homogenized' macro problem
  // ( effect of the elliptic hmm operator on a certain discrete function )
  typename HMM::EllipticHMMOperatorType discrete_elliptic_hmm_op(discreteFunctionSpace,
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

  // ! solution vector
  // solution of the heterogeneous multiscale finite element method, where we used the Newton method to solve the
  // non-linear system of equations
  typename HMM::DiscreteFunctionType hmm_solution(filename + " HMM (+Newton) Solution", discreteFunctionSpace);
  hmm_solution.clear();

  #ifdef ADAPTIVE
  // one fot the discreteFunctionSpace
  RestrictProlongOperatorType rp(hmm_solution);
  AdaptationManagerType adaptationManager(grid, rp);
  #endif // ifdef ADAPTIVE

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
    #ifdef STOCHASTIC_PERTURBATION
    double hmm_tolerance = 1e-01 * VARIANCE;
    #else
    double hmm_tolerance = 1e-05;
    #endif // ifdef STOCHASTIC_PERTURBATION
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

  // ! ---------- Error Estimation ----------
  #ifdef ERRORESTIMATION
  {
    // Notation: u_H = hmm_solution
    // to load the solutions of the cell problems:
    // location of the solutions of the cell problems for the base function set:
    std::string cell_solution_location_baseSet;
    // location of the solutions of the cell problems for the discrete function u_H:
    std::string cell_solution_location_discFunc;

    cell_solution_location_baseSet = "data/HMM/" + filename + "/cell_problems/_cellSolutions_baseSet";
    cell_solution_location_discFunc = "data/HMM/" + filename + "/cell_problems/_cellSolutions_discFunc";

    // reader for the cell problem data file (for tha macro base set):
    DiscreteFunctionReader discrete_function_reader_baseSet( (cell_solution_location_baseSet).c_str() );
    discrete_function_reader_baseSet.open();

    // reader for the cell problem data file (for u_H):
    DiscreteFunctionReader discrete_function_reader_discFunc( (cell_solution_location_discFunc).c_str() );
    discrete_function_reader_discFunc.open();

    ErrorEstimatorType error_estimator(periodicDiscreteFunctionSpace,
                                       discreteFunctionSpace,
                                       auxiliaryDiscreteFunctionSpace,
                                       diffusion_op);

    RangeType estimated_source_error = 0.0;
    RangeType estimated_approximation_error = 0.0;
    RangeType estimated_residual_error = 0.0;

    // estimated_residual_error = sqrt( estimated_residual_error_micro_jumps^2 + estimated_residual_error_macro_jumps^2 )
    RangeType estimated_residual_error_micro_jumps = 0.0;
    RangeType estimated_residual_error_macro_jumps = 0.0;
    #ifdef TFR
    RangeType estimated_tfr_error = 0.0;
    #endif

    RangeType local_error_indicator[discreteFunctionSpace.grid().size(0)];
    RangeType minimal_loc_indicator = 10000.0;
    RangeType maximal_loc_indicator = 0.0;
    RangeType average_loc_indicator = 0.0;

    int element_number = 0;
    typedef DiscreteFunctionSpaceType::IteratorType IteratorType;
    IteratorType endit = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
    {
      const EntityType& entity = *it;
      // corrector of u_H^(n-1) \approx u_H on the macro element T
      PeriodicDiscreteFunctionType corrector_u_H_on_entity("Corrector of u_H", periodicDiscreteFunctionSpace);
      corrector_u_H_on_entity.clear();

      // in the linear case, we still need to compute the corrector of u_H:
      if (DSC_CONFIG.get("problem.linear", true)) {
        PeriodicDiscreteFunctionType corrector_of_base_func("Corrector of macro base function",
                                                            periodicDiscreteFunctionSpace);
        corrector_of_base_func.clear();

        DiscreteFunctionType::LocalFunctionType local_hmm_solution = hmm_solution.localFunction(entity);

        const BaseFunctionSetType& baseSet = discreteFunctionSpace.baseFunctionSet(entity);
        const unsigned int numMacroBaseFunctions = baseSet.numBaseFunctions();
        int cell_problem_id[numMacroBaseFunctions];
        for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
        {
          cell_problem_id[i] = cp_num_manager.get_number_of_cell_problem(it, i);
          discrete_function_reader_baseSet.read(cell_problem_id[i], corrector_of_base_func);
          corrector_of_base_func *= local_hmm_solution[i];
          corrector_u_H_on_entity += corrector_of_base_func;
          corrector_of_base_func.clear();
        }
      } else {
        // in the nonlinear case this corrector is already available
        discrete_function_reader_discFunc.read(element_number, corrector_u_H_on_entity);
      }

      // contribution of the local source error

      RangeType local_source_indicator = error_estimator.indicator_f(f, entity);
      estimated_source_error += local_source_indicator;

      // contribution of the local approximation error
      RangeType local_approximation_indicator = error_estimator.indicator_app_1(entity,
                                                                                hmm_solution,
                                                                                corrector_u_H_on_entity);

      local_approximation_indicator += error_estimator.indicator_app_2(entity, hmm_solution, corrector_u_H_on_entity);
      estimated_approximation_error += local_approximation_indicator;

      // contribution of the local residual error
      RangeType local_residual_indicator = error_estimator.indicator_res_T(entity, hmm_solution, corrector_u_H_on_entity);
      estimated_residual_error_micro_jumps += local_residual_indicator;

      typedef GridPartType::IntersectionIteratorType IntersectionIteratorType;
      IntersectionIteratorType endnit = gridPart.iend(entity);
      for (IntersectionIteratorType nit = gridPart.ibegin(entity); nit != endnit; ++nit)
      {
        if ( nit->neighbor() )           // if there is a neighbor entity
        {
          // corrector of u_H^(n-1) \approx u_H on the neighbor element
          PeriodicDiscreteFunctionType corrector_u_H_on_neighbor_entity("Corrector of u_H", periodicDiscreteFunctionSpace);
          corrector_u_H_on_neighbor_entity.clear();

          EntityPointerType it_outside = nit->outside();
          const EntityType& entity_outside = *it_outside;

          // in the linear case, we still need to compute the corrector of u_H:
          if (DSC_CONFIG.get("problem.linear", true)) {
            PeriodicDiscreteFunctionType corrector_of_base_func_neighbor("Corrector of macro base function",
                                                                         periodicDiscreteFunctionSpace);
            corrector_of_base_func_neighbor.clear();

            DiscreteFunctionType::LocalFunctionType local_hmm_solution_neighbor = hmm_solution.localFunction(entity_outside);

            const BaseFunctionSetType& baseSet_neighbor = discreteFunctionSpace.baseFunctionSet(entity_outside);
            const unsigned int numMacroBaseFunctions_neighbor = baseSet_neighbor.numBaseFunctions();
            int cell_problem_id_neighbor[numMacroBaseFunctions_neighbor];
            for (unsigned int i = 0; i < numMacroBaseFunctions_neighbor; ++i)
            {
              cell_problem_id_neighbor[i] = cp_num_manager.get_number_of_cell_problem(it_outside, i);
              discrete_function_reader_baseSet.read(cell_problem_id_neighbor[i], corrector_of_base_func_neighbor);
              corrector_of_base_func_neighbor *= local_hmm_solution_neighbor[i];
              corrector_u_H_on_neighbor_entity += corrector_of_base_func_neighbor;
              corrector_of_base_func_neighbor.clear();
            }
          } else {
            int neighbor_element_number = cp_num_manager.get_number_of_cell_problem(it_outside);
            // in the nonlinear case this corrector is already available
            discrete_function_reader_discFunc.read(neighbor_element_number, corrector_u_H_on_neighbor_entity);
          }

          RangeType val = error_estimator.indicator_res_E(*nit,
                                                          hmm_solution,
                                                          corrector_u_H_on_entity,
                                                          corrector_u_H_on_neighbor_entity);
          local_residual_indicator += val;
          estimated_residual_error_macro_jumps += val;
        }
      }

      estimated_residual_error += local_residual_indicator;

      #ifdef TFR
      // use 'indicator_effective_tfr' or 'indicator_tfr_1'
      // contribution of the local tfr error:
      RangeType local_tfr_indicator = error_estimator.indicator_tfr_1(entity, hmm_solution, corrector_u_H_on_entity);
      estimated_tfr_error += local_tfr_indicator;
      #endif // ifdef TFR

      local_error_indicator[element_number] = local_source_indicator
                                              + local_approximation_indicator
                                              + local_residual_indicator;

      #ifdef TFR
      local_error_indicator[element_number] += local_tfr_indicator;
      #endif

      if (local_error_indicator[element_number] < minimal_loc_indicator)
      {
        minimal_loc_indicator = local_error_indicator[element_number];
      }

      if (local_error_indicator[element_number] > maximal_loc_indicator)
      {
        maximal_loc_indicator = local_error_indicator[element_number];
      }

      average_loc_indicator += local_error_indicator[element_number];

      element_number += 1;
    }       // endfor

    average_loc_indicator /= discreteFunctionSpace.grid().size(0);

    estimated_source_error = sqrt(estimated_source_error);
    estimated_approximation_error = sqrt(estimated_approximation_error);
    estimated_residual_error = sqrt(estimated_residual_error);

    estimated_residual_error_micro_jumps = sqrt(estimated_residual_error_micro_jumps);
    estimated_residual_error_macro_jumps = sqrt(estimated_residual_error_macro_jumps);

    RangeType estimated_error = estimated_source_error + estimated_approximation_error + estimated_residual_error;

    #ifdef TFR
    estimated_tfr_error = sqrt(estimated_tfr_error);
    estimated_error += estimated_tfr_error;
    #endif // ifdef TFR

    #ifdef ADAPTIVE
    // maximum variation (up) from average
    double max_variation = average_loc_indicator / maximal_loc_indicator;
    double min_variation = average_loc_indicator / minimal_loc_indicator;
    #endif // ifdef ADAPTIVE
           // ! -------- End Error Estimation --------
  }
  #endif   // #ifdef ERRORESTIMATION

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

  RangeType hmm_error = impL2error.norm_adaptive_grids_2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(
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

  RangeType hmm_vs_hmm_ref_error = impL2error.norm_adaptive_grids_2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(
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
  RangeType hom_error = l2error.norm2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(homogenized_solution,
                                                                                           fem_newton_solution);

  DSC_LOG_INFO << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std::endl << std::endl;
  if ( data_file.is_open() )
  {
    data_file << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std::endl;
  }
  #endif // ifdef FINE_SCALE_REFERENCE

  RangeType hom_hmm_error = l2error.norm2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(hmm_solution,
                                                                                               homogenized_solution);

  DSC_LOG_INFO << "|| u_hom - u_hmm ||_L2 =  " << hom_hmm_error << std::endl << std::endl;
  if ( data_file.is_open() )
  {
    data_file << "|| u_hom - u_hmm ||_L2 =  " << hom_hmm_error << std::endl;
  }
  #endif // ifdef HOMOGENIZEDSOL_AVAILABLE

  if (ProblemDataType::has_exact_solution)
  {
    typename HMM::ExactSolutionType u;
    typename HMM::RangeType exact_hmm_error = l2error.template norm< typename HMM::ExactSolutionType >(u,
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

  if (ProblemDataType::has_exact_solution) {
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

  #ifdef ADAPTIVE
  int default_refinement = 0;

  // double error_tolerance_ = 0.2;

  // int(...) rundet ab zum naechsten Integer
  // assuming we had a quadratic order of convergence of the error estimator and
  // that we have a certain estimated error for the current uniform grid, than we can compute how many additional
  // uniform refinements are required to get under the error (estimator) tolerance:
  // number_of_uniform_refinments = int( sqrt( (\eta_have) / (\eta_want) ) )
  // (obtained from the EOC formula)
  // uniform contribution only for the first loop cycle
  if (loop_cycle == 1)
  {
    // "divided by 2.0" we go half the way with a uniform computation
    int number_of_uniform_refinements = 2 * int(int( sqrt(estimated_error / error_tolerance_) ) / 2.0);

    if ( data_file.is_open() )
    {
      data_file << std::endl << "Uniform default refinement:" << std::endl << std::endl;
      data_file << "sqrt( estimated_error / error_tolerance_ ) = "
                << sqrt(estimated_error / error_tolerance_) << std::endl;
      data_file << "number_of_uniform_refinements = " << number_of_uniform_refinements << std::endl;
      data_file << "***************" << std::endl << std::endl;
    }

    default_refinement = number_of_uniform_refinements;
  }

  int number_of_areas;
  if (loop_cycle == 1)
  {
    number_of_areas = 1;
  } else {
    if (loop_cycle == 3)
    {
      number_of_areas = 2;
    } else {
      number_of_areas = 2;
    }
  }

  double border[number_of_areas - 1];
  border[0] = 0.5;
  for (int bo = 1; bo < (number_of_areas - 1); ++bo)
  {
    border[bo] = border[bo - 1] + ( (1.0 - border[bo - 1]) / 2.0 );
  }

  // 3 areas: 1: |0-30%| 2: |30-80%| 3: |80-100%|
  // border[0] = 0.3;
  // border[1] = 0.8;
  // border[2] = 0.95;

  int refinements_in_area[number_of_areas];
  for (int bo = 0; bo < number_of_areas; ++bo)
  {
    refinements_in_area[bo] = default_refinement + bo + 1;
  }

  if ( data_file.is_open() )
  {
    data_file << "Adaption strategy:" << std::endl << std::endl;
    data_file
    << "Define 'variation = (indicator_on_element - average_indicator) / (maximum_indicator - average_indicator)'"
    << std::endl;
    data_file << "Subdivide the region [average_indicator,maximum_indicator] into " << number_of_areas << " areas."
              << std::endl;
    if (number_of_areas == 1)
    {
      data_file << "1.: [average_indicator,maximum_indicator]. Mark elements for " << refinements_in_area[0]
                << " refinements." << std::endl;
    } else {
      data_file << "1.: [average_indicator," << border[0]
                << "*maximum_indicator]. If 'variance' in area: mark elements for " << refinements_in_area[0]
                << " refinements." << std::endl;
      for (int bo = 1; bo < (number_of_areas - 1); ++bo)
        data_file << bo + 1 << ".: ["
                  << border[bo
                  - 1] << "*average_indicator," << border[bo]
                  << "*maximum_indicator]. If 'variance' in area: mark elements for "
                  << refinements_in_area[bo] << " refinements." << std::endl;
      data_file << number_of_areas << ".: ["
                << border[number_of_areas
                - 2] << "*average_indicator,maximum_indicator]. If 'variance' in area: mark elements for "
                << refinements_in_area[number_of_areas - 1] << " refinements." << std::endl;
    }
    data_file << "Default refinement for elements with 'variance <= 0 ': " << default_refinement << std::endl;
  }

  if (estimated_error < error_tolerance_)
  {
    repeat = false;
    DSC_LOG_INFO << "Total HMM time = " << total_hmm_time << "s." << std::endl;
    data_file << std::endl << std::endl << "Total HMM time = " << total_hmm_time << "s." << std::endl << std::endl;
  } else {
    element_number = 0;
    typedef DiscreteFunctionSpaceType::IteratorType IteratorType;
    IteratorType endit_test = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit_test; ++it)
    {
      int additional_refinement;
      if (local_error_indicator[element_number] <= average_loc_indicator)
      {
        additional_refinement = default_refinement;
      } else {
        double variation = (local_error_indicator[element_number] - average_loc_indicator)
                           / (maximal_loc_indicator - average_loc_indicator);

        if (number_of_areas == 1)
        {
          additional_refinement = refinements_in_area[0];
        } else {
          if (variation <= border[0])
          {
            additional_refinement = refinements_in_area[0];
          }

          for (int bo = 1; bo <= (number_of_areas - 2); ++bo)
            if ( (variation > border[bo - 1]) && (variation <= border[bo]) )
            {
              additional_refinement = refinements_in_area[bo];
            }


          if (variation > border[number_of_areas - 2])
          {
            additional_refinement = refinements_in_area[number_of_areas - 1];
          }
        }
      }

      grid.mark(additional_refinement, *it);
      element_number += 1;
    }
    adaptationManager.adapt();
  }

  if ( data_file.is_open() )
  {
    data_file << std::endl << "#########################################################################"
              << std::endl << std::endl << std::endl;
  }

  loop_cycle += 1;
}     // end of repeat loop (for the adaptive cicles)

  #endif   // #ifdef ADAPTIVE
}

#endif // DUNE_MS_HMM_ALGORITHM_HH
