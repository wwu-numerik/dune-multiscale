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
#include <dune/multiscale/tools/homogenizer/elliptic_homogenizer.hh>
#include <dune/stuff/common/logging.hh>
#include "algorithm_step.hh"

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

template <class HMM>
void fsr_compute(typename HMM::DiscreteFunctionType& fem_newton_solution,
                 const typename HMM::DiscreteFunctionSpaceType& finerDiscreteFunctionSpace,
                 const typename HMM::EllipticOperatorType& discrete_elliptic_op,
                 const std::string& filename,
                 const Dune::RightHandSideAssembler< typename HMM::DiscreteFunctionType >& rhsassembler)
{
  static const int hmm_polorder = 2* HMM::DiscreteFunctionSpaceType::polynomialOrder + 2;

  //! *************************** Assembling the reference problem ****************************
  // ( fine scale reference solution = fem_newton_solution )

  //! (stiffness) matrix
  typename HMM::FEMMatrix fem_newton_matrix("FEM Newton stiffness matrix", finerDiscreteFunctionSpace, finerDiscreteFunctionSpace);

  //! right hand side vector
  // right hand side for the finite element method with Newton solver:
  // ( also right hand side for the finer discrete function space )
  typename HMM::DiscreteFunctionType fem_newton_rhs("fem newton rhs", finerDiscreteFunctionSpace);
  fem_newton_rhs.clear();

  const typename HMM::FirstSourceType f;   // standard source f

  if (DSC_CONFIG_GET("problem.linear", true))
  {
    DSC_LOG_INFO << "Solving linear problem." << std::endl;
    DSC_LOG_INFO << "Solving linear problem with standard FEM and resolution level "
              << typename HMM::ModelProblemDataType().getRefinementLevelReferenceProblem() << "." << std::endl;
    DSC_LOG_INFO << "------------------------------------------------------------------------------" << std::endl;

    // to assemble the computational time
    Dune::Timer assembleTimer;

    // assemble the stiffness matrix
    discrete_elliptic_op.assemble_matrix(fem_newton_matrix);

    DSC_LOG_INFO << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

    // assemble right hand side
    rhsassembler.template assemble< hmm_polorder >(f, fem_newton_rhs);

    // set Dirichlet Boundary to zero
    boundaryTreatment(fem_newton_rhs);

    const typename HMM::InverseFEMMatrix fem_biCGStab(fem_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE);
    fem_biCGStab(fem_newton_rhs, fem_newton_solution);

    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
    DSC_LOG_INFO << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl
              << std::endl;
  } else {
    DSC_LOG_INFO << "Solving non-linear problem." << std::endl;
    DSC_LOG_INFO << "Solving nonlinear problem with FEM + Newton-Method. Resolution level of grid = "
              << typename HMM::ModelProblemDataType().getRefinementLevelReferenceProblem() << "." << std::endl;
    DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;

    Dune::Timer assembleTimer;
    //! residual vector
    // current residual
    typename HMM::DiscreteFunctionType fem_newton_residual(filename + "FEM Newton Residual", finerDiscreteFunctionSpace);
    fem_newton_residual.clear();

    typename HMM::RangeType relative_newton_error_finescale = 10000.0;
    typename HMM::RangeType rhs_L2_norm = 10000.0;

    int iteration_step = 1;
    // the Newton step for the FEM reference problem (solved with Newton Method):
    // L2-Norm of residual < tolerance ?
    double tolerance = 1e-06;
    while (relative_newton_error_finescale > tolerance)
    {
      // (here: fem_newton_solution = solution from the last iteration step)
      DSC_LOG_INFO << "Newton iteration " << iteration_step << ":" << std::endl;
      Dune::Timer stepAssembleTimer;
      // assemble the stiffness matrix
      discrete_elliptic_op.assemble_jacobian_matrix(fem_newton_solution, fem_newton_matrix);

      DSC_LOG_INFO << "Time to assemble FEM Newton stiffness matrix for current iteration: "
                   << stepAssembleTimer.elapsed() << "s" << std::endl;

      // assemble right hand side
      const typename HMM::DiffusionType diffusion_op;
      rhsassembler.template assemble_for_Newton_method< hmm_polorder >(f,
                                                                       diffusion_op,
                                                                       fem_newton_solution,
                                                                       fem_newton_rhs);

      const Dune::L2Norm< typename HMM::DiscreteFunctionType::GridPartType > l2norm(fem_newton_rhs.gridPart());
      rhs_L2_norm = l2norm.norm(fem_newton_rhs);
      if (rhs_L2_norm < 1e-10)
      {
        // residual solution almost identical to zero: break
        DSC_LOG_INFO << "Residual solution almost identical to zero. Therefore: break loop." << std::endl;
        DSC_LOG_INFO << "(L^2-Norm of current right hand side = " << rhs_L2_norm << " < 1e-10)" << std::endl;
        break;
      }
      // set Dirichlet Boundary to zero
      boundaryTreatment(fem_newton_rhs);

      const typename HMM::InverseFEMMatrix fem_newton_biCGStab(fem_newton_matrix, 1e-8, 1e-8, 20000, true);
      fem_newton_biCGStab(fem_newton_rhs, fem_newton_residual);

      if ( fem_newton_residual.dofsValid() )
      {
        fem_newton_solution += fem_newton_residual;
        relative_newton_error_finescale = l2norm.norm(fem_newton_residual);
        relative_newton_error_finescale /= l2norm.norm(fem_newton_solution);

        DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
        // residual solution almost identical to zero: break
        DSC_LOG_INFO << "Relative L2-Newton Error = " << relative_newton_error_finescale << std::endl;
        if (relative_newton_error_finescale <= tolerance)
        {
          DSC_LOG_INFO << "Since tolerance = " << tolerance << ": break loop." << std::endl;
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
    DSC_LOG_INFO << seperator_line;
    DSC_LOG_INFO << "Problem with FEM + Newton-Method solved in " << assembleTimer.elapsed() << "s." << std::endl
              << std::endl << std::endl;
  }// end 'problem.linear <-> else'

  //! ********************** End of assembling the reference problem ***************************
}

template < class HMM >
void fsr_load(typename HMM::DiscreteFunctionType& fem_newton_solution)
{
  fem_newton_solution.clear();
  const boost::filesystem::path modeprob_s(
      (boost::format("Model_Problem_%d") % DSC_CONFIG_GET("problem.numbder", 0u)).str());

  const int refinement_level_referenceprob_ =
      typename HMM::ModelProblemDataType().getRefinementLevelReferenceProblem();
  const boost::filesystem::path reference_solution_directory_s(
        (boost::format("reference_solution_ref_%d") % refinement_level_referenceprob_).str());

  const std::string reference_solution_name_s(
        (boost::format("finescale_solution_discFunc_refLevel_%d") % refinement_level_referenceprob_).str());

  boost::filesystem::path location_fine_scale_ref("HMM");
  location_fine_scale_ref /= modeprob_s
                             / reference_solution_directory_s
                             / reference_solution_name_s;

  // reader for the cell problem data file:
  DiscreteFunctionReader(location_fine_scale_ref).read(0, fem_newton_solution);
  DSC_LOG_INFO << "fine scale reference read." << std::endl;
}


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
void print_info(const ProblemDataType& info, std::ostream& out)
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
  const int refinement_level_macrogrid_ = DSC_CONFIG_GET("grid.refinement_level_macrogrid", 0);
  out << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << "." << std::endl
            << std::endl;
  if (DSC_CONFIG_GET("problem.linear", true))
    out << "Problem is declared as being LINEAR." << std::endl;
  else
    out << "Problem is declared as being NONLINEAR." << std::endl;

  if (ProblemDataType::has_exact_solution) {
    out << "Exact solution is available." << std::endl << std::endl;
  } else {
    out << "Exact solution is not available." << std::endl << std::endl;
  }
  out << "Computations were made for:" << std::endl << std::endl;
  out << "Refinement Level for (uniform) Macro Grid = " << refinement_level_macrogrid_ << std::endl;
  const int refinement_level_cellgrid = DSC_CONFIG_GET("grid.refinement_level_cellgrid", 1);
  out << "Refinement Level for Periodic Micro Grid = " << refinement_level_cellgrid << std::endl << std::endl;
  #ifdef TFR
  out << "We use TFR-HMM (HMM with test function reconstruction)." << std::endl;
  #else
  out << "We use HMM without test function reconstruction (NO TFR)." << std::endl;
  #endif // ifdef TFR
  if (DSC_CONFIG_GET("AD_HOC_COMPUTATION", false)) {
    out << "Cell problems are solved ad hoc (where required)." << std::endl << std::endl;
  } else {
    out << "Cell problems are solved and saved (in a pre-process)." << std::endl << std::endl;
    #ifdef ERRORESTIMATION
    out << "Error estimation activated!" << std::endl << std::endl;
    #endif
  }
  out << "Epsilon = " << epsilon_ << std::endl;
  out << "Estimated Epsilon = " << epsilon_est_ << std::endl;
  out << "Delta (edge length of cell-cube) = " << delta_ << std::endl;
  if (DSC_CONFIG_GET("problem.stochastic_pertubation", false))
    out << std::endl << "Stochastic perturbation added. Variance = " << DSC_CONFIG_GET("problem.stochastic_variance", 0.01) << std::endl;
  if (DSC_CONFIG_GET("hmm.adaptive", true)) {
    //only used in adaptive config
    const double error_tolerance_ = DSC_CONFIG_GET("problem.error_tolerance", 1e-6);
    out << std::endl << "Adaptive computation. Global error tolerance for program abort = "
              << error_tolerance_ << std::endl;
  }
  out << std::endl << std::endl;
}

template < class HMM >
bool adapt(const HMMResult<HMM>& result,
           const int loop_cycle,
           const double error_tolerance_,
           const typename HMM::DiscreteFunctionSpaceType& discreteFunctionSpace,
           typename HMM::AdaptationManagerType& adaptationManager
           )
{
  bool repeat = true;
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
    const int number_of_uniform_refinements = 2 * int(int( sqrt(result.estimated_error / error_tolerance_) ) / 2.0);
    DSC_LOG_INFO << std::endl << "Uniform default refinement:" << std::endl << std::endl;
    DSC_LOG_INFO << "sqrt( estimated_error / error_tolerance_ ) = "
              << sqrt(result.estimated_error / error_tolerance_) << std::endl;
    DSC_LOG_INFO << "number_of_uniform_refinements = " << number_of_uniform_refinements << std::endl;
    DSC_LOG_INFO << "***************" << std::endl << std::endl;

    default_refinement = number_of_uniform_refinements;
  }

  const int number_of_areas = (loop_cycle == 1) ? 1 : 2;

  std::vector<double> border(number_of_areas - 1);
  border[0] = 0.5;
  for (int bo = 1; bo < (number_of_areas - 1); ++bo)
  {
    border[bo] = border[bo - 1] + ( (1.0 - border[bo - 1]) / 2.0 );
  }

  // 3 areas: 1: |0-30%| 2: |30-80%| 3: |80-100%|
  // border[0] = 0.3;
  // border[1] = 0.8;
  // border[2] = 0.95;

  std::vector<int> refinements_in_area(number_of_areas);
  for (int bo = 0; bo < number_of_areas; ++bo)
  {
    refinements_in_area[bo] = default_refinement + bo + 1;
  }

  DSC_LOG_INFO << "Adaption strategy:" << std::endl << std::endl;
  DSC_LOG_INFO << "Define 'variation = (indicator_on_element - average_indicator) / (maximum_indicator - average_indicator)'"
               << std::endl;
  DSC_LOG_INFO << "Subdivide the region [average_indicator,maximum_indicator] into " << number_of_areas << " areas."
               << std::endl;
  if (number_of_areas == 1)
  {
    DSC_LOG_INFO << "1.: [average_indicator,maximum_indicator]. Mark elements for " << refinements_in_area[0]
              << " refinements." << std::endl;
  } else {
    DSC_LOG_INFO << "1.: [average_indicator," << border[0]
              << "*maximum_indicator]. If 'variance' in area: mark elements for " << refinements_in_area[0]
              << " refinements." << std::endl;
    for (int bo = 1; bo < (number_of_areas - 1); ++bo)
      DSC_LOG_INFO << bo + 1 << ".: ["
                << border[bo
                - 1] << "*average_indicator," << border[bo]
                << "*maximum_indicator]. If 'variance' in area: mark elements for "
                << refinements_in_area[bo] << " refinements." << std::endl;
    DSC_LOG_INFO << number_of_areas << ".: ["
              << border[number_of_areas
              - 2] << "*average_indicator,maximum_indicator]. If 'variance' in area: mark elements for "
              << refinements_in_area[number_of_areas - 1] << " refinements." << std::endl;
  }
  DSC_LOG_INFO << "Default refinement for elements with 'variance <= 0 ': " << default_refinement << std::endl;


  if (result.estimated_error < error_tolerance_)
  {
    repeat = false;
    //!TODO profiler
//    DSC_LOG_INFO << "Total HMM time = " << total_hmm_time << "s." << std::endl;
//    DSC_LOG_INFO << std::endl << std::endl << "Total HMM time = " << total_hmm_time << "s." << std::endl << std::endl;
  } else {

    int element_number = 0;
    typedef typename HMM::DiscreteFunctionSpaceType::IteratorType IteratorType;
    const IteratorType endit_test = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit_test; ++it)
    {
      int additional_refinement;
      if (result.local_error_indicator[element_number] <= result.average_loc_indicator)
      {
        additional_refinement = default_refinement;
      } else
      {
        const double variation = (result.local_error_indicator[element_number] - result.average_loc_indicator)
                           / (result.maximal_loc_indicator - result.average_loc_indicator);

        if (number_of_areas == 1)
        {
          additional_refinement = refinements_in_area[0];
        } else
        {
          if (variation <= border[0])
          {
            additional_refinement = refinements_in_area[0];
          }
          for (int bo = 1; bo <= (number_of_areas - 2); ++bo) {
            if ( (variation > border[bo - 1]) && (variation <= border[bo]) )
            {
              additional_refinement = refinements_in_area[bo];
            }
          }
          if (variation > border[number_of_areas - 2])
          {
            additional_refinement = refinements_in_area[number_of_areas - 1];
          }
        }
      }
      discreteFunctionSpace.gridPart().grid().mark(additional_refinement, *it);
      element_number += 1;
    }
    adaptationManager.adapt();
  }
  return repeat;
}

//! the main hmm computation
template < class HMMTraits >
void algorithm(const std::string& /*UnitCubeName*/,
               typename HMMTraits::GridPointerType& macro_grid_pointer,   // grid pointer that belongs to the macro grid
               typename HMMTraits::GridPointerType& fine_macro_grid_pointer,   // grid pointer that belongs to the fine macro grid (for
                                                           // reference computations)
               typename HMMTraits::GridPointerType& periodic_grid_pointer,   // grid pointer that belongs to the periodic micro grid
               int /*refinement_difference*/,   // refinement difference for the macro grid (problem-to-solve vs. reference
                                            // problem)
               const std::string filename) {
  typedef HMMTraits HMM;
  using namespace Dune;

  const typename HMM::ModelProblemDataType problem_data;
  print_info(problem_data, DSC_LOG_INFO);
  //! ---- tools ----
  // model problem data
// UNUSED  Problem::ModelProblemData problem_data;
// set of hmm parameters/information
// UNUSED  Multiscale::HMMParameters method_info;

  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for HMM-macro-problem
  typename HMM::GridPartType gridPart(*macro_grid_pointer);
  // grid part for the periodic function space, required for HMM-cell-problems
  typename HMM::PeriodicGridPartType periodicGridPart(*periodic_grid_pointer);
  // grid part for the global function space, required for the detailed fine-scale computation (very high resolution)
  typename HMM::GridPartType gridPartFine(*fine_macro_grid_pointer);

  typename HMM::GridType& grid = gridPart.grid();
  typename HMM::GridType& gridFine = gridPartFine.grid();
  //! --------------------------------------------------------------------------------------

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  typename HMM::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  // the global-problem function space for the reference computation:
  typename HMM::DiscreteFunctionSpaceType finerDiscreteFunctionSpace(gridPartFine);
  // the local-problem function space (containing periodic functions):
  typename HMM::PeriodicDiscreteFunctionSpaceType periodicDiscreteFunctionSpace(periodicGridPart);
  //! --------------------------------------------------------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const typename HMM::DiffusionType diffusion_op;

  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  Dune::RightHandSideAssembler< typename HMM::DiscreteFunctionType > rhsassembler;
  const typename HMM::FirstSourceType f;   // standard source f

  //! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
 const typename HMM::EllipticOperatorType discrete_elliptic_op(finerDiscreteFunctionSpace, diffusion_op);

// UNUSED  RangeType size_of_domain = get_size_of_domain(discreteFunctionSpace);
  static const int hmm_polorder = 2* HMM::DiscreteFunctionSpaceType::polynomialOrder + 2;
  if (DSC_CONFIG_GET("HOMOGENIZEDSOL_AVAILABLE", false)) {
    if (DSC_CONFIG_GET("problem.linear", true)) {
      const std::string unit_cell_location = "../dune/multiscale/grids/cell_grids/unit_cube.dgf";
      // descretized homogenizer:

      const typename HMM::HomogenizerType disc_homogenizer(unit_cell_location);
      const typename HMM::HomogenizerType::HomTensorType A_hom = disc_homogenizer.getHomTensor(diffusion_op);
      const typename HMM::HomDiffusionType hom_diffusion_op(A_hom);

      //!TODO check: hatte nur 2 tmp parameter, Masse hinzugefUGT
      typedef DiscreteEllipticOperator< typename HMM::DiscreteFunctionType,
                                        typename HMM::HomDiffusionType, typename HMM::MassTermType > HomEllipticOperatorType;

      HomEllipticOperatorType hom_discrete_elliptic_op(finerDiscreteFunctionSpace, hom_diffusion_op);

      typename HMM::FEMMatrix hom_stiff_matrix("homogenized stiffness matrix", finerDiscreteFunctionSpace, finerDiscreteFunctionSpace);

      typename HMM::DiscreteFunctionType hom_rhs("homogenized rhs", finerDiscreteFunctionSpace);
      hom_rhs.clear();

      typename HMM::DiscreteFunctionType homogenized_solution(filename + " Homogenized Solution", finerDiscreteFunctionSpace);
      homogenized_solution.clear();
      hom_discrete_elliptic_op.assemble_matrix(hom_stiff_matrix);

      rhsassembler.template assemble < hmm_polorder >(f, hom_rhs);

      // set Dirichlet Boundary to zero
      boundaryTreatment(hom_rhs);

      const typename HMM::InverseFEMMatrix hom_biCGStab(hom_stiff_matrix, 1e-8, 1e-8, 20000, VERBOSE);
      hom_biCGStab(hom_rhs, homogenized_solution);
    } else {
      //!TODO nonlinear Homogenizer??
    }
  }

  //! solution vector
  // - By fem_newton_solution, we denote the "fine scale reference solution" (used for comparison)
  //    ( if the elliptic problem is linear, the 'fem_newton_solution' is determined without the Newton method )
  // - solution of the finite element method, where we used the Newton method to solve the non-linear system of equations
  //   in general this will be an accurate approximation of the exact solution, that is why we it also called reference solution
  typename HMM::DiscreteFunctionType fem_newton_solution(filename + " Reference (FEM Newton) Solution", finerDiscreteFunctionSpace);
  fem_newton_solution.clear();

  if (DSC_CONFIG_GET("fsr", true))
  {
    if (DSC_CONFIG_GET("fsr_compute", true))
    {
      fsr_compute<HMM>(fem_newton_solution, finerDiscreteFunctionSpace,
                       discrete_elliptic_op, filename, rhsassembler);
    }
    //! load und compute sollten sich ausschliessen??
    if (DSC_CONFIG_GET("fsr_load", false))
    {
      fsr_load<HMM>(fem_newton_solution);
    }
  }

  // noch per Hand die Daten eingetragen:
  const std::string macroGridName_refHMM = problem_data.getMacroGridFile();
  DSC_LOG_INFO << "loading dgf: " << macroGridName_refHMM << std::endl;
  typename HMM::GridPointerType macro_grid_pointer_refHMM(macroGridName_refHMM);
  typename HMM::GridPartType gridPart_refHMM(*macro_grid_pointer_refHMM);
  typename HMM::DiscreteFunctionSpaceType discreteFunctionSpace_refHMM(gridPart_refHMM);
  typename HMM::DiscreteFunctionType hmm_reference_solution(filename + " Reference (HMM) Solution", discreteFunctionSpace_refHMM);
  if (DSC_CONFIG_GET("hmm_reference", false))
  {
    //! Macro_'gridLevel'
    const int gridLevel_refHMM = DSC_CONFIG_GET("grid.ref_hmm", 10);
    macro_grid_pointer_refHMM->globalRefine(gridLevel_refHMM);
    hmm_reference_solution.clear();
    //!TODO this being hardcoded makes no sense
    const std::string location_hmm_ref = "HMM/Model_Problem_1/Macro_10_Micro_8/hmm_solution_discFunc_refLevel_10";
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader_hmm_ref(location_hmm_ref);
    discrete_function_reader_hmm_ref.read(0, hmm_reference_solution);
    DSC_LOG_INFO << "HMM reference read." << std::endl;
  }

  //! ************************* Assembling and solving the HMM problem ****************************
  // number of the loop cycle of the while-loop
  int loop_cycle = 1;
  const double error_tolerance_ = DSC_CONFIG_GET("problem.error_tolerance", 1e-6);
  bool repeat = true;
  while (repeat == true)
  {
    DSC_LOG_INFO << "########################### LOOP CYCLE " << loop_cycle << " ###########################"
                << std::endl << std::endl << std::endl;

    //! solution vector
    // solution of the heterogeneous multiscale finite element method, where we used the Newton method to solve the
    // non-linear system of equations
    typename HMM::DiscreteFunctionType hmm_solution(filename + " HMM (+Newton) Solution", discreteFunctionSpace);
    hmm_solution.clear();

    typename HMM::RestrictProlongOperatorType rp(hmm_solution);
    typename HMM::AdaptationManagerType adaptationManager(grid, rp);
    const auto result = single_step<HMM>(gridPart, gridPartFine, discreteFunctionSpace, periodicDiscreteFunctionSpace,
                diffusion_op, rhsassembler, filename, hmm_solution, fem_newton_solution);

    if (!DSC_CONFIG_GET("hmm.adaptive", true))
      break;

    if (!adapt<HMM>(result, loop_cycle, error_tolerance_, discreteFunctionSpace, adaptationManager))
      break;

    DSC_LOG_INFO << std::endl << "#########################################################################"
              << std::endl << std::endl << std::endl;
    loop_cycle += 1;
  } // end while (repeat) of repeat loop (for the adaptive cicles)
}

#endif // DUNE_MS_HMM_ALGORITHM_HH
