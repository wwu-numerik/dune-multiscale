#include <config.h>
#include <config.h>
#include <dune/multiscale/hmm/algorithm.hh>

#include <dune/multiscale/hmm/algorithm_step.hh>
#include <dune/multiscale/hmm/result.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/hmm/cell_problem_solver.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/fem/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/problems/selector.hh>

#include <dune/stuff/grid/information.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

//! \todo DOCME
template <class Stream, class DiscFunc>
void oneLinePrint(Stream& stream, const DiscFunc& func) {
  auto it = func.dbegin();
  stream << "\n" << func.name() << ": [ ";
  for (; it != func.dend(); ++it)
    stream << std::setw(5) << *it << "  ";

  stream << " ] " << std::endl;
} // oneLinePrint

//! loads a reference solution from disk
void load_reference(typename CommonTraits::DiscreteFunctionType& reference_solution) {
  reference_solution.clear();

  const std::string location_reference_solution = DSC_CONFIG_GET("problem.rs_file", "");
  if (location_reference_solution == "")
    DUNE_THROW(Dune::InvalidStateException, "Key 'problem.reference_solution' is activated, but no location for "
                                            "reference solution specified. Please specify 'problem.rs_file' or "
                                            "deactivate 'problem.reference_solution'.");

  boost::filesystem::path ref_sol_file(location_reference_solution);

  // reader for the cell problem data file:
  DiscreteFunctionReader(ref_sol_file).read(0, reference_solution);
  DSC_LOG_INFO << "Reference solution successfully read from file." << std::endl;
}

//! \todo replace me with Stuff::something
template <class DiscreteFunctionSpaceType>
typename DiscreteFunctionSpaceType::RangeType
get_size_of_domain(const DiscreteFunctionSpaceType& discreteFunctionSpace) {
  return DSG::dimensions(discreteFunctionSpace.gridPart().grid()).sum();
} // get_size_of_domain

//! outputs Problem info to output stream
void print_info(const CommonTraits::ModelProblemDataType& info, std::ostream& out) {
  // epsilon is specified in the parameter file
  // 'epsilon' in for instance A^{epsilon}(x) = A(x,x/epsilon)
  const double epsilon_ = DSC_CONFIG_GET("problem.epsilon", 1.0f);
  // estimated epsilon (specified in ModelProblemData)
  // estimated epsilon in case epsilon is unknown
  const double epsilon_est_ = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);
  // edge length of the cells in the cells, belonging to the cell problems
  // note that (delta/epsilon_est) needs to be a positive integer!
  // edge length of the cells in the cell proplems,
  const double delta_ = DSC_CONFIG_GET("hmm.delta", 1.0f);
  const int refinement_level_macrogrid_ = DSC_CONFIG_GET("hmm.coarse_grid_level", 0);
  out << "Error File for Elliptic Model Problem " << Problem::name() << "." << std::endl << std::endl;
  if (DSC_CONFIG_GET("problem.linear", true))
    out << "Problem is declared as being LINEAR." << std::endl;
  else
    out << "Problem is declared as being NONLINEAR." << std::endl;

  if (info.hasExactSolution()) {
    out << "Exact solution is available." << std::endl << std::endl;
  } else {
    out << "Exact solution is not available." << std::endl << std::endl;
  }
  out << "Computations were made for:" << std::endl << std::endl;
  out << "Refinement Level for (uniform) Macro Grid = " << refinement_level_macrogrid_ << std::endl;
  const int refinement_level_cellgrid = DSC_CONFIG_GET("hmm.cell_grid_level", 1);
  out << "Refinement Level for Periodic Micro Grid = " << refinement_level_cellgrid << std::endl << std::endl;
  if (DSC_CONFIG_GET("hmm.petrov_galerkin", true))
    out << "We use the HMM in Petrov Galerkin formulation." << std::endl;
  else
    out << "We use the HMM in classical 'symmetric' formulation (non-Petrov-Galerkin)." << std::endl;

  out << "Cell problems are solved and saved (in a pre-process)." << std::endl << std::endl;
  if (DSC_CONFIG_GET("hmm.error_estimation", false))
    out << "Error estimation activated!" << std::endl << std::endl;
  out << "Epsilon = " << epsilon_ << std::endl;
  out << "Estimated Epsilon = " << epsilon_est_ << std::endl;
  out << "Delta (edge length of cell-cube) = " << delta_ << std::endl;
  if (DSC_CONFIG_GET("problem.stochastic_pertubation", false))
    out << std::endl << "Stochastic perturbation added. Variance = " << DSC_CONFIG_GET("problem.stochastic_variance",
                                                                                       0.01) << std::endl;
  if (DSC_CONFIG_GET("hmm.adaptivity", false)) {
    // only used in adaptive config
    const double error_tolerance_ = DSC_CONFIG_GET("hmm.error_tolerance", 1e-6);
    out << std::endl << "Adaptive computation. Global error tolerance for program abort = " << error_tolerance_
        << std::endl;
  }
  out << std::endl << std::endl;
}

//! \TODO docme
bool adapt(const HMMResult& result, const int loop_cycle, const double error_tolerance_,
           const typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
           typename CommonTraits::AdaptationManagerType& adaptationManager) {
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
  if (loop_cycle == 1) {
    // "divided by 2.0" we go half the way with a uniform computation
    const int number_of_uniform_refinements = 2 * int(int(sqrt(result.estimated_error / error_tolerance_)) / 2.0);
    DSC_LOG_INFO << std::endl << "Uniform default refinement:" << std::endl << std::endl;
    DSC_LOG_INFO << "sqrt( estimated_error / error_tolerance_ ) = " << sqrt(result.estimated_error / error_tolerance_)
                 << std::endl;
    DSC_LOG_INFO << "number_of_uniform_refinements = " << number_of_uniform_refinements << std::endl;
    DSC_LOG_INFO << "***************" << std::endl << std::endl;

    default_refinement = number_of_uniform_refinements;
  }

  const int number_of_areas = (loop_cycle == 1) ? 1 : 2;

  std::vector<double> border(number_of_areas);
  border[0] = 0.5;
  for (std::size_t bo = 1; bo < border.size(); ++bo) {
    border[bo] = border[bo - 1] + ((1.0 - border[bo - 1]) / 2.0);
  }

  // 3 areas: 1: |0-30%| 2: |30-80%| 3: |80-100%|
  // border[0] = 0.3;
  // border[1] = 0.8;
  // border[2] = 0.95;

  std::vector<int> refinements_in_area(number_of_areas);
  for (int bo = 0; bo < number_of_areas; ++bo) {
    refinements_in_area[bo] = default_refinement + bo + 1;
  }

  DSC_LOG_INFO << "Adaption strategy:" << std::endl << std::endl;
  DSC_LOG_INFO
  << "Define 'variation = (indicator_on_element - average_indicator) / (maximum_indicator - average_indicator)'"
  << std::endl;
  DSC_LOG_INFO << "Subdivide the region [average_indicator,maximum_indicator] into " << number_of_areas << " areas."
               << std::endl;
  if (number_of_areas == 1) {
    DSC_LOG_INFO << "1.: [average_indicator,maximum_indicator]. Mark elements for " << refinements_in_area[0]
                 << " refinements." << std::endl;
  } else {
    DSC_LOG_INFO << "1.: [average_indicator," << border[0]
                 << "*maximum_indicator]. If 'variance' in area: mark elements for " << refinements_in_area[0]
                 << " refinements." << std::endl;
    for (int bo = 1; bo < (number_of_areas - 1); ++bo)
      DSC_LOG_INFO << bo + 1 << ".: [" << border[bo - 1] << "*average_indicator," << border[bo]
                   << "*maximum_indicator]. If 'variance' in area: mark elements for " << refinements_in_area[bo]
                   << " refinements." << std::endl;
    DSC_LOG_INFO << number_of_areas << ".: [" << border[number_of_areas - 2]
                 << "*average_indicator,maximum_indicator]. If 'variance' in area: mark elements for "
                 << refinements_in_area[number_of_areas - 1] << " refinements." << std::endl;
  }
  DSC_LOG_INFO << "Default refinement for elements with 'variance <= 0 ': " << default_refinement << std::endl;

  if (result.estimated_error < error_tolerance_) {
    repeat = false;
    //!TODO profiler
    //    DSC_LOG_INFO << "Total HMM time = " << total_hmm_time << "s." << std::endl;
    //    DSC_LOG_INFO << std::endl << std::endl << "Total HMM time = " << total_hmm_time << "s." << std::endl <<
    // std::endl;
  } else {
    int element_number = 0;
    for (const auto& entity : discreteFunctionSpace) {
      int additional_refinement;
      if (result.local_error_indicator[element_number] <= result.average_loc_indicator) {
        additional_refinement = default_refinement;
      } else {
        const double variation = (result.local_error_indicator[element_number] - result.average_loc_indicator) /
                                 (result.maximal_loc_indicator - result.average_loc_indicator);

        if (number_of_areas == 1) {
          additional_refinement = refinements_in_area[0];
        } else {
          if (variation <= border[0]) {
            additional_refinement = refinements_in_area[0];
          }
          for (int bo = 1; bo <= (number_of_areas - 2); ++bo) {
            if ((variation > border[bo - 1]) && (variation <= border[bo])) {
              additional_refinement = refinements_in_area[bo];
            }
          }
          if (variation > border[number_of_areas - 2]) {
            additional_refinement = refinements_in_area[number_of_areas - 1];
          }
        }
      }
      discreteFunctionSpace.gridPart().grid().mark(additional_refinement, entity);
      element_number += 1;
    }
    adaptationManager.adapt();
  }
  return repeat;
}

//! the main hmm computation
void
algorithm(typename CommonTraits::GridPointerType& macro_grid_pointer, // grid pointer that belongs to the macro grid
          typename CommonTraits::GridPointerType& fine_macro_grid_pointer, // grid pointer that belongs to the fine
                                                                           // macro grid (for
                                                                           // reference computations)
          typename CommonTraits::GridPointerType& periodic_grid_pointer,   // grid pointer that belongs to the periodic
                                                                           // micro grid
          const std::string filename) {
  using namespace Dune;

  auto problem_data = Problem::getModelData();
  print_info(*problem_data, DSC_LOG_INFO);
  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for HMM-macro-problem
  typename CommonTraits::GridPartType gridPart(*macro_grid_pointer);
  // grid part for the periodic function space, required for HMM-cell-problems
  typename HMMTraits::PeriodicGridPartType periodicGridPart(*periodic_grid_pointer);
  // grid part for the global function space, required for the detailed fine-scale computation (very high resolution)
  typename CommonTraits::GridPartType gridPartFine(*fine_macro_grid_pointer);

  auto& grid = gridPart.grid();
  //! --------------------------------------------------------------------------------------

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  typename CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  // the global-problem function space for the reference computation:
  typename CommonTraits::DiscreteFunctionSpaceType finerDiscreteFunctionSpace(gridPartFine);
  // the local-problem function space (containing periodic functions):
  typename HMMTraits::PeriodicDiscreteFunctionSpaceType periodicDiscreteFunctionSpace(periodicGridPart);
  //! --------------------------------------------------------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const auto diffusion_op = Problem::getDiffusion();

  //! solution vector
  // - By reference_solution, we denote an (possibly accurate) approximation of the exact solution (used for comparison)
  // - if the elliptic problem is linear, the reference solution can be either determined with a fine scale FEM
  // computation
  //   or - if the structure is additionally periodic -- with a FEM computation for the homogenized problem
  // - if the elliptic problem is nonlinear, the reference was most likely determined with a finite element method on
  // the fine scale,
  //   where we used the Newton method to solve the non-linear system of equations
  // - in general 'reference_solution' will be an accurate approximation of the exact solution, that is why we it also
  // called reference solution
  typename CommonTraits::DiscreteFunctionType reference_solution(filename + " Reference Solution",
                                                                 finerDiscreteFunctionSpace);
  reference_solution.clear();

  if (DSC_CONFIG_GET("problem.reference_solution", false))
    load_reference(reference_solution);

  //! ************************* Assembling and solving the HMM problem ****************************
  // number of the loop cycle of the while-loop
  int loop_cycle = 1;
  const double error_tolerance_ = DSC_CONFIG_GET("hmm.error_tolerance", 1e-6);
  bool repeat = true;
  while (repeat) {
    DSC_LOG_INFO << std::endl << "########################### LOOP CYCLE " << loop_cycle
                 << " ###########################" << std::endl << std::endl << std::endl;

    //! solution vector
    // solution of the heterogeneous multiscale finite element method, where we used the Newton method to solve the
    // non-linear system of equations
    typename CommonTraits::DiscreteFunctionType hmm_solution(" HMM (+Newton) Solution", discreteFunctionSpace);
    hmm_solution.clear();

    typename CommonTraits::RestrictProlongOperatorType rp(hmm_solution);
    typename CommonTraits::AdaptationManagerType adaptationManager(grid, rp);
    const auto result = single_step(gridPart, gridPartFine, discreteFunctionSpace, periodicDiscreteFunctionSpace,
                                    *diffusion_op, hmm_solution, reference_solution, loop_cycle);
    // call of 'single_step': 'reference_solution' only required for error computation

    if (!DSC_CONFIG_GET("hmm.adaptivity", false))
      break;

    if (!adapt(result, loop_cycle, error_tolerance_, discreteFunctionSpace, adaptationManager))
      break;

    DSC_LOG_INFO << std::endl << "#########################################################################"
                 << std::endl << std::endl << std::endl;
    loop_cycle += 1;
  } // end while (repeat) of repeat loop (for the adaptive cicles)
    //! ******************** End of assembling and solving the HMM problem ***************************
}

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {
