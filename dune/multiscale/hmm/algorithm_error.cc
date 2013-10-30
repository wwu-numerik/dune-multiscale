#include <config.h>
#include <config.h>
#include <dune/multiscale/hmm/algorithm_error.hh>

#include <string>
#include <vector>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/result.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/error_estimator.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/problems/selector.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

HMMResult estimate_error(const typename CommonTraits::GridPartType& gridPart,
                         const typename CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
                         const typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                         const typename CommonTraits::DiffusionType& diffusion_op,
                         const CellProblemNumberingManager& cp_num_manager,
                         const typename CommonTraits::DiscreteFunction_ptr& hmm_solution) {
  using namespace Dune::Stuff;
  // auxiliary grid part for the periodic function space, required for HMM-cell-problems
  typename CommonTraits::GridPartType auxiliaryGridPart(periodicDiscreteFunctionSpace.gridPart().grid());
  // and the corresponding auxiliary one:
  typename CommonTraits::DiscreteFunctionSpaceType auxiliaryDiscreteFunctionSpace(auxiliaryGridPart);
  // auxiliaryGridPart for the error estimator (the auxiliaryGridPart yields an intersection iterator, which can not be
  // done by the periodicGridPart)
  // Notation: u_H = hmm_solution
  // to load the solutions of the cell problems:
  // location of the solutions of the cell problems for the base function set:
  const std::string cell_solution_location_baseSet = "/cell_problems/_cellSolutions_baseSet";
  // location of the solutions of the cell problems for the discrete function u_H:
  const std::string cell_solution_location_discFunc = "/cell_problems/_cellSolutions_discFunc";

  auto discrete_function_reader_baseSet = DiscreteFunctionIO<typename HMMTraits::PeriodicDiscreteFunctionType>::disk(cell_solution_location_baseSet);
  auto discrete_function_reader_discFunc = DiscreteFunctionIO<typename HMMTraits::PeriodicDiscreteFunctionType>::disk(cell_solution_location_discFunc);

  const typename HMMTraits::ErrorEstimatorType error_estimator(periodicDiscreteFunctionSpace, discreteFunctionSpace,
                                                               auxiliaryDiscreteFunctionSpace, diffusion_op);

  HMMResult result(discreteFunctionSpace.grid().size(0));
  auto f = Problem::getFirstSource();
  int element_number = 0;
  for (const auto& entity : discreteFunctionSpace) {
    // corrector of u_H^(n-1) \approx u_H on the macro element T
    auto corrector_u_H_on_entity = make_df_ptr<typename HMMTraits::PeriodicDiscreteFunctionType>("Corrector of u_H",
                                                                             periodicDiscreteFunctionSpace);
    corrector_u_H_on_entity->clear();

    // in the linear case, we still need to compute the corrector of u_H:
    if (DSC_CONFIG_GET("problem.linear", true)) {
      auto corrector_of_base_func = make_df_ptr<typename HMMTraits::PeriodicDiscreteFunctionType>("Corrector of macro base function",
                                                                              periodicDiscreteFunctionSpace);
      corrector_of_base_func->clear();
      const auto local_hmm_solution = hmm_solution->localFunction(entity);
      const auto& baseSet = discreteFunctionSpace.basisFunctionSet(entity);
      const auto numMacroBaseFunctions = baseSet.size();
      std::vector<std::size_t> cell_problem_id(numMacroBaseFunctions);
      for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
        const typename CommonTraits::EntityType::EntityPointer p(entity);
        cell_problem_id[i] = cp_num_manager.get_number_of_cell_problem(p, i);
        discrete_function_reader_baseSet.read(cell_problem_id[i], corrector_of_base_func);
        *corrector_of_base_func *= local_hmm_solution[i];
        *corrector_u_H_on_entity += *corrector_of_base_func;
      }
    } else {
      // in the nonlinear case this corrector is already available
      discrete_function_reader_discFunc.read(element_number, corrector_u_H_on_entity);
    }

    // contribution of the local source error
    auto local_source_indicator = error_estimator.indicator_f(*f, entity);
    result.estimated_source_error += local_source_indicator;

    // contribution of the local approximation error
    auto local_approximation_indicator = error_estimator.indicator_app_1(entity, *hmm_solution, *corrector_u_H_on_entity);

    local_approximation_indicator += error_estimator.indicator_app_2(entity, *hmm_solution, *corrector_u_H_on_entity);
    result.estimated_approximation_error += local_approximation_indicator;

    // contribution of the local residual error
    auto local_residual_indicator = error_estimator.indicator_res_T(entity, *hmm_solution, *corrector_u_H_on_entity);
    result.estimated_residual_error_micro_jumps += local_residual_indicator;

    for (const auto& intersection : DSC::intersectionRange(gridPart, entity)) {
      if (intersection.neighbor()) // if there is a neighbor entity
      {
        // corrector of u_H^(n-1) \approx u_H on the neighbor element
        auto corrector_u_H_on_neighbor_entity = make_df_ptr<typename HMMTraits::PeriodicDiscreteFunctionType>(
            "Corrector of u_H", periodicDiscreteFunctionSpace);
        corrector_u_H_on_neighbor_entity->clear();

        auto it_outside = intersection.outside();
        const auto& entity_outside = *it_outside;

        // in the linear case, we still need to compute the corrector of u_H:
        if (DSC_CONFIG_GET("problem.linear", true)) {
          auto corrector_of_base_func_neighbor = make_df_ptr<typename HMMTraits::PeriodicDiscreteFunctionType>(
              "Corrector of macro base function", periodicDiscreteFunctionSpace);
          corrector_of_base_func_neighbor->clear();

          auto local_hmm_solution_neighbor = hmm_solution->localFunction(entity_outside);

          const auto& baseSet_neighbor = discreteFunctionSpace.basisFunctionSet(entity_outside);
          const auto numMacroBaseFunctions_neighbor = baseSet_neighbor.size();
          std::vector<std::size_t> cell_problem_id_neighbor(numMacroBaseFunctions_neighbor);
          for (unsigned int i = 0; i < numMacroBaseFunctions_neighbor; ++i) {
            cell_problem_id_neighbor[i] = cp_num_manager.get_number_of_cell_problem(it_outside, i);
            discrete_function_reader_baseSet.read(cell_problem_id_neighbor[i], corrector_of_base_func_neighbor);
            *corrector_of_base_func_neighbor *= local_hmm_solution_neighbor[i];
            *corrector_u_H_on_neighbor_entity += *corrector_of_base_func_neighbor;
          }
        } else {
          auto neighbor_element_number = cp_num_manager.get_number_of_cell_problem(it_outside);
          // in the nonlinear case this corrector is already available
          discrete_function_reader_discFunc.read(neighbor_element_number, corrector_u_H_on_neighbor_entity);
        }

        auto val = error_estimator.indicator_res_E(intersection, *hmm_solution, *corrector_u_H_on_entity,
                                                   *corrector_u_H_on_neighbor_entity);
        local_residual_indicator += val;
        result.estimated_residual_error_macro_jumps += val;
      }
    }

    result.estimated_residual_error += local_residual_indicator;

    // tfr = test function reconstruction ( = non-Petrov-Galerkin )
    typename CommonTraits::RangeType local_tfr_indicator(0);
    if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
      // use 'indicator_effective_tfr' or 'indicator_tfr_1'
      // contribution of the local tfr error:
      local_tfr_indicator = error_estimator.indicator_tfr_1(entity, *hmm_solution, *corrector_u_H_on_entity);
      result.estimated_tfr_error += local_tfr_indicator;
    }

    result.local_error_indicator[element_number] =
        local_source_indicator + local_approximation_indicator + local_residual_indicator;
    result.local_error_indicator[element_number] += local_tfr_indicator;

    //!TODO minmaxaverge
    if (result.local_error_indicator[element_number] < result.minimal_loc_indicator) {
      result.minimal_loc_indicator = result.local_error_indicator[element_number];
    }

    if (result.local_error_indicator[element_number] > result.maximal_loc_indicator) {
      result.maximal_loc_indicator = result.local_error_indicator[element_number];
    }

    result.average_loc_indicator += result.local_error_indicator[element_number];

    element_number += 1;
  } // endfor

  result.average_loc_indicator /= double(discreteFunctionSpace.grid().size(0));

  result.estimated_source_error = sqrt(result.estimated_source_error);
  result.estimated_approximation_error = sqrt(result.estimated_approximation_error);
  result.estimated_residual_error = sqrt(result.estimated_residual_error);

  result.estimated_residual_error_micro_jumps = sqrt(result.estimated_residual_error_micro_jumps);
  result.estimated_residual_error_macro_jumps = sqrt(result.estimated_residual_error_macro_jumps);

  result.estimated_error =
      result.estimated_source_error + result.estimated_approximation_error + result.estimated_residual_error;

  if (!DSC_CONFIG_GET("hmm.petrov_galerkin", true)) {
    result.estimated_tfr_error = sqrt(result.estimated_tfr_error);
    result.estimated_error += result.estimated_tfr_error;
  }

  if (DSC_CONFIG_GET("hmm.adaptivity", false)) {
    // maximum variation (up) from average
    result.max_variation = result.average_loc_indicator / result.maximal_loc_indicator;
    result.min_variation = result.average_loc_indicator / result.minimal_loc_indicator;
  }
  return result;
} //! -------- End Error Estimation --------

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {
