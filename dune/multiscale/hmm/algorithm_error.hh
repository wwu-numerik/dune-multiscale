#ifndef ALGORITHM_ERROR_HH
#define ALGORITHM_ERROR_HH

#include <string>
#include <vector>

template <class HMMTraits>
struct HMMResult {
    typedef HMMTraits HMM;
    typename HMM::RangeType estimated_source_error;
    typename HMM::RangeType estimated_approximation_error;
    typename HMM::RangeType estimated_residual_error;
    typename HMM::RangeType estimated_residual_error_micro_jumps;
    typename HMM::RangeType estimated_residual_error_macro_jumps;
    typename HMM::RangeType estimated_tfr_error;
    std::vector<typename HMM::RangeType> local_error_indicator;
    typename HMM::RangeType minimal_loc_indicator;
    typename HMM::RangeType maximal_loc_indicator;
    typename HMM::RangeType average_loc_indicator;
    typename HMM::RangeType estimated_error;
    double max_variation;
    double min_variation;
    HMMResult(std::size_t codim0_count)
        : estimated_source_error(0.0)
        , estimated_approximation_error(0.0)
        , estimated_residual_error(0.0)
        , estimated_residual_error_micro_jumps(0.0)
        , estimated_residual_error_macro_jumps(0.0)
        , estimated_tfr_error(0.0)
        , local_error_indicator(codim0_count, 0.0)
        , minimal_loc_indicator(10000.0)
        , maximal_loc_indicator(0.0)
        , average_loc_indicator(0.0)
        , estimated_error(0.0)
        , max_variation(0.0)
        , min_variation(0.0)
    {}
};


//! ---------- Error Estimation ----------
template < class ProblemDataType, class HMMTraits >
HMMResult<HMMTraits>  estimate_error(
        const typename HMMTraits::GridPartType& gridPart,
        const typename HMMTraits::GridPartType& /*gridPartFine*/,
        const typename HMMTraits::DiscreteFunctionSpaceType& discreteFunctionSpace,
        const typename HMMTraits::PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
        const typename HMMTraits::DiffusionType& diffusion_op,
        const Dune::RightHandSideAssembler< typename HMMTraits::DiscreteFunctionType >& /*rhsassembler*/,
        const std::ofstream& /*data_file*/,
        const std::string filename,
        const typename HMMTraits::CellProblemNumberingManagerType& cp_num_manager,
        const typename HMMTraits::DiscreteFunctionType& hmm_solution
         )
{
typedef HMMTraits HMM;
// auxiliary grid part for the periodic function space, required for HMM-cell-problems
typename HMM::GridPartType auxiliaryGridPart(periodicDiscreteFunctionSpace.gridPart().grid());
// and the corresponding auxiliary one:
  typename HMM::DiscreteFunctionSpaceType auxiliaryDiscreteFunctionSpace(auxiliaryGridPart);
    // auxiliaryGridPart for the error estimator (the auxiliaryGridPart yields an intersection iterator, which can not be
    // done by the periodicGridPart)
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

  const typename HMM::ErrorEstimatorType error_estimator(periodicDiscreteFunctionSpace,
                                     discreteFunctionSpace,
                                     auxiliaryDiscreteFunctionSpace,
                                     diffusion_op);

  HMMResult<HMM> result(discreteFunctionSpace.grid().size(0));
  const typename HMM::FirstSourceType f;   // standard source f
  int element_number = 0;
  typedef typename HMM::DiscreteFunctionSpaceType::IteratorType IteratorType;
  IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    const typename HMM::EntityType& entity = *it;
    // corrector of u_H^(n-1) \approx u_H on the macro element T
    typename HMM::PeriodicDiscreteFunctionType corrector_u_H_on_entity("Corrector of u_H", periodicDiscreteFunctionSpace);
    corrector_u_H_on_entity.clear();

    // in the linear case, we still need to compute the corrector of u_H:
    if (DSC_CONFIG.get("problem.linear", true)) {
      typename HMM::PeriodicDiscreteFunctionType corrector_of_base_func("Corrector of macro base function",
                                                          periodicDiscreteFunctionSpace);
      corrector_of_base_func.clear();
      typename HMM::DiscreteFunctionType::LocalFunctionType local_hmm_solution = hmm_solution.localFunction(entity);
      const typename HMM::BaseFunctionSetType& baseSet = discreteFunctionSpace.baseFunctionSet(entity);
      const unsigned int numMacroBaseFunctions = baseSet.size();
      int cell_problem_id[numMacroBaseFunctions];
      for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
      {
        const typename HMM::EntityType::EntityPointer p(*it);
        cell_problem_id[i] = cp_num_manager.get_number_of_cell_problem(p, i);
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
    typename HMM::RangeType local_source_indicator = error_estimator.indicator_f(f, entity);
    result.estimated_source_error += local_source_indicator;

    // contribution of the local approximation error
    typename HMM::RangeType local_approximation_indicator = error_estimator.indicator_app_1(entity,
                                                                              hmm_solution,
                                                                              corrector_u_H_on_entity);

    local_approximation_indicator += error_estimator.indicator_app_2(entity, hmm_solution, corrector_u_H_on_entity);
    result.estimated_approximation_error += local_approximation_indicator;

    // contribution of the local residual error
    typename HMM::RangeType local_residual_indicator = error_estimator.indicator_res_T(entity, hmm_solution, corrector_u_H_on_entity);
    result.estimated_residual_error_micro_jumps += local_residual_indicator;

    typedef typename HMM::GridPartType::IntersectionIteratorType IntersectionIteratorType;
    IntersectionIteratorType endnit = gridPart.iend(entity);
    for (IntersectionIteratorType nit = gridPart.ibegin(entity); nit != endnit; ++nit)
    {
      if ( nit->neighbor() )           // if there is a neighbor entity
      {
        // corrector of u_H^(n-1) \approx u_H on the neighbor element
        typename HMM::PeriodicDiscreteFunctionType corrector_u_H_on_neighbor_entity("Corrector of u_H", periodicDiscreteFunctionSpace);
        corrector_u_H_on_neighbor_entity.clear();

        typename HMM::EntityPointerType it_outside = nit->outside();
        const typename HMM::EntityType& entity_outside = *it_outside;

        // in the linear case, we still need to compute the corrector of u_H:
        if (DSC_CONFIG.get("problem.linear", true)) {
          typename HMM::PeriodicDiscreteFunctionType corrector_of_base_func_neighbor("Corrector of macro base function",
                                                                       periodicDiscreteFunctionSpace);
          corrector_of_base_func_neighbor.clear();

          typename HMM::DiscreteFunctionType::LocalFunctionType local_hmm_solution_neighbor = hmm_solution.localFunction(entity_outside);

          const typename HMM::BaseFunctionSetType& baseSet_neighbor = discreteFunctionSpace.baseFunctionSet(entity_outside);
          const unsigned int numMacroBaseFunctions_neighbor = baseSet_neighbor.size();
          std::vector<int> cell_problem_id_neighbor(numMacroBaseFunctions_neighbor);
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

        typename HMM::RangeType val = error_estimator.indicator_res_E(*nit,
                                                        hmm_solution,
                                                        corrector_u_H_on_entity,
                                                        corrector_u_H_on_neighbor_entity);
        local_residual_indicator += val;
        result.estimated_residual_error_macro_jumps += val;
      }
    }

    result.estimated_residual_error += local_residual_indicator;

    #ifdef TFR
    // use 'indicator_effective_tfr' or 'indicator_tfr_1'
    // contribution of the local tfr error:
    typename HMM::RangeType local_tfr_indicator = error_estimator.indicator_tfr_1(entity, hmm_solution, corrector_u_H_on_entity);
    estimated_tfr_error += local_tfr_indicator;
    #endif // ifdef TFR

    result.local_error_indicator[element_number] = local_source_indicator
                                            + local_approximation_indicator
                                            + local_residual_indicator;

    #ifdef TFR
    local_error_indicator[element_number] += local_tfr_indicator;
    #endif

    if (result.local_error_indicator[element_number] < result.minimal_loc_indicator)
    {
      result.minimal_loc_indicator = result.local_error_indicator[element_number];
    }

    if (result.local_error_indicator[element_number] > result.maximal_loc_indicator)
    {
      result.maximal_loc_indicator = result.local_error_indicator[element_number];
    }

    result.average_loc_indicator += result.local_error_indicator[element_number];

    element_number += 1;
  }       // endfor

  result.average_loc_indicator /= double(discreteFunctionSpace.grid().size(0));

  result.estimated_source_error = sqrt(result.estimated_source_error);
  result.estimated_approximation_error = sqrt(result.estimated_approximation_error);
  result.estimated_residual_error = sqrt(result.estimated_residual_error);

  result.estimated_residual_error_micro_jumps = sqrt(result.estimated_residual_error_micro_jumps);
  result.estimated_residual_error_macro_jumps = sqrt(result.estimated_residual_error_macro_jumps);

  result.estimated_error = result.estimated_source_error +
          result.estimated_approximation_error + result.estimated_residual_error;

  #ifdef TFR
  result.estimated_tfr_error = sqrt(result.estimated_tfr_error);
  result.estimated_error += result.estimated_tfr_error;
  #endif // ifdef TFR

  #ifdef ADAPTIVE
  // maximum variation (up) from average
  result.max_variation = result.average_loc_indicator / result.maximal_loc_indicator;
  result.min_variation = result.average_loc_indicator / result.minimal_loc_indicator;
  #endif // ifdef ADAPTIVE
  return result;
}// ! -------- End Error Estimation --------


#endif // ALGORITHM_ERROR_HH
