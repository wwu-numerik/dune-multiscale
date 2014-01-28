// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_MSFEM_ALGORITHM_HH
#define DUNE_MULTISCALE_MSFEM_ALGORITHM_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <string>
#include <vector>

namespace Dune {
namespace Multiscale {

struct OutputParameters;

namespace MsFEM {

class MacroMicroGridSpecifier;
class LocalGridList;

typedef ErrorEstimator<typename CommonTraits::DiscreteFunctionType, typename CommonTraits::DiffusionType,
                            typename CommonTraits::FirstSourceType, MacroMicroGridSpecifier,
                            LocalGridList> ErrorEstimatorType;
//! \TODO docme
void adapt(CommonTraits::GridType& grid, CommonTraits::GridType& grid_coarse, const int loop_number,
           int& total_refinement_level_, int& coarse_grid_level_, int& number_of_layers_,
           const std::vector<CommonTraits::RangeVectorVector*>& locals,
           const std::vector<CommonTraits::RangeVector*>& totals,
           const CommonTraits::RangeVector& total_estimated_H1_error_);

//! \TODO docme
void solution_output(const CommonTraits::DiscreteFunction_ptr &msfem_solution,
                     const CommonTraits::DiscreteFunction_ptr &coarse_part_msfem_solution,
                     const CommonTraits::DiscreteFunction_ptr &fine_part_msfem_solution,
                     Dune::Multiscale::OutputParameters& outputparam, const int loop_number);

//! \TODO docme
void data_output(const CommonTraits::GridPartType& gridPart,
                 const CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace_coarse,
                 Dune::Multiscale::OutputParameters& outputparam, const int loop_number);

//! \TODO docme
bool error_estimation(const CommonTraits::DiscreteFunctionType& msfem_solution,
                      const CommonTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                      const CommonTraits::DiscreteFunctionType& fine_part_msfem_solution,
                      ErrorEstimatorType& estimator,
                      MacroMicroGridSpecifier& specifier, const int loop_number,
                      std::vector<CommonTraits::RangeVectorVector*>& locals,
                      std::vector<CommonTraits::RangeVector*>& totals,
                      CommonTraits::RangeVector& total_estimated_H1_error_);

//! \TODO docme
void algorithm(const std::string& macroGridName, const int loop_number, int& total_refinement_level_,
               int& coarse_grid_level_, int& number_of_layers_, std::vector<CommonTraits::RangeVectorVector*>& locals,
               std::vector<CommonTraits::RangeVector*>& totals, CommonTraits::RangeVector& total_estimated_H1_error_);

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_ALGORITHM_HH
