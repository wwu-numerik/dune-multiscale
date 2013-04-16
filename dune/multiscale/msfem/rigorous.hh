// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_RIGOROUS_HH
#define DUNE_MULTISCALE_RIGOROUS_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! \TODO docme
void solution_output(const CommonTraits::DiscreteFunctionType& msfem_solution,
                     const CommonTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                     const CommonTraits::DiscreteFunctionType& fine_part_msfem_solution,
                     Dune::myDataOutputParameters& outputparam,
                     int& total_refinement_level_,
                     int& coarse_grid_level_);
//! \TODO docme
void data_output(const CommonTraits::GridPartType& gridPart,
                 const CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace_coarse,
                 Dune::myDataOutputParameters& outputparam );

//! \TODO docme
void algorithm(const std::string& macroGridName,
               int& total_refinement_level_,
               int& coarse_grid_level_,
               int& number_of_layers_ );

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // DUNE_MULTISCALE_RIGOROUS_HH
