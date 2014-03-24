// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_FEM_ALGORITHM_HH
#define DUNE_FEM_ALGORITHM_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/fem/fem_traits.hh>

namespace Dune {
namespace Multiscale {
namespace FEM {

template< class GridViewType >
class ProblemNineDiffusion;

template< class GridViewType >
class ProblemNineForce;

template< class GridViewType >
class ProblemNineDirichlet;

template< class GridViewType >
class ProblemNineNeumann;

template< class GridViewType >
class ProblemNineExactSolution;

//! the main FEM computation
void algorithm(const std::shared_ptr<CommonTraits::GridType>& macro_grid_pointer, const std::string filename);

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_FEM_ALGORITHM_HH
