// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_FEM_ALGORITHM_HH
#define DUNE_FEM_ALGORITHM_HH

#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/spaces/interface.hh>

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace FEM {

template <class GridViewType>
class ProblemNineDiffusion;

template <class GridViewType>
class ProblemNineForce;

template <class GridViewType>
class ProblemNineDirichlet;

template <class GridViewType>
class ProblemNineNeumann;

template <class GridViewType>
class ProblemNineExactSolution;

template <class GridType, Stuff::Grid::ChooseLayer grid_layer = Stuff::Grid::ChooseLayer::leaf,
          GDT::ChooseSpaceBackend space_backend = GDT::ChooseSpaceBackend::pdelab,
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::ChooseBackend::istl_sparse>
class EllipticDuneGdtDiscretization;

//! the main FEM computation
void algorithm(const std::shared_ptr<const CommonTraits::GridType>& macro_grid_pointer, const std::string filename);

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_FEM_ALGORITHM_HH
