#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <string>
#include <fstream>
#include <cmath>

#include <dune/common/timer.hh>

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/fem/print_info.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/global.hh>
#include <dune/stuff/grid/provider.hh>

#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/gdt/discretefunction/default.hh>


namespace Dune {
namespace Multiscale {
namespace FEM {

//! the main FEM computation
void algorithm(const std::shared_ptr< const CommonTraits::GridType >& macro_grid_pointer, const std::string /*filename*/) {

  const auto& source = Problem::getSource();
  const auto& diffusion = Problem::getDiffusion();

  Stuff::Grid::Providers::ConstDefault< CommonTraits::GridType > grid_provider(macro_grid_pointer);
  const CommonTraits::GdtSpaceType space = CommonTraits::SpaceProviderType::create(grid_provider, CommonTraits::st_gdt_grid_level);
  CommonTraits::GdtVectorType solution_vector(space.mapper().size());
  CommonTraits::DiscreteFunctionType solution(space, solution_vector);
  Elliptic_FEM_Solver(space).apply(*diffusion, *source, solution);

  solution.visualize("solution.vtk");
  ErrorCalculator(nullptr, &solution).print(DSC_LOG_INFO);
} // ... algorithm(...)

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
