#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <assert.h>
#include <boost/format.hpp>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/multiscale/common/error_calc.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/fem/print_info.hh>

#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/output/entity_visualization.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/gdt/discretefunction/default.hh>

#include <cmath>
#include <iterator>
#include <memory>
#include <sstream>

#include "algorithm.hh"
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! \TODO docme
void solution_output(const std::unique_ptr<LocalsolutionProxy>& msfem_solution) {
  using namespace Dune;

  Dune::Multiscale::OutputParameters outputparam;
  outputparam.set_prefix("msfem_solution_");
  msfem_solution->visualize(outputparam.fullpath(msfem_solution->name()));
}

//! \TODO docme
void data_output(const CommonTraits::GridViewType& gridPart,
                 const CommonTraits::DiscreteFunctionSpaceType& coarseSpace) {
  using namespace Dune;
  Dune::Multiscale::OutputParameters outputparam;

  if (Problem::getModelData()->hasExactSolution()) {
    const auto& u = Dune::Multiscale::Problem::getExactSolution();
    outputparam.set_prefix("exact_solution");
    u->visualize(gridPart, outputparam.fullpath(u->name()));
  }

  CommonTraits::DiscreteFunctionType coarse_grid_visualization(coarseSpace, "Visualization_of_the_coarse_grid");
  coarse_grid_visualization.visualize(outputparam.fullpath(coarse_grid_visualization.name()));
}

//! algorithm
std::map<std::string, double> algorithm() {
  using namespace Dune;

  auto grids = make_grids();
  CommonTraits::GridProviderType coarse_grid_provider(*grids.first);
  CommonTraits::GridProviderType fine_grid_provider(*grids.second);
  const CommonTraits::GdtSpaceType coarseSpace =
      CommonTraits::SpaceProviderType::create(coarse_grid_provider, CommonTraits::st_gdt_grid_level);
  const CommonTraits::GdtSpaceType fineSpace =
      CommonTraits::SpaceProviderType::create(fine_grid_provider, CommonTraits::st_gdt_grid_level);

  std::unique_ptr<LocalsolutionProxy> msfem_solution(nullptr);

  Elliptic_MsFEM_Solver().apply(coarseSpace, msfem_solution);

  if (DSC_CONFIG_GET("global.vtk_output", false)) {
    DSC_LOG_INFO_0 << "Solution output for MsFEM Solution." << std::endl;
    data_output(*fineSpace.grid_view(), coarseSpace);
    solution_output(msfem_solution);
  }

  std::unique_ptr<CommonTraits::DiscreteFunctionType> fem_solution(nullptr);

  if (DSC_CONFIG_GET("msfem.fem_comparison", false)) {
    fem_solution = DSC::make_unique<CommonTraits::DiscreteFunctionType>(fineSpace, "fem_solution");
    const Dune::Multiscale::Elliptic_FEM_Solver fem_solver(fineSpace);
    fem_solver.apply(*fem_solution);
    if (DSC_CONFIG_GET("global.vtk_output", false)) {
      Dune::Multiscale::FEM::write_discrete_function(*fem_solution, "fem_solution");
    }
  }

  return ErrorCalculator(msfem_solution, fem_solution.get()).print(DSC_LOG_INFO_0);

} // function algorithm

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
