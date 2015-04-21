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

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/output/entity_visualization.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/gdt/discretefunction/default.hh>

#include <cmath>
#include <iterator>
#include <memory>
#include <sstream>

#include "algorithm.hh"
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>

namespace Dune {
namespace Multiscale {

//! algorithm
std::map<std::string, double> msfem_algorithm() {
  using namespace Dune;

  DSC::ScopedTiming algo("msfem.algorithm");
  const MPIHelper::MPICommunicator& comm = Dune::MPIHelper::getCommunicator();
  DMP::ProblemContainer problem(comm, comm, DSC_CONFIG);
  auto grid = make_coarse_grid(problem);
  DSC_CONFIG.set("grids.dim", CommonTraits::world_dim);
  problem.getMutableModelData().prepare_new_evaluation(problem);

  const CommonTraits::SpaceType coarseSpace(
      CommonTraits::SpaceChooserType::PartViewType::create(*grid, CommonTraits::st_gdt_grid_level));
  std::unique_ptr<LocalsolutionProxy> msfem_solution(nullptr);

  LocalGridList localgrid_list(problem, coarseSpace);
  Elliptic_MsFEM_Solver().apply(problem, coarseSpace, msfem_solution, localgrid_list);

  if (problem.config().get("global.vtk_output", false)) {
    CommonTraits::DiscreteFunctionType coarse_grid_visualization(coarseSpace, "Visualization_of_the_coarse_grid");
    coarse_grid_visualization.visualize(OutputParameters(problem.config().get("global.datadir", "data")).fullpath(coarse_grid_visualization.name()));
  }

  if (!problem.config().get("global.skip_error", false))
    return ErrorCalculator(problem, msfem_solution).print(DSC_LOG_INFO_0);
  return decltype(ErrorCalculator(problem, msfem_solution).print(DSC_LOG_INFO_0))();
} // function algorithm

} // namespace Multiscale {
} // namespace Dune {
