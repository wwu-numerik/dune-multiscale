#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <string>

#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/grid/provider.hh>

namespace Dune {
namespace Multiscale {

//! the main FEM computation
void cgfem_algorithm() {
  const auto& comm = Dune::MPIHelper::getCommunicator();
  DSC_CONFIG.set("grids.dim", CommonTraits::world_dim);
  DMP::ProblemContainer problem(comm, comm, DSC_CONFIG);
  problem.getMutableModelData().prepare_new_evaluation(problem);

  Elliptic_FEM_Solver solver(problem);
  auto& solution = solver.solve();

  if (DSC_CONFIG_GET("global.vtk_output", false)) {
    DSC_LOG_INFO_0 << "Solution output for FEM Solution." << std::endl;
    Dune::Multiscale::OutputParameters outputparam;
    outputparam.set_prefix("fine-cg-fem_solution_");
    solution.visualize(outputparam.fullpath(solution.name()));
  }
  if (!DSC_CONFIG_GET("global.skip_error", false))
    ErrorCalculator(problem, nullptr, &solution).print(DSC_LOG_INFO_0);
} // ... algorithm(...)

} // namespace Multiscale {
} // namespace Dune {
