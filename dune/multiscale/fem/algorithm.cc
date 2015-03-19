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
  const MPIHelper::MPICommunicator& comm = Dune::MPIHelper::getCommunicator();
  Problem::getMutableModelData().problem_init(comm, comm);
  Problem::getMutableModelData().prepare_new_evaluation();

  Elliptic_FEM_Solver solver;
  auto& solution = solver.solve();

  if (DSC_CONFIG_GET("global.vtk_output", false)) {
    DSC_LOG_INFO_0 << "Solution output for FEM Solution." << std::endl;
    Dune::Multiscale::OutputParameters outputparam;
    outputparam.set_prefix("fine-cg-fem_solution_");
    solution.visualize(outputparam.fullpath(solution.name()));
  }
  if (!DSC_CONFIG_GET("global.skip_error", false))
    ErrorCalculator(nullptr, &solution).print(DSC_LOG_INFO_0);
} // ... algorithm(...)

} // namespace Multiscale {
} // namespace Dune {
