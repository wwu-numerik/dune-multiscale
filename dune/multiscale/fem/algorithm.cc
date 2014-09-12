#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <string>

#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/msfem/fem_solver.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/grid/provider.hh>

namespace Dune {
namespace Multiscale {

//! the main FEM computation
void cgfem_algorithm() {
  Elliptic_FEM_Solver solver;
  auto& solution = solver.solve();

  if (DSC_CONFIG_GET("global.vtk_output", false)) {
    DSC_LOG_INFO_0 << "Solution output for FEM Solution." << std::endl;
    solution.visualize("solution.vtk");
  }
  ErrorCalculator(nullptr, &solution).print(DSC_LOG_INFO_0);
} // ... algorithm(...)

} // namespace Multiscale {
} // namespace Dune {
