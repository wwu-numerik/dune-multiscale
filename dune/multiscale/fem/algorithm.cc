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

#include <dune/xt/common/logging.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/stuff/grid/provider.hh>

namespace Dune {
namespace Multiscale {

//! the main FEM computation
std::map<std::string, double> cgfem_algorithm()
{
  const auto& comm = Dune::MPIHelper::getCommunicator();
  DXTC_CONFIG.set("grids.dim", CommonTraits::world_dim, true);
  DMP::ProblemContainer problem(comm, comm, DXTC_CONFIG);
  problem.getMutableModelData().prepare_new_evaluation(problem);

  Elliptic_FEM_Solver solver(problem);
  auto& solution = solver.solve();

  if (problem.config().get("global.vtk_output", false)) {
    MS_LOG_INFO_0 << "Solution output for FEM Solution." << std::endl;
    Dune::Multiscale::OutputParameters outputparam(problem.config());
    problem.getDiffusion().visualize(solution.space().grid_view(),
                                     OutputParameters(problem.config()).fullpath("msfem_diffusion"));
    outputparam.set_prefix("fine-cg-fem_solution_");
    solution.visualize(outputparam.fullpath(solution.name()));
  }
  if (!problem.config().get("msfem.skip_error", false))
    return ErrorCalculator(problem, nullptr, &solution).print(MS_LOG_INFO_0);
  return std::map<std::string, double>();
} // ... algorithm(...)

} // namespace Multiscale {
} // namespace Dune {
