#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/output_traits.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/fem/print_info.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

#include <string>
#include <fstream>

#include "fem_traits.hh"

namespace Dune {
namespace Multiscale {
namespace FEM {

//! the main FEM computation
void algorithm(const std::shared_ptr<CommonTraits::GridType>& macro_grid_pointer, const std::string /*filename*/) {
  using namespace Dune;
  print_info(*Problem::getModelData(), DSC_LOG_INFO);

  const auto& grid_view = macro_grid_pointer->leafGridView();
  CommonTraits::FEMapType fe_map(grid_view);
  CommonTraits::GridFunctionSpaceType space(grid_view, fe_map);

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const auto& diffusion_op = Problem::getDiffusion();
  CommonTraits::PdelabVectorType solution(space, 0.0);

  const Dune::Multiscale::Elliptic_FEM_Solver fem_solver(space);
  const auto& f_ptr = Dune::Multiscale::Problem::getSource();
  fem_solver.apply(*diffusion_op, *f_ptr, solution);

  // write FEM solution to a file and produce a VTK output
  if (DSC_CONFIG_GET("global.vtk_output", false)) {
    DSC_LOG_INFO_0 << "Solution output for FEM Solution." << std::endl;
    write_discrete_function(solution, "fem");
  }
  ErrorCalculator(nullptr, &solution).print(DSC_LOG_INFO_0);

  DiscreteFunctionIO<CommonTraits::DiscreteFunctionType>::clear();
}

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
