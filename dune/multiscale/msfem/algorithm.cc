#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <assert.h>
#include <boost/format.hpp>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/common/output_traits.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>

#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/output/entity_visualization.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>

#include <cmath>
#include <iterator>
#include <memory>
#include <sstream>

#include "algorithm.hh"
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/common/grid_creation.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! \TODO docme
void solution_output(const CommonTraits::DiscreteFunctionType& msfem_solution,
                     const CommonTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                     const CommonTraits::DiscreteFunctionType& fine_part_msfem_solution) {
  using namespace Dune;
  const auto& grid = msfem_solution.space().gridPart().grid();
  Dune::Multiscale::OutputParameters outputparam;

  OutputTraits::IOTupleType msfem_solution_series(&msfem_solution);
  const auto msfem_fname_s = std::string("msfem_solution_");
  outputparam.set_prefix(msfem_fname_s);
  OutputTraits::DataOutputType(grid, msfem_solution_series, outputparam).writeData(1.0 /*dummy*/, msfem_fname_s);

  OutputTraits::IOTupleType coarse_msfem_solution_series(&coarse_part_msfem_solution);
  const auto coarse_msfem_fname_s = std::string("coarse_part_msfem_solution_");
  outputparam.set_prefix(coarse_msfem_fname_s);
  OutputTraits::DataOutputType(grid, coarse_msfem_solution_series, outputparam)
      .writeData(1.0 /*dummy*/, coarse_msfem_fname_s);

  OutputTraits::IOTupleType fine_msfem_solution_series(&fine_part_msfem_solution);
  const auto fine_msfem_fname_s = std::string("fine_part_msfem_solution_");
  outputparam.set_prefix(fine_msfem_fname_s);
  OutputTraits::DataOutputType(grid, fine_msfem_solution_series, outputparam)
      .writeData(1.0 /*dummy*/, fine_msfem_fname_s);

  DSG::ElementVisualization::all(fine_part_msfem_solution.gridPart().grid(), outputparam.path());
}

//! \TODO docme
void data_output(const CommonTraits::GridPartType& gridPart,
                 const CommonTraits::DiscreteFunctionSpaceType& coarse_discreteFunctionSpace) {
  using namespace Dune;
  Dune::Multiscale::OutputParameters outputparam;

  if (Problem::getModelData()->hasExactSolution()) {
    auto u_ptr = Dune::Multiscale::Problem::getExactSolution();
    const auto& u = *u_ptr;
    const OutputTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPart);
    OutputTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("exact_solution");
    OutputTraits::ExSolDataOutputType(gridPart.grid(), exact_solution_series, outputparam)
        .writeData(1.0 /*dummy*/, "exact-solution");
  }

  CommonTraits::DiscreteFunctionType coarse_grid_visualization("Visualization of the coarse grid",
                                                               coarse_discreteFunctionSpace);
  coarse_grid_visualization.clear();
  OutputTraits::IOTupleType coarse_grid_series(&coarse_grid_visualization);
  const std::string coarse_grid_name("coarse_grid_visualization_");
  outputparam.set_prefix(coarse_grid_name);
  OutputTraits::DataOutputType(coarse_discreteFunctionSpace.gridPart().grid(), coarse_grid_series, outputparam)
      .writeData(1.0 /*dummy*/, coarse_grid_name);
}

//! algorithm
void algorithm() {
  using namespace Dune;

  auto grids = make_grids();
  CommonTraits::GridType& coarse_grid = *grids.first;
  CommonTraits::GridPartType coarse_gridPart(coarse_grid);
  CommonTraits::GridType& fine_grid = *grids.second;
  CommonTraits::GridPartType fine_gridPart(fine_grid);

  CommonTraits::DiscreteFunctionSpaceType fine_discreteFunctionSpace(fine_gridPart);
  CommonTraits::DiscreteFunctionSpaceType coarse_discreteFunctionSpace(coarse_gridPart);

  auto diffusion_op_ptr = Dune::Multiscale::Problem::getDiffusion();
  const auto& diffusion_op = *diffusion_op_ptr;
  auto f_ptr = Dune::Multiscale::Problem::getFirstSource();
  const auto& f = *f_ptr;

  CommonTraits::DiscreteFunctionType msfem_solution("MsFEM Solution", fine_discreteFunctionSpace);
  msfem_solution.clear();
  CommonTraits::DiscreteFunctionType coarse_part_msfem_solution("Coarse Part MsFEM Solution",
                                                                fine_discreteFunctionSpace);
  coarse_part_msfem_solution.clear();
  CommonTraits::DiscreteFunctionType fine_part_msfem_solution("Fine Part MsFEM Solution", fine_discreteFunctionSpace);
  fine_part_msfem_solution.clear();

  Elliptic_MsFEM_Solver().apply(coarse_discreteFunctionSpace, diffusion_op, f, coarse_part_msfem_solution,
                                fine_part_msfem_solution, msfem_solution);

  if (DSC_CONFIG_GET("global.vtk_output", false)) {
    DSC_LOG_INFO_0 << "Solution output for MsFEM Solution." << std::endl;
    data_output(fine_gridPart, coarse_discreteFunctionSpace);
    solution_output(msfem_solution, coarse_part_msfem_solution, fine_part_msfem_solution);
  }

  std::unique_ptr<CommonTraits::DiscreteFunctionType> fem_solution(nullptr);
  if (DSC_CONFIG_GET("msfem.fem_comparison", false)) {
    fem_solution = DSC::make_unique<CommonTraits::DiscreteFunctionType>("FEM Solution", fine_discreteFunctionSpace);
    fem_solution->clear();
    const Dune::Multiscale::Elliptic_FEM_Solver fem_solver(fine_discreteFunctionSpace);
    const auto l_ptr = Dune::Multiscale::Problem::getLowerOrderTerm();
    fem_solver.apply(diffusion_op, l_ptr, f, *fem_solution);
    if (DSC_CONFIG_GET("global.vtk_output", false)) {
      DSC_LOG_INFO_0 << "Data output for FEM Solution." << std::endl;
      Dune::Multiscale::OutputParameters outputparam;
      OutputTraits::IOTupleType fem_solution_series(fem_solution.get());
      outputparam.set_prefix("fem_solution");
      OutputTraits::DataOutputType fem_dataoutput(fine_grid, fem_solution_series, outputparam);
      fem_dataoutput.writeData(1.0 /*dummy*/, "fem_solution");
    }
  }

  ErrorCalculator(&msfem_solution, fem_solution.get()).print(DSC_LOG_INFO_0);

  if (DSC_CONFIG_GET("adaptive", false))
    DSC_LOG_INFO_0 << "\n\n---------------------------------------------" << std::endl;
} // function algorithm

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
