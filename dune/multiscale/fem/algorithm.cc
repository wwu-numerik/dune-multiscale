#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/common/newton_rhs.hh>
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

#include <dune/fem/misc/l2norm.hh>

#include <string>
#include <fstream>

#include "fem_traits.hh"

namespace Dune {
namespace Multiscale {
namespace FEM {


//! the main FEM computation
void algorithm(const std::shared_ptr<CommonTraits::GridType> &macro_grid_pointer, const std::string filename) {
  using namespace Dune;
  const auto problem_data = Problem::getModelData();
  print_info(*problem_data, DSC_LOG_INFO);
  typename CommonTraits::GridPartType gridPart(*macro_grid_pointer);
  typename CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const auto diffusion_op = Problem::getDiffusion();
  auto discrete_solution = make_df_ptr<typename CommonTraits::DiscreteFunctionType>(filename + " FEM(-Newton) Solution",
                                                                discreteFunctionSpace);
  discrete_solution->clear();

  const Dune::Multiscale::Elliptic_FEM_Solver fem_solver(discreteFunctionSpace);
  const auto l_ptr = Dune::Multiscale::Problem::getLowerOrderTerm();
  const auto f_ptr = Dune::Multiscale::Problem::getFirstSource();
  fem_solver.apply(*diffusion_op, l_ptr, *f_ptr, *discrete_solution, true /*use_smp*/);

  // write FEM solution to a file and produce a VTK output
  write_discrete_function(discrete_solution, "fem");
  ErrorCalculator(nullptr, discrete_solution.get()).print(DSC_LOG_INFO_0);
  DiscreteFunctionIO<CommonTraits::DiscreteFunctionType>::clear();
}


} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
