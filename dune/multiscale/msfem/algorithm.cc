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
void adapt(CommonTraits::GridType& grid, CommonTraits::GridType& coarse_grid, const int loop_number,
           int& total_refinement_level_, int& coarse_grid_level_, int& number_of_layers_,
           const std::vector<CommonTraits::RangeVectorVector*>& locals,
           const std::vector<CommonTraits::RangeVector*>& totals,
           const CommonTraits::RangeVector& total_estimated_H1_error_) {
  using namespace Dune;
  typedef CommonTraits::GridType::LeafGridView GridView;
  typedef GridView::Codim<0>::Iterator ElementLeafIterator;
  typedef CommonTraits::GridType::Traits::LeafIndexSet GridLeafIndexSet;

  bool coarse_scale_error_dominant = false;
  // identify the dominant contribution:
  const double average_est_error = total_estimated_H1_error_[loop_number - 1] / 6.0; // 6 contributions

  const auto& total_approximation_error_ = *totals[4];
  const auto& total_fine_grid_jumps_ = *totals[5];
  if ((total_approximation_error_[loop_number - 1] >= average_est_error) ||
      (total_fine_grid_jumps_[loop_number - 1] >= average_est_error)) {
    total_refinement_level_ += 2; // 'the fine grid level'
    DSC_LOG_INFO_0 << "Fine scale error identified as being dominant. Decrease the number of global refinements by 2."
                   << std::endl;
    DSC_LOG_INFO_0 << "NEW: Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
    DSC_LOG_INFO_0 << "Note that this means: the fine grid is " << total_refinement_level_ - coarse_grid_level_
                   << " refinement levels finer than the coarse grid." << std::endl;
  }

  const auto& total_projection_error_ = *totals[1];
  const auto& total_conservative_flux_jumps_ = *totals[3];
  if ((total_projection_error_[loop_number - 1] >= average_est_error) ||
      (total_conservative_flux_jumps_[loop_number - 1] >= average_est_error)) {
    number_of_layers_ += 1;
    DSC_LOG_INFO_0
    << "Oversampling error identified as being dominant. Increase the number of layers for each subgrid by 5."
    << std::endl;
    DSC_LOG_INFO_0 << "NEW: Number of layers = " << number_of_layers_ << std::endl;
  }

  const auto& total_coarse_residual_ = *totals[0];
  const auto& total_coarse_grid_jumps_ = *totals[2];
  if ((total_coarse_residual_[loop_number - 1] >= average_est_error) ||
      (total_coarse_grid_jumps_[loop_number - 1] >= average_est_error)) {
    DSC_LOG_INFO_0 << "Coarse scale error identified as being dominant. Start adaptive coarse grid refinment."
                   << std::endl;
    coarse_scale_error_dominant = true; /* mark elementwise for 2 refinments */
  }

  const auto& loc_coarse_residual_ = *locals[0];
  const auto& loc_coarse_grid_jumps_ = *totals[1];
  std::vector<CommonTraits::RangeType> average_coarse_error_indicator(loop_number); // arithmetic average
  for (int l = 0; l < loop_number; ++l) {
    assert(loc_coarse_residual_[l].size() == loc_coarse_grid_jumps_[l].size());
    average_coarse_error_indicator[l] = 0.0;
    for (size_t i = 0; i < loc_coarse_grid_jumps_[l].size(); ++i) {
      average_coarse_error_indicator[l] += loc_coarse_residual_[l][i] + loc_coarse_grid_jumps_[l][i];
    }
    average_coarse_error_indicator[l] = average_coarse_error_indicator[l] / loc_coarse_residual_[l].size();
  }

  if (coarse_scale_error_dominant) {
    int number_of_refinements = 4;
    // allowed varianve from average ( in percent )
    double variance = 1.1; // = 110 % of the average
    for (int l = 0; l < loop_number; ++l) {
      const auto& gridLeafIndexSet = grid.leafIndexSet();
      auto gridView = grid.leafView();

      int total_number_of_entities = 0;
      int number_of_marked_entities = 0;
      for (ElementLeafIterator it = gridView.begin<0>(); it != gridView.end<0>(); ++it) {
        CommonTraits::RangeType loc_coarse_error_indicator = loc_coarse_grid_jumps_[l][gridLeafIndexSet.index(*it)] +
                                                             loc_coarse_residual_[l][gridLeafIndexSet.index(*it)];
        total_number_of_entities += 1;

        if (loc_coarse_error_indicator >= variance * average_coarse_error_indicator[l]) {
          grid.mark(number_of_refinements, *it);
          number_of_marked_entities += 1;
        }
      }

      if (l == (loop_number - 1)) {
        DSC_LOG_INFO_0 << number_of_marked_entities << " of " << total_number_of_entities
                       << " coarse grid entities marked for mesh refinement." << std::endl;
      }

      grid.preAdapt();
      grid.adapt();
      grid.postAdapt();

      GridView gridView_coarse = coarse_grid.leafView();
      const GridLeafIndexSet& gridLeafIndexSet_coarse = coarse_grid.leafIndexSet();

      for (ElementLeafIterator it = gridView_coarse.begin<0>(); it != gridView_coarse.end<0>(); ++it) {
        CommonTraits::RangeType loc_coarse_error_indicator =
            loc_coarse_grid_jumps_[l][gridLeafIndexSet_coarse.index(*it)] +
            loc_coarse_residual_[l][gridLeafIndexSet_coarse.index(*it)];

        if (loc_coarse_error_indicator >= variance * average_coarse_error_indicator[l]) {
          coarse_grid.mark(number_of_refinements, *it);
        }
      }
      coarse_grid.preAdapt();
      coarse_grid.adapt();
      coarse_grid.postAdapt();
    }
  }
}

//! \TODO docme
void solution_output(const CommonTraits::DiscreteFunctionType& msfem_solution,
                     const CommonTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                     const CommonTraits::DiscreteFunctionType& fine_part_msfem_solution,
                     const int loop_number) {
  using namespace Dune;
  const auto& grid = msfem_solution.space().gridPart().grid();
  Dune::Multiscale::OutputParameters outputparam;

  OutputTraits::IOTupleType msfem_solution_series(&msfem_solution);
  const auto msfem_fname_s = (boost::format("msfem_solution_%d_") % loop_number).str();
  outputparam.set_prefix(msfem_fname_s);
  OutputTraits::DataOutputType(grid, msfem_solution_series, outputparam).writeData(1.0 /*dummy*/, msfem_fname_s);

  OutputTraits::IOTupleType coarse_msfem_solution_series(&coarse_part_msfem_solution);
  const auto coarse_msfem_fname_s = (boost::format("coarse_part_msfem_solution_%d_") % loop_number).str();
  outputparam.set_prefix(coarse_msfem_fname_s);
  OutputTraits::DataOutputType(grid, coarse_msfem_solution_series, outputparam).writeData(1.0 /*dummy*/, coarse_msfem_fname_s);

  OutputTraits::IOTupleType fine_msfem_solution_series(&fine_part_msfem_solution);
  const auto fine_msfem_fname_s = (boost::format("fine_part_msfem_solution_%d_") % loop_number).str();
  outputparam.set_prefix(fine_msfem_fname_s);
  OutputTraits::DataOutputType(grid, fine_msfem_solution_series, outputparam).writeData(1.0 /*dummy*/, fine_msfem_fname_s);

  DSG::ElementVisualization::all(fine_part_msfem_solution.gridPart().grid(), Dune::Fem::MPIManager::helper(),
                                 outputparam.path());
}

//! \TODO docme
void data_output(const CommonTraits::GridPartType& gridPart,
                 const CommonTraits::DiscreteFunctionSpaceType& coarse_discreteFunctionSpace,
                 const int loop_number) {
  using namespace Dune;
  Dune::Multiscale::OutputParameters outputparam;

  if (Problem::getModelData()->hasExactSolution()) {
    auto u_ptr = Dune::Multiscale::Problem::getExactSolution();
    const auto& u = *u_ptr;  
    const OutputTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPart);
    OutputTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("exact_solution");
    OutputTraits::ExSolDataOutputType(gridPart.grid(), exact_solution_series, outputparam).writeData(1.0 /*dummy*/, "exact-solution");
  }
  
  CommonTraits::DiscreteFunctionType coarse_grid_visualization("Visualization of the coarse grid",
                                                               coarse_discreteFunctionSpace);
  coarse_grid_visualization.clear();
  OutputTraits::IOTupleType coarse_grid_series(&coarse_grid_visualization);
  const auto coarse_grid_fname = (boost::format("coarse_grid_visualization_%d_") % loop_number).str();
  outputparam.set_prefix(coarse_grid_fname);
  OutputTraits::DataOutputType(coarse_discreteFunctionSpace.gridPart().grid(),
                               coarse_grid_series, outputparam).writeData(1.0 /*dummy*/, coarse_grid_fname);
}


//! algorithm
void algorithm(const std::string& /*macroGridName*/, const int loop_number, int& total_refinement_level_,
               int& coarse_grid_level_, int& number_of_layers_, std::vector<CommonTraits::RangeVectorVector*>& locals,
               std::vector<CommonTraits::RangeVector*>& totals, CommonTraits::RangeVector& total_estimated_H1_error_) {
  using namespace Dune;

  auto grids = make_grids();
  CommonTraits::GridType& coarse_grid = *grids.first;
  CommonTraits::GridPartType coarse_gridPart(coarse_grid);
  CommonTraits::GridType& fine_grid = *grids.second;
  CommonTraits::GridPartType fine_gridPart(fine_grid);

  // maybe adapt if this is our second iteration
  static bool local_indicators_available_ = false;
  if (DSC_CONFIG_GET("adaptive", false) && local_indicators_available_)
    adapt(fine_grid, coarse_grid, loop_number, total_refinement_level_, coarse_grid_level_, number_of_layers_, locals,
          totals, total_estimated_H1_error_);

  CommonTraits::DiscreteFunctionSpaceType fine_discreteFunctionSpace(fine_gridPart);
  CommonTraits::DiscreteFunctionSpaceType coarse_discreteFunctionSpace(coarse_gridPart);

  auto diffusion_op_ptr = Dune::Multiscale::Problem::getDiffusion();
  const auto& diffusion_op = *diffusion_op_ptr;
  auto f_ptr = Dune::Multiscale::Problem::getFirstSource();
  const auto& f = *f_ptr;

  CommonTraits::DiscreteFunctionType msfem_solution("MsFEM Solution", fine_discreteFunctionSpace);
  msfem_solution.clear();
  CommonTraits::DiscreteFunctionType coarse_part_msfem_solution("Coarse Part MsFEM Solution", fine_discreteFunctionSpace);
  coarse_part_msfem_solution.clear();
  CommonTraits::DiscreteFunctionType fine_part_msfem_solution("Fine Part MsFEM Solution", fine_discreteFunctionSpace);
  fine_part_msfem_solution.clear();


  Elliptic_MsFEM_Solver().apply(coarse_discreteFunctionSpace, diffusion_op, f, coarse_part_msfem_solution,
                     fine_part_msfem_solution, msfem_solution);

  if (DSC_CONFIG_GET("msfem.vtkOutput", false)) {
    DSC_LOG_INFO_0 << "Solution output for MsFEM Solution." << std::endl;
    data_output(fine_gridPart, coarse_discreteFunctionSpace, loop_number);
    solution_output(msfem_solution, coarse_part_msfem_solution, fine_part_msfem_solution, loop_number);
  }

  CommonTraits::DiscreteFunctionType fem_solution("FEM Solution", fine_discreteFunctionSpace);
  if (DSC_CONFIG_GET("msfem.fem_comparison", false)) {
    fem_solution.clear();
    const Dune::Multiscale::Elliptic_FEM_Solver fem_solver(fine_discreteFunctionSpace);
    const auto l_ptr = Dune::Multiscale::Problem::getLowerOrderTerm();
    fem_solver.apply(diffusion_op, l_ptr, f, fem_solution);
    if (DSC_CONFIG_GET("msfem.vtkOutput", false)) {
      DSC_LOG_INFO_0 << "Data output for FEM Solution." << std::endl;
      Dune::Multiscale::OutputParameters outputparam;
      OutputTraits::IOTupleType fem_solution_series(&fem_solution);
      outputparam.set_prefix("fem_solution");
      OutputTraits::DataOutputType fem_dataoutput(fine_grid, fem_solution_series, outputparam);
      fem_dataoutput.writeData(1.0 /*dummy*/, "fem_solution");
    }
  }

  ErrorCalculator(&msfem_solution, &fem_solution).print(DSC_LOG_INFO_0);

  if (DSC_CONFIG_GET("adaptive", false))
    DSC_LOG_INFO_0 << "\n\n---------------------------------------------" << std::endl;
} // function algorithm

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
