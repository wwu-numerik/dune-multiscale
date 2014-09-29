#include <config.h>

#include "error_calc.hh"

#include <assert.h>
#include <boost/filesystem/fstream.hpp>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/proxygridview.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>

#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>

using namespace Dune::Multiscale;

void solution_output(const CommonTraits::DiscreteFunctionType& msfem_solution) {
  using namespace Dune;

  Dune::Multiscale::OutputParameters outputparam;
  outputparam.set_prefix("msfem_solution_");
  msfem_solution.visualize(outputparam.fullpath(msfem_solution.name()));
}

void data_output(const CommonTraits::GridViewType& gridPart) {
  using namespace Dune;
  Dune::Multiscale::OutputParameters outputparam;

  if (Problem::getModelData()->hasExactSolution()) {
    const auto& u = Dune::Multiscale::Problem::getExactSolution();
    outputparam.set_prefix("exact_solution");
    u->visualize(gridPart, outputparam.fullpath(u->name()));
  }
}

Dune::Multiscale::ErrorCalculator::ErrorCalculator(const std::unique_ptr<LocalsolutionProxy>& msfem_solution,
                                                   CommonTraits::ConstDiscreteFunctionType* fem_solution)
  : msfem_solution_(msfem_solution)
  , fem_solution_(fem_solution)
{
  assert(fem_solution_);
}

ErrorCalculator::ErrorCalculator(const std::unique_ptr<LocalsolutionProxy> &msfem_solution)
  : msfem_solution_(msfem_solution)
  , fem_solution_(nullptr)
{
  assert(msfem_solution_);
  if (DSC_CONFIG_GET("msfem.fem_comparison", false)) {
    fem_solver_ = DSC::make_unique<Elliptic_FEM_Solver>();
    fem_solution_ = &fem_solver_->solve();
  }
}

std::map<std::string, double> Dune::Multiscale::ErrorCalculator::print(std::ostream& out) {
  using namespace Dune::GDT::Products;
  out << std::endl << "The L2 errors:" << std::endl << std::endl;

  const size_t over_integrate = 0; // <- would let the product use a higher quadrature order than needed

  typedef Stuff::Functions::Difference<Problem::ExactSolutionType, CommonTraits::ConstDiscreteFunctionType>
  DifferenceType;
  /// TODO only call assemble once
  const auto fine_grid = make_grids().second;
  CommonTraits::GridProviderType fine_grid_provider(*fine_grid);
  const auto fine_space_ptr = //fem_solution_->space();
      CommonTraits::SpaceProviderType::create(fine_grid_provider, CommonTraits::st_gdt_grid_level);
  const auto& fine_space = fem_solution_ ? fem_solution_->space() : *fine_space_ptr;
  GDT::SystemAssembler<CommonTraits::SpaceType> system_assembler(fine_space);
  const auto& grid_view = fine_space.grid_view();

  std::map<std::string, double> csv;

  CommonTraits::DiscreteFunctionType fine_msfem_solution(fine_space, "MsFEM_Solution");
  if (msfem_solution_) {
    DS::MsFEMProjection::project(*msfem_solution_, fine_msfem_solution, msfem_solution_->search());
    if (DSC_CONFIG_GET("global.vtk_output", false)) {
      DSC_LOG_INFO_0 << "Solution output for MsFEM Solution." << std::endl;
      data_output(*fine_space.grid_view());
      solution_output(fine_msfem_solution);
    }
  }

  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::getModelData()->hasExactSolution()) {
    const auto& u = *DMP::getExactSolution();

    if (msfem_solution_) {
      const DifferenceType difference(u, fine_msfem_solution);
      L2Localizable<CommonTraits::GridViewType, DifferenceType> l2_error_product(*grid_view, difference,
                                                                                 over_integrate);
      system_assembler.add(l2_error_product);
      H1SemiLocalizable<CommonTraits::GridViewType, DifferenceType> h1_semi_error_product(*grid_view, difference,
                                                                                          over_integrate);
      system_assembler.add(h1_semi_error_product);
      system_assembler.tbb_assemble();

      const auto msfem_error = std::sqrt(l2_error_product.apply2());
      out << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl;
      const auto h1_msfem_error = std::sqrt(h1_semi_error_product.apply2());
      out << "|| u_msfem - u_exact ||_H1s =  " << h1_msfem_error << std::endl << std::endl;

      csv["msfem_exact_L2"] = msfem_error;
      csv["msfem_exact_H1s"] = h1_msfem_error;
    }

    if (fem_solution_) {
      const DifferenceType difference(u, *fem_solution_);
      L2Localizable<CommonTraits::GridViewType, DifferenceType> l2_error_product(*grid_view, difference,
                                                                                 over_integrate);
      system_assembler.add(l2_error_product);
      H1SemiLocalizable<CommonTraits::GridViewType, DifferenceType> h1_semi_error_product(*grid_view, difference,
                                                                                          over_integrate);
      system_assembler.add(h1_semi_error_product);
      system_assembler.tbb_assemble();

      const auto fem_error = std::sqrt(l2_error_product.apply2());
      out << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl;
      const auto h1_fem_error = std::sqrt(h1_semi_error_product.apply2());
      out << "|| u_fem - u_exact ||_H1s =  " << h1_fem_error << std::endl << std::endl;

      csv["fem_exact_L2"] = fem_error;
      csv["fem_exact_H1s"] = h1_fem_error;
    }
  }

  if (msfem_solution_ && fem_solution_) {

    typedef Stuff::Functions::Difference<CommonTraits::ConstDiscreteFunctionType,
                                         CommonTraits::ConstDiscreteFunctionType> DiscreteDifferenceType;
    const DiscreteDifferenceType difference(fine_msfem_solution, *fem_solution_);

    L2Localizable<CommonTraits::GridViewType, DiscreteDifferenceType> l2_error_product(*grid_view, difference,
                                                                                       over_integrate);
    system_assembler.add(l2_error_product);
    L2Localizable<CommonTraits::GridViewType, CommonTraits::ConstDiscreteFunctionType> l2_msfem(
        *grid_view, fine_msfem_solution, over_integrate);
    system_assembler.add(l2_msfem);
    H1SemiLocalizable<CommonTraits::GridViewType, DiscreteDifferenceType> h1_semi_error_product(*grid_view, difference,
                                                                                                over_integrate);
    system_assembler.add(h1_semi_error_product);
    system_assembler.tbb_assemble();

    const auto approx_msfem_error = std::sqrt(l2_error_product.apply2());
    const auto no = std::sqrt(l2_msfem.apply2());
    if (std::abs(no) > 1e-12)
      out << "|| u_msfem - u_fem ||_L2 / || u_msfem ||_L2 =  " << approx_msfem_error / no << std::endl;
    else
      out << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl;

    const auto h1_approx_msfem_error = std::sqrt(h1_semi_error_product.apply2());
    out << "|| u_msfem - u_fem ||_H1s =  " << h1_approx_msfem_error << std::endl << std::endl;

    csv["msfem_fem_L2"] = approx_msfem_error;
    csv["msfem_fem_H1s"] = h1_approx_msfem_error;
  }

  std::unique_ptr<boost::filesystem::ofstream> csvfile(
      DSC::make_ofstream(std::string(DSC_CONFIG_GET("global.datadir", "data/")) + std::string("/errors.csv")));
  const std::string sep(",");
  for (const auto& key_val : csv) {
    *csvfile << key_val.first << sep;
  }
  *csvfile << std::endl;
  for (const auto& key_val : csv) {
    *csvfile << key_val.second << sep;
  }
  *csvfile << std::endl;

  return csv;
}
