#include <config.h>
#include <assert.h>
#include <boost/filesystem/fstream.hpp>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/functions/norm.hh>
#include <dune/stuff/functions/femadapter.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/problems/base.hh"
#include "error_calc.hh"

Dune::Multiscale::ErrorCalculator::ErrorCalculator(const CommonTraits::DiscreteFunctionType* const msfem_solution,
                                                   const CommonTraits::ConstDiscreteFunctionType * const fem_solution)
  : msfem_solution_(msfem_solution)
  , fem_solution_(fem_solution) {}

void Dune::Multiscale::ErrorCalculator::print(std::ostream& out) {
  assert(msfem_solution_ || fem_solution_);
  out << std::endl << "The L2 errors:" << std::endl << std::endl;

  const size_t over_integrate = 0; // <- would let the product use a higher quadrature oder than needed

  typedef Stuff::Functions::Difference< CommonTraits::ExactSolutionType, CommonTraits::ConstDiscreteFunctionType > DifferenceType;
  /// TODO this should actually select the space from either non-null solution, once msfem is gdt too
  /// also only call assemble once
  GDT::SystemAssembler< CommonTraits::GdtSpaceType > system_assembler(fem_solution_->space());
  const auto& grid_view = fem_solution_->space().grid_view();

  std::map<std::string, double> csv;

  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::getModelData()->hasExactSolution()) {
    const auto& u = *Dune::Multiscale::Problem::getExactSolution();


    #if 0 // "quality assurance" for discrete exact solution
    {
      typedef Dune::Fem::LagrangeDiscreteFunctionSpace<CommonTraits::FunctionSpaceType,
          CommonTraits::GridPartType, CommonTraits::exact_solution_space_order>
      ExactSpaceType;
      typedef BackendChooser<ExactSpaceType>::DiscreteFunctionType ExactDiscreteFunctionType;
      const int experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit = 6;
      const Dune::Fem::GridFunctionAdapter<CommonTraits::ExactSolutionType, CommonTraits::GridPartType> u_disc_adapter(
          "", u, gridPart, experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit);

      const auto ana_error = DS::l2distance(u_disc_adapter, u_disc);
      out << "|| u_exact_ana - u_exact ||_L2 =  " << ana_error << std::endl;
      const auto h1_ana_error = DS::h1distance(u_disc_adapter, u_disc);
      out << "|| u_exact_ana - u_exact ||_H2 =  " << h1_ana_error << "\n\n";

      csv["msfem_exact_L2"] = ana_error;
      csv["msfem_exact_H1"] = h1_ana_error;
    }
    #endif

    if (msfem_solution_) {
      const DifferenceType difference(u, *msfem_solution_);
      GDT::Products::L2Localizable< CommonTraits::GridViewType, DifferenceType >
          l2_error_product(*grid_view, difference, over_integrate);
      system_assembler.add(l2_error_product);
      GDT::Products::H1SemiLocalizable< CommonTraits::GridViewType, DifferenceType >
          h1_semi_error_product(*grid_view, difference, over_integrate);
      system_assembler.add(h1_semi_error_product);
      system_assembler.assemble();

      const auto msfem_error = std::sqrt(l2_error_product.apply2());
      out << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl;
      const auto h1_msfem_error = std::sqrt(h1_semi_error_product.apply2());
      out << "|| u_msfem - u_exact ||_H1s =  " << h1_msfem_error << std::endl << std::endl;

      csv["msfem_exact_L2"] = msfem_error;
      csv["msfem_exact_H1s"] = h1_msfem_error;
    }

    if (fem_solution_) {
      const DifferenceType difference(u, *fem_solution_);
      GDT::Products::L2Localizable< CommonTraits::GridViewType, DifferenceType >
          l2_error_product(*grid_view, difference, over_integrate);
      system_assembler.add(l2_error_product);
      GDT::Products::H1SemiLocalizable< CommonTraits::GridViewType, DifferenceType >
          h1_semi_error_product(*grid_view, difference, over_integrate);
      system_assembler.add(h1_semi_error_product);
      system_assembler.assemble();

      const auto fem_error = std::sqrt(l2_error_product.apply2());
      out << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl;
      const auto h1_fem_error = std::sqrt(h1_semi_error_product.apply2());
      out << "|| u_fem - u_exact ||_H1s =  " << h1_fem_error << std::endl << std::endl;

      csv["fem_exact_L2"] = fem_error;
      csv["fem_exact_H1s"] = h1_fem_error;
    }
  }
  if (msfem_solution_ && fem_solution_) {

    typedef Stuff::Functions::Difference< CommonTraits::ConstDiscreteFunctionType,
        CommonTraits::ConstDiscreteFunctionType > DiscreteDifferenceType;
    const DiscreteDifferenceType difference(*msfem_solution_, *fem_solution_);

    GDT::Products::L2Localizable< CommonTraits::GridViewType, DiscreteDifferenceType >
        l2_error_product(*grid_view, difference, over_integrate);
    system_assembler.add(l2_error_product);
    GDT::Products::L2Localizable< CommonTraits::GridViewType, CommonTraits::ConstDiscreteFunctionType >
        l2_msfem(*grid_view, *msfem_solution_, over_integrate);
    system_assembler.add(l2_msfem);
    GDT::Products::H1SemiLocalizable< CommonTraits::GridViewType, DiscreteDifferenceType >
        h1_semi_error_product(*grid_view, difference, over_integrate);
    system_assembler.add(h1_semi_error_product);
    system_assembler.assemble();

    const auto approx_msfem_error = std::sqrt(l2_error_product.apply2());
    const auto no = std::sqrt(l2_msfem.apply2());
    if (std::abs(no)>1e-12)
      out << "|| u_msfem - u_fem ||_L2 / || u_msfem ||_L2 =  " << approx_msfem_error/no << std::endl;
    else
      out << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl;

    const auto h1_approx_msfem_error = std::sqrt(h1_semi_error_product.apply2());
    out << "|| u_msfem - u_fem ||_H1s =  " << h1_approx_msfem_error << std::endl << std::endl;

    csv["msfem_fem_L2"] = approx_msfem_error;
    csv["msfem_fem_H1s"] = h1_approx_msfem_error;
  }

  std::unique_ptr<boost::filesystem::ofstream> csvfile(
      DSC::make_ofstream(DSC_CONFIG_GET("global.datadir", "data") + "/errors.csv"));
  const std::string sep(",");
  for (const auto& key_val : csv) {
    *csvfile << key_val.first << sep;
  }
  *csvfile << std::endl;
  for (const auto& key_val : csv) {
    *csvfile << key_val.second << sep;
  }
  *csvfile << std::endl;
}
