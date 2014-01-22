#include <config.h>
#include <assert.h>
#include <boost/filesystem/fstream.hpp>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/functions/time.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/problems/base.hh"
#include "error_calc.hh"

Dune::Multiscale::ErrorCalculator::ErrorCalculator(const CommonTraits::DiscreteFunctionType* const msfem_solution,
                                                   const CommonTraits::DiscreteFunctionType* const fem_solution)
  : msfem_solution_(msfem_solution)
  , fem_solution_(fem_solution) {}

void Dune::Multiscale::ErrorCalculator::print(std::ostream& out) {
  assert(msfem_solution_ || fem_solution_);
  out << std::endl << "The L2 errors:" << std::endl << std::endl;

  auto gridPart = fem_solution_ ? fem_solution_->gridPart() : msfem_solution_->gridPart();
  Dune::Fem::H1Norm<CommonTraits::GridPartType> h1norm(gridPart);
  Dune::Fem::L2Error<typename CommonTraits::DiscreteFunctionType> l2error;

  std::map<std::string, double> csv;

  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::getModelData()->hasExactSolution()) {
    auto u_ptr = Dune::Multiscale::Problem::getExactSolution();
    const auto& u = *u_ptr;
    const int experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit = 6;
    const Dune::Fem::GridFunctionAdapter<CommonTraits::ExactSolutionType, CommonTraits::GridPartType> u_disc(
        "", u, gridPart, experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit);

    if (msfem_solution_) {
      auto msfem_error = l2error.norm(timefunctionAdapted(u), *msfem_solution_);
      out << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl;

      auto h1_msfem_error = h1norm.distance(u_disc, *msfem_solution_);
      out << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;

      csv["msfem_exact_L2"] = msfem_error;
      csv["msfem_exact_H1"] = h1_msfem_error;
    }

    if (fem_solution_) {
      auto fem_error = l2error.norm(timefunctionAdapted(u), *fem_solution_);
      out << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl;

      auto h1_fem_error = h1norm.distance(u_disc, *fem_solution_);
      out << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;

      csv["fem_exact_L2"] = fem_error;
      csv["fem_exact_H1"] = h1_fem_error;
    }
  }
  if (msfem_solution_ && fem_solution_) {
    auto approx_msfem_error = l2error.norm2<2 * CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2>(
        *fem_solution_, *msfem_solution_);
    out << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl;

    auto h1_approx_msfem_error = h1norm.distance(*fem_solution_, *msfem_solution_);
    out << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;

    csv["msfem_fem_L2"] = approx_msfem_error;
    csv["msfem_fem_H1"] = h1_approx_msfem_error;
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
