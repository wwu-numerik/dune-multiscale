#include <config.h>
#include <assert.h>
#include <boost/filesystem/fstream.hpp>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/functions/time.hh>
#include <dune/stuff/functions/norm.hh>
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
                                                   const CommonTraits::PdelabVectorType * const fem_solution)
  : msfem_solution_(msfem_solution)
  , fem_solution_(fem_solution) {}


/** Compute and print errors between exact, msfem and fem solution
 *
 *  This computes the errors between the exact solution, the msfem solution
 *  and the fem solution (if they exist respectively). The results are printed
 *  to the given ostream and returned in a std::map.
 *
 *  \param[in] out The ostream for error printing
 *  \return Returns the computed errors as std::map with the following entries:
 *          - "msfem_exact_L2" The L2 error between the msfem solution and the
 *            exact solution
 *          - "msfem_exact_H1" The H1 error between the msfem solution and the
 *            exact solution
 *          - "fem_exact_L2" The L2 error between the fem solution and the
 *            exact solution
 *          - "fem_exact_H1" The H1 error between the fem solution and the
 *            exact solution
 *          - "msfem_fem_L2" The L2 error between the msfem solution and the
 *            fem solution
 *          Here, each entry will only be present if the corresponding error was 
 *          computed.
 */
std::map<std::string, double> Dune::Multiscale::ErrorCalculator::print(std::ostream& out) {
  assert(msfem_solution_ || fem_solution_);
  out << std::endl << "The L2 errors:" << std::endl << std::endl;

  std::map<std::string, double> csv;

  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::getModelData()->hasExactSolution()) {
    const auto& u = *Dune::Multiscale::Problem::getExactSolution();

    typedef Dune::Fem::LagrangeDiscreteFunctionSpace<CommonTraits::FunctionSpaceType,
        CommonTraits::GridPartType, CommonTraits::exact_solution_space_order>
    ExactSpaceType;
    typedef BackendChooser<ExactSpaceType>::DiscreteFunctionType ExactDiscreteFunctionType;

    #if 0 // "quality assurance" for discrete exact solution
    {
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
      auto gridPart = msfem_solution_->gridPart();
      ExactSpaceType u_disc_space(gridPart);
      ExactDiscreteFunctionType u_disc("", u_disc_space);
      Fem::LagrangeInterpolation<CommonTraits::ExactSolutionType, ExactDiscreteFunctionType>::apply(u, u_disc);

      const auto msfem_error = DS::l2distance(u_disc, *msfem_solution_);
      out << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl;

      const auto h1_msfem_error = DS::h1distance(u_disc, *msfem_solution_);
      out << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;

      csv["msfem_exact_L2"] = msfem_error;
      csv["msfem_exact_H1"] = h1_msfem_error;
    }

    if (fem_solution_) {
      const auto fem_error = DS::l2distance(u, *fem_solution_);
      out << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl;

      const auto h1_fem_error = DS::h1distance(u, *fem_solution_);
      out << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;

      csv["fem_exact_L2"] = fem_error;
      csv["fem_exact_H1"] = h1_fem_error;
    }
  }
  if (msfem_solution_ && fem_solution_) {
    const auto approx_msfem_error = DS::l2distance(*msfem_solution_,*fem_solution_);
    const auto no = DS::l2norm(*msfem_solution_);
    if (std::abs(no)>1e-12)
      out << "|| u_msfem - u_fem ||_L2 / || u_msfem ||_L2 =  " << approx_msfem_error/no << std::endl;
    else
      out << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl;

//    const auto h1_approx_msfem_error = DS::h1distance(*msfem_solution_,*fem_solution_);
//    out << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;

    csv["msfem_fem_L2"] = approx_msfem_error;
//    csv["msfem_fem_H1"] = h1_approx_msfem_error;
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
  
  return csv;
}
