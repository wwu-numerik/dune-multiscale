#include <config.h>

#include "print_info.hh"

#include <boost/format.hpp>
#include <boost/optional/optional.hpp>

#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/gdt/discretefunction/default.hh>

#include <memory>
#include <sstream>


namespace Dune {
namespace Multiscale {
namespace FEM {

void write_discrete_function(typename CommonTraits::DiscreteFunction_ptr& discrete_solution, const std::string prefix) {
  // writing paraview data output
  // general output parameters
  Dune::Multiscale::OutputParameters outputparam;
  outputparam.set_prefix((boost::format("%s_solution") % prefix).str());
  discrete_solution->visualize(outputparam.prefix());

  if (Problem::getModelData()->hasExactSolution()) {
    const auto& u = *Dune::Multiscale::Problem::getExactSolution();
    outputparam.set_prefix("exact_solution");
    u.visualize(discrete_solution->space(), outputparam.prefix());
  }
}

void print_info(const CommonTraits::ModelProblemDataType& info, std::ostream& out) {
  // epsilon is specified in the parameter file
  // 'epsilon' in for instance A^{epsilon}(x) = A(x,x/epsilon)
  const double epsilon_ = DSC_CONFIG_GET("problem.epsilon", 1.0f);
  const int refinement_level_ = DSC_CONFIG_GET("fem.grid_level", 4);
  out << "Log-File for Elliptic Model Problem " << Problem::name() << "." << std::endl << std::endl;
  if (Problem::getModelData()->linear())
    out << "Problem is declared as being LINEAR." << std::endl;
  else
    out << "Problem is declared as being NONLINEAR." << std::endl;

  if (info.hasExactSolution()) {
    out << "Exact solution is available." << std::endl << std::endl;
  } else {
    out << "Exact solution is not available." << std::endl << std::endl;
  }
  out << "Computations were made for:" << std::endl << std::endl;
  out << "Refinement Level for Grid = " << refinement_level_ << std::endl << std::endl;

  out << "Epsilon = " << epsilon_ << std::endl << std::endl;
}

void write_discrete_function(CommonTraits::DiscreteFunctionType &discrete_solution, const std::string prefix)
{
  discrete_solution.visualize(prefix, true, Dune::VTK::appendedbase64);
}

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {
