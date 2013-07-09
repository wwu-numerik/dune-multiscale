// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_PROBLEMS_SELECTOR_HH
#define DUNE_MS_PROBLEMS_SELECTOR_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG


#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {

std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getFirstSource();
std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getExactSolution();
std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getSecondSource();
std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getMassTerm();
std::unique_ptr<Dune::Multiscale::CommonTraits::DiffusionType> getDiffusion();
std::unique_ptr<const Dune::Multiscale::CommonTraits::LowerOrderTermType> getLowerOrderTerm();
std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getDefaultDummyFunction();
std::unique_ptr<Dune::Multiscale::CommonTraits::ModelProblemDataType> getModelData();

std::string name();

} //! @} namespace Problem
} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
