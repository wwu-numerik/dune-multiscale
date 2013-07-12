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

typedef std::unique_ptr<const CommonTraits::FunctionBaseType> BasePtr;
BasePtr getFirstSource();
BasePtr getExactSolution();
BasePtr getSecondSource();
BasePtr getMassTerm();
std::unique_ptr<const CommonTraits::DiffusionType> getDiffusion();
std::unique_ptr<const CommonTraits::LowerOrderTermType> getLowerOrderTerm();
std::unique_ptr<const CommonTraits::DirichletBCType> getDirichletBC();
std::unique_ptr<const CommonTraits::NeumannBCType> getNeumannBC();
BasePtr getDefaultDummyFunction();
std::unique_ptr<const CommonTraits::ModelProblemDataType> getModelData();

std::string name();

} //! @} namespace Problem
} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
