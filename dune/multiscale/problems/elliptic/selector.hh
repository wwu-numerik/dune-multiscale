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

//for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "eight.hh"
#include "eleven.hh"
#include "example.hh"
#include "five.hh"
#include "four.hh"
#include "nine.hh"
#include "one.hh"
#include "seven.hh"
#include "six.hh"
#include "ten.hh"
#include "three.hh"
#include "toy.hh"
#include "two.hh"

#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {
namespace Problem {
//this pulls everything from the subnamespace into Problem and should only be done for the "active" problem
//using namespace Dune::Multiscale::Problem::PROBLEM_NAME;
typedef Dune::Multiscale::Problem::PROBLEM_NAME::Diffusion Diffusion;

typedef Dune::Multiscale::Problem::PROBLEM_NAME::DefaultDummyFunction DefaultDummyFunction;
typedef Dune::Multiscale::Problem::PROBLEM_NAME::MassTerm MassTerm;
typedef Dune::Multiscale::Problem::PROBLEM_NAME::LowerOrderTerm LowerOrderTerm;
typedef Dune::Multiscale::Problem::PROBLEM_NAME::ModelProblemData ModelProblemData;

#define STRIN(x) #x
#define STR(x) STRIN(x)
static const std::string name = std::string(STRIN(PROBLEM_NAME));
#undef STR
#undef STRIN

std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getFirstSource();
std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getExactSolution();
std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> getSecondSource();

} //! @} namespace Problem
} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
