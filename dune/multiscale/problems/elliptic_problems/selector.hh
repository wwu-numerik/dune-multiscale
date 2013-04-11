#ifndef DUNE_MS_PROBLEMS_SELECTOR_HH
#define DUNE_MS_PROBLEMS_SELECTOR_HH

//for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "eight.hh"
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

namespace Problem {
//this pulls everything from the subnamespace into Problem and should only be done for the "active" problem
#ifndef PROBLEM_NAME
  #define PROBLEM_NAME Nine
#endif

using namespace Problem::PROBLEM_NAME;

#define STRIN(x) #x
#define STR(x) STRIN(x)
static const std::string name = std::string(STRIN(PROBLEM_NAME));
#undef STR
#undef STRIN
} //! @} namespace Problem

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
