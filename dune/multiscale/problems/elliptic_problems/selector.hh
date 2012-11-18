#ifndef DUNE_MS_PROBLEMS_SELECTOR_HH
#define DUNE_MS_PROBLEMS_SELECTOR_HH

//for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "easy.hh"
#include "eight.hh"
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
  using namespace Problem::Nine;
#else
  using namespace Problem::PROBLEM_NAME;
#endif
} // namespace Problem

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
