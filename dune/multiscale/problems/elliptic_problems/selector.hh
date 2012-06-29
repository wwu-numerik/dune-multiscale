#ifndef DUNE_MS_PROBLEMS_SELECTOR_HH
#define DUNE_MS_PROBLEMS_SELECTOR_HH

#include "one.hh"
#include "two.hh"

namespace Problem {
//this pulls everything from the subnamespace into Problem and should only be done for the "active" problem
using namespace Problem::One;
} // namespace Problem

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
