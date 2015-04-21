// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_SRC_COMMON_HH
#define DUNE_MULTISCALE_SRC_COMMON_HH

namespace tbb {
class tbb_exception;
}

namespace Dune {

class Exception;

namespace Multiscale {

//! setup code common to fem/msfem/hmm
void init(int argc, char** argv); // init

int handle_exception(const Dune::Exception& exp);
int handle_exception(const std::exception& exp);
int handle_exception(const tbb::tbb_exception& exp);

void mem_usage();
void dump_environment();

} // namespace Dune {
} // namespace Multiscale {

#endif // DUNE_MULTISCALE_SRC_COMMON_HH
