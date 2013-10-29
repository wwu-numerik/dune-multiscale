#include <config.h>
#include <config.h>

#include <dune/stuff/test/test_common.hh>
#include <gtest.h>

#include <string>
#include <array>
#include <initializer_list>
#include <vector>
#include <dune/stuff/common/fixed_map.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/ranges.hh>

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/common/error_container.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/multiscale/problems/selector.hh>

#include <boost/filesystem.hpp>

using namespace Dune::Stuff::Common;

void set_param();

std::string prepend_test_dir(std::string fn)
{
  boost::filesystem::path path(st_testdata_directory);
  path /= fn;
  return path.string();
}

TEST(MSFEM, All) {
  using namespace Dune::Multiscale;
  using namespace Dune::Multiscale::MsFEM;

  int coarse_grid_level_ = 3;
  int number_of_layers_ = 4;
  int total_refinement_level_ = 5;
  const std::string macroGridName = prepend_test_dir("compare_msfem.dgf");
  set_param();

  //! ---------------------- local error indicators --------------------------------
  // ----- local error indicators (for each coarse grid element T) -------------
  const int max_loop_number = 10;
  ErrorContainer errors(max_loop_number);
  unsigned int loop_number = 0;
  while (algorithm(macroGridName, loop_number++, total_refinement_level_, coarse_grid_level_, number_of_layers_,
                   errors.locals, errors.totals, errors.total_estimated_H1_error_)) {
  }

  for (const auto total : errors.totals) {
    for (auto error : *total) {
      // fails iff 1e-2 <= error
      EXPECT_GT(double(1e-2), error);
    }
  }
}

int main(int argc, char** argv) {
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}

void set_param() {}
