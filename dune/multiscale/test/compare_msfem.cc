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

  set_param();

  const auto errors = algorithm();
  for (const auto& err : errors) {
    // fails iff 1e-2 <= error
    EXPECT_GT(double(1e-2), err.second);
  }
}

int main(int argc, char** argv) {
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}

void set_param() {
  DSC_CONFIG.set("grids.macro_cells_per_dim", "[5;5;5]");
  DSC_CONFIG.set("micro_cells_per_macrocell_dim", "[10;10;10]");
  DSC_CONFIG.set("msfem.oversampling_layers", 4);
  DSC_CONFIG.set("global.vtk_output", 0);
  DSC_CONFIG.set("problem.name", "Nine");
}
