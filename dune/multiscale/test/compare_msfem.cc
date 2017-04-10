#include <dune/multiscale/test/test_common.hxx>

#include <dune/xt/common/fixed_map.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/configuration.hh>

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/xt/common/fixed_map.hh>

#include <boost/filesystem.hpp>

TEST(MsFemCompare, All)
{
  using namespace Dune::Multiscale;
  auto errorsMap = msfem_algorithm();

  EXPECT_GT(DXTC_CONFIG.get("expected_errors.msfem_exact_L2", -1.), errorsMap["msfem_exact_L2"]);
  EXPECT_GT(DXTC_CONFIG.get("expected_errors.msfem_exact_H1s", -1.), errorsMap["msfem_exact_H1s"]);

  auto second_run = msfem_algorithm();
  for (auto error_pair : errorsMap) {
    EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(error_pair.second, second_run[error_pair.first]))
        << "FloatCmp " << error_pair.second << "\n      != " << second_run[error_pair.first];
  }
}
