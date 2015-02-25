#include <dune/multiscale/test/test_common.hxx>

#include <dune/stuff/common/fixed_map.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/configuration.hh>

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/fixed_map.hh>

#include <boost/filesystem.hpp>

//          (Config key -> value, error name -> max value)
typedef pair< map<string,string>, FixedMap<string,double,2>> TestArgs;


const TestArgs m_small{p_small, {{"msfem_exact_L2", 0.251}, {"msfem_exact_H1s", 2.67}}};
const TestArgs m_large{p_large, {{"msfem_exact_L2", 0.07}, {"msfem_exact_H1s", 1.15}}};
const TestArgs m_minimal{p_minimal, {{"msfem_exact_L2", 0.57}, {"msfem_exact_H1s", 2.3}}};

struct MsFemCompare : public ::testing::TestWithParam<TestArgs> {

  MsFemCompare()
  {
    set_param(this->GetParam().first);
    set_param({{"msfem.fem_comparison", "0"},
        {"global.vtk_output", "0"},
        {"problem.name", "Nine"}});
  }
};

TEST_P(MsFemCompare, All) {
  using namespace Dune::Multiscale;
  const auto errorsMap = msfem_algorithm();
  
  const auto expected_errors = this->GetParam().second;
  const auto found = errorsMap.find("msfem_exact_L2");
  const auto end = errorsMap.end();
  EXPECT_NE(found, end);
  EXPECT_GT(expected_errors["msfem_exact_L2"], found->second);
  const auto found2 = errorsMap.find("msfem_exact_H1s");
  EXPECT_NE(found2, end);
  EXPECT_GT(expected_errors["msfem_exact_H1s"], found->second);
}

//template < int foo >
//struct VC {
//  static auto values() -> decltype(testing::Values(/*m_small, */m_large)) {
//    return testing::Values(/*m_small,*/ m_large);
//  }
//};

template < int foo >
struct VC {
  static auto values() -> decltype(testing::Values(/*m_small, */m_large)) {
    return testing::Values(m_large);
  }
};

template <>
struct VC<3> {
  static auto values() -> decltype(testing::Values(m_small)) {
    return testing::Values(m_small);
  }
};

INSTANTIATE_TEST_CASE_P( MsFemComparisons, MsFemCompare, VC<CommonTraits::world_dim>::values());


