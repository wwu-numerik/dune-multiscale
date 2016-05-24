#ifndef DUNE_MULTISCALE_TEST_COMMON_HH
#define DUNE_MULTISCALE_TEST_COMMON_HH

#include <config.h>

#include <dune/xt/common/test/main.hxx>
#include <dune/xt/common/test/gtest/gtest.h>

#include <string>
#include <array>
#include <initializer_list>
#include <vector>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/common/configuration.hh>

#include <unordered_set>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/common/grid_creation.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/configuration.hh>

#include <boost/filesystem.hpp>


using namespace Dune::Stuff::Common;
using namespace Dune::Multiscale;
using namespace std;


string prepend_test_dir(string fn)
{
  boost::filesystem::path path(st_testdata_directory);
  path /= fn;
  return path.string();
}

void set_param(map<string, string> params) {
  for( auto vp : params) {
    DXTC_CONFIG.set(vp.first, vp.second, true);
  }
}

class GridTestBase : public ::testing::Test {

public:

 GridTestBase() {
    problem_ = DSC::make_unique<DMP::ProblemContainer>(Dune::MPIHelper::getCommunicator(), Dune::MPIHelper::getCommunicator(), DXTC_CONFIG);
    grids_ = make_grids(*problem_);
    DXTC_LOG_DEBUG << "Instantiating tests for dimension " << CommonTraits::world_dim << std::endl;
    constexpr bool sp_grid = std::is_same<CommonTraits::GridType,
        Dune::SPGrid<double, CommonTraits::world_dim, Dune::SPIsotropicRefinement>>::value;
    if (sp_grid && DXTC_CONFIG_GET("threading.max_count", 1) > 1) {
      DUNE_THROW(Dune::InvalidStateException, "SPGRID currently fails with > 1 threads");
    }
  }
 virtual ~GridTestBase() {
 }

protected:
 std::unique_ptr<DMP::ProblemContainer> problem_;
 std::pair<std::shared_ptr<CommonTraits::GridType>, std::shared_ptr<CommonTraits::GridType>> grids_;
};


struct GridAndSpaces : public GridTestBase {
public:

 GridAndSpaces()
   : GridTestBase()
   , coarseSpace(CommonTraits::SpaceChooserType::make_space(*grids_.first))
   , fineSpace(CommonTraits::SpaceChooserType::make_space(*grids_.second))
 {}


protected:
  const CommonTraits::SpaceType coarseSpace;
  const CommonTraits::SpaceType fineSpace;
};



//static const auto default_common_values = CommonTraits::world_dim < 3
//                                  // Values need to have number of elements
//#ifndef NDEBUG
//                                  ? testing::Values(p_small, p_small_aniso, p_small_wover)
//                                  : testing::Values(p_small, p_small_aniso, p_small_wover);
//#else
//                                  ? testing::Values(p_large, p_aniso, p_wover)
//                                  : testing::Values(p_small, p_aniso, p_wover);
//#endif

#endif // DUNE_MULTISCALE_TEST_COMMON_HH
