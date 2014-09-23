#ifndef DUNE_MULTISCALE_TEST_COMMON_HH
#define DUNE_MULTISCALE_TEST_COMMON_HH

#include <config.h>

#include <dune/stuff/test/main.hxx>
#include <dune/stuff/test/gtest/gtest.h>

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
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
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
    DSC_CONFIG.set(vp.first, vp.second, true);
  }
}

class GridTestBase : public ::testing::TestWithParam<map<string,string>> {

public:

 GridTestBase() {
    set_param(GetParam());
    grids_ = make_grids();
  }
 virtual ~GridTestBase() {
   grids_ = decltype(grids_)();
 }

protected:
 decltype(make_grids()) grids_;
};


struct GridAndSpaces : public GridTestBase {
public:

 GridAndSpaces()
   : GridTestBase()
   , coarse_grid_provider(*grids_.first)
   , fine_grid_provider(*grids_.second)
   , coarseSpace(CommonTraits::SpaceProviderType::create(coarse_grid_provider, CommonTraits::st_gdt_grid_level))
   , fineSpace(CommonTraits::SpaceProviderType::create(fine_grid_provider, CommonTraits::st_gdt_grid_level))
 {}


protected:
  CommonTraits::GridProviderType coarse_grid_provider;
  CommonTraits::GridProviderType fine_grid_provider;
  const CommonTraits::SpaceType coarseSpace;
  const CommonTraits::SpaceType fineSpace;
};

static const map<string, string> p_small = {{"grids.macro_cells_per_dim", "[4 4 4]"}
                                           ,{"grids.micro_cells_per_macrocell_dim", "[8 8 8]"}
                                           ,{"msfem.oversampling_layers", "0"}
                                           };
static const map<string, string> p_large = {{"grids.macro_cells_per_dim", "[20 20 20]"}
                                           ,{"grids.micro_cells_per_macrocell_dim", "[40 40 40]"}
                                           ,{"msfem.oversampling_layers", "0"}
                                           };
static const map<string, string> p_aniso = {{"grids.macro_cells_per_dim", "[14 4 6]"}
                                           ,{"grids.micro_cells_per_macrocell_dim", "[3 32 8]"}
                                           ,{"msfem.oversampling_layers", "0"}
                                           };
static const map<string, string> p_wover = {{"grids.macro_cells_per_dim", "[14 14 14]"}
                                           ,{"grids.micro_cells_per_macrocell_dim", "[18 18 18]"}
                                           ,{"msfem.oversampling_layers", "6"}
                                           };
static const map<string, string> p_huge  = {{"grids.macro_cells_per_dim", "[40 40 40]"}
                                           ,{"grids.micro_cells_per_macrocell_dim", "[120 120 120]"}
                                           ,{"msfem.oversampling_layers", "0"}
                                           };
static const map<string, string> p_fail  = {{"grids.macro_cells_per_dim", "[12 15 10]"}
                                           ,{"grids.micro_cells_per_macrocell_dim", "[6 7 10]"}
                                           ,{"msfem.oversampling_layers", "0"}
                                           };
static const map<string, string> p_minimal  = {{"grids.macro_cells_per_dim", "[1 1 1]"}
                                           ,{"grids.micro_cells_per_macrocell_dim", "[1 1 1]"}
                                           ,{"msfem.oversampling_layers", "0"}
                                           };
#endif // DUNE_MULTISCALE_TEST_COMMON_HH
