#include <dune/multiscale/test/test_common.hh>

#include <unordered_set>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>

template <typename T>
std::vector<typename T::GlobalCoordinate> corners(const T& geo) {
  std::vector<typename T::GlobalCoordinate> ret;
  for (auto c : DSC::valueRange(geo.corners())) {
    auto co = geo.corner(c);
    ret.push_back(co);
  }
  return ret;
}

template <typename T>
std::vector<typename T::GlobalCoordinate> cornersA(const T& geo) {
  std::vector<typename T::GlobalCoordinate> ret;
  auto rg = DSC::cornerRange(geo);
  auto end = rg.end();
  for (auto c = rg.begin(); c != end; ++c) {
    typename T::GlobalCoordinate p = *c;
    ret.push_back(p);
  }
  return ret;
}

struct GridMatch : public GridTestBase {

  void check_dimensions() {
    const auto dimensions = make_pair(DSG::dimensions(grids_.first->leafGridView()),
                                           DSG::dimensions(grids_.second->leafGridView()));
    const auto limits = make_pair(dimensions.first.coord_limits, dimensions.second.coord_limits);
    for ( auto d : DSC::valueRange(CommonTraits::dimDomain)) {
      EXPECT_DOUBLE_EQ(limits.first[d].min(), limits.second[d].min());
      EXPECT_DOUBLE_EQ(limits.first[d].max(), limits.second[d].max());
    }
    EXPECT_GT(dimensions.first.entity_volume.max(), dimensions.second.entity_volume.max());
    // second max not a typo
    EXPECT_GT(dimensions.first.entity_volume.min(), dimensions.second.entity_volume.max());
  }

  void check_unique_corners() {
    for(auto& grid : {grids_.first, grids_.second}) {
      for (const auto& ent : DSC::viewRange(grid->leafGridView())) {
        const auto& geo = ent.geometry();
        const auto cor = corners(geo);
        EXPECT_EQ(cor.size(), geo.corners());
        EXPECT_EQ(cor.size(), 4);
//        EXPECT_EQ(cor, cornersA(geo));
         for (auto i : DSC::valueRange(geo.corners()))
           for (auto j : DSC::valueRange(geo.corners()))
             if (i!=j) {
               EXPECT_TRUE(DSC::FloatCmp::vec_ne(geo.corner(i), geo.corner(j)));
               EXPECT_NE(geo.corner(i), geo.corner(j));
             } else {
               EXPECT_TRUE(DSC::FloatCmp::eq(geo.corner(i), geo.corner(j)));
               EXPECT_EQ(geo.corner(i), geo.corner(j));
             }
      }
    }
  }

//  void project() {
//    using namespace Dune::Multiscale::MsFEM;
//    CommonTraits::GdtSpaceType coarse_space(grids_.first->leafGridView());
//    LocalGridList subgrid_list(coarse_space);
//    LocalsolutionProxy::CorrectionsMapType local_corrections;
//    LocalGridSearch search(coarse_space, subgrid_list);
//    auto& coarse_indexset = coarse_space.grid_view()->grid().leafIndexSet();
//    LocalsolutionProxy proxy(local_corrections, coarse_indexset, search);

//    CommonTraits::DiscreteFunctionType fine_scale_part;
//    DS::MsFEMProjection::project(proxy, fine_scale_part, search);

//  }

  void check_inside() {
    for (const auto& ent : DSC::viewRange(grids_.second->leafGridView())) {
      const auto& geo = ent.geometry();
      for(auto corner : corners(ent.geometry())) {
        bool found = false;
        for (const auto& coarse_ent : DSC::viewRange(grids_.first->leafGridView())) {
          found = found || DSG::reference_element(coarse_ent).checkInside(geo.local(corner));
        }
        EXPECT_TRUE(found);
      }
    }
  }

  void check_lagrange_points() {

    for(auto& grid : {grids_.first, grids_.second}) {
      CommonTraits::GridProviderType grid_provider(grid);
      const CommonTraits::GdtSpaceType space =
          CommonTraits::SpaceProviderType::create(grid_provider, CommonTraits::st_gdt_grid_level);

      for (const auto& ent : DSC::viewRange(grid->leafGridView())) {
        const auto& geo = ent.geometry();
        const auto cor = corners(geo);
        const auto lp = DS::global_evaluation_points(space, ent);
        EXPECT_EQ(cor, lp);
      }
    }
  }

  void check_fine_lp_in_coarse() {

    for(auto& grid : {grids_.second}) {
      CommonTraits::GridProviderType grid_provider(grid);
      const CommonTraits::GdtSpaceType space =
          CommonTraits::SpaceProviderType::create(grid_provider, CommonTraits::st_gdt_grid_level);

      for (const auto& ent : DSC::viewRange(grid->leafGridView())) {
        const auto lg_points = DS::global_evaluation_points(space, ent);
        for(auto  lg : lg_points) {
          bool found = false;
          for (const auto& coarse_ent : DSC::viewRange(grids_.first->leafGridView())) {
            found = found || DSG::reference_element(coarse_ent).checkInside(ent.geometry().local(lg));
          }
          EXPECT_TRUE(found);
        }
      }
    }
  }

};

TEST_P(GridMatch, Match) {
  this->check_dimensions();
  this->check_unique_corners();
  this->check_inside();
  this->check_lagrange_points();
  this->check_fine_lp_in_coarse();
}


INSTANTIATE_TEST_CASE_P(
    TestName,
    GridMatch,
    testing::Values(p_small, p_large, p_aniso, p_wover, p_huge));

int main(int argc, char** argv) {
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}


