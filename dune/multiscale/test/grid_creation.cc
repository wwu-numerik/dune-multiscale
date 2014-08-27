#include <dune/multiscale/test/test_common.hh>

#include <unordered_set>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/common/configuration.hh>
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

struct PointsAndStuff : public GridAndSpaces {

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
    MsFEM::LocalGridList subgrid_list(coarseSpace);
    MsFEM::LocalGridSearch search(coarseSpace, subgrid_list);

    for (const auto& ent : fineSpace) {
      const auto lg_points = DS::global_evaluation_points(fineSpace, ent);
      for(auto  lg : lg_points) {
        bool found = false;
        for (const auto& coarse_ent : DSC::viewRange(grids_.first->leafGridView())) {
          found = found || DSG::reference_element(coarse_ent).checkInside(ent.geometry().local(lg));
        }
        EXPECT_TRUE(found);
      }
    }
  }

  void check_search() {
    MsFEM::LocalGridList subgrid_list(coarseSpace);
    MsFEM::LocalGridSearch search(coarseSpace, subgrid_list);

    for (const auto& ent : fineSpace) {
      const auto lg_points = DS::global_evaluation_points(fineSpace, ent);
      const auto evaluation_entity_ptrs = search(lg_points);
      EXPECT_GE(evaluation_entity_ptrs.size(), lg_points.size());
    }
  }


  void check_local_grids() {
    MsFEM::LocalGridList subgrid_list(coarseSpace);
    EXPECT_EQ(subgrid_list.size(), grids_.first->size(0));
  }


};

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
    const auto dim_world = CommonTraits::GridType::dimensionworld;
    const auto sub = DSC_CONFIG.sub("grids");
    const auto microPerMacro = sub.get<CommonTraits::DomainType>("micro_cells_per_macrocell_dim", CommonTraits::DomainType(-1), dim_world);
    const auto coarse_cells = sub.get<CommonTraits::DomainType>("macro_cells_per_dim", CommonTraits::DomainType(-1), dim_world);
    for(auto c: coarse_cells)
      EXPECT_GT(c, 0);
    for(auto c: microPerMacro)
      EXPECT_GT(c, 0);
    long expected_coarse = 1;
    for(auto i : DSC::valueRange(dim_world))
      expected_coarse *= coarse_cells[i];
    auto expected_fine = expected_coarse;
    for(auto i : DSC::valueRange(dim_world))
      expected_fine *= microPerMacro[i];
    EXPECT_EQ(grids_.first->leafGridView().size(0), expected_coarse);
    EXPECT_EQ(grids_.second->leafGridView().size(0), expected_fine);

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

  void check_inside() {
    for (const auto& ent : DSC::viewRange(grids_.second->leafGridView())) {
      for(auto corner : corners(ent.geometry())) {
        bool found = false;
        for (const auto& coarse_ent : DSC::viewRange(grids_.first->leafGridView())) {
          found = found || DSG::reference_element(coarse_ent).checkInside(coarse_ent.geometry().local(corner));
        }
        EXPECT_TRUE(found);
      }
    }
  }

};

TEST_P(GridMatch, Match) {
  this->check_dimensions();
  this->check_unique_corners();
  this->check_inside();
}

TEST_P(PointsAndStuff, LP) {
  this->check_lagrange_points();
  this->check_fine_lp_in_coarse();
  this->check_local_grids();
  this->check_search();
}

static const auto common_values = testing::Values(p_small, p_large, p_aniso, p_wover, /*p_huge,*/ p_fail);
//static const auto common_values = testing::Values(p_large);

INSTANTIATE_TEST_CASE_P( TestNameA, GridMatch, common_values);
INSTANTIATE_TEST_CASE_P( TestNameB, PointsAndStuff, common_values);

int main(int argc, char** argv) {
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}


