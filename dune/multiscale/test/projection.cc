#include <dune/multiscale/test/test_common.hh>


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
};

TEST_P(GridMatch, Match) {

  this->check_dimensions();

}


INSTANTIATE_TEST_CASE_P(
    TestName,
    GridMatch,
    testing::Values(p_small, p_large, p_aniso, p_wover, p_huge));

int main(int argc, char** argv) {
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}


