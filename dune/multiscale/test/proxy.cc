#include <dune/multiscale/test/test_common.hxx>

#include <unordered_set>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>


struct PointsAndStuff : public GridAndSpaces {

  void check_fine_lp_in_coarse() {
    MsFEM::LocalGridList subgrid_list(coarseSpace);
    MsFEM::LocalGridSearch search(coarseSpace, subgrid_list);
    auto proxy_view = msfem_solution_->grid_view();
  }

};

TEST_P(PointsAndStuff, LP) {
  this->check_fine_lp_in_coarse();
}

static const auto common_values = CommonTraits::world_dim < 3
                                  // Values need to have number of elements
                                  ? testing::Values(p_small/*, p_large, p_aniso, p_wover, p_fail*/)
                                  : testing::Values(p_small/*, p_minimal, p_minimal, p_minimal, p_minimal*/);

INSTANTIATE_TEST_CASE_P( TestNameB, PointsAndStuff, common_values);



