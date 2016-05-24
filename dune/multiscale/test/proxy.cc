#include <dune/multiscale/test/test_common.hxx>

#include <unordered_set>
#include <dune/xt/common/float_cmp.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/common/df_io.hh>

struct PointsAndStuff : public GridAndSpaces {

  void check_fine_lp_in_coarse() {
    const auto clearGuard = Dune::Multiscale::DiscreteFunctionIO::clear_guard();
    LocalGridList localgrid_list(*problem_, coarseSpace);
    LocalGridSearch search(coarseSpace, localgrid_list);
//    auto proxy_view = msfem_solution_->grid_view();
  }

};

TEST_F(PointsAndStuff, LP) {
  this->check_fine_lp_in_coarse();
}

