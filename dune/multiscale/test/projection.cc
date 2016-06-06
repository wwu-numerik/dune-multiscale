#include <dune/multiscale/test/test_common.hxx>

#include <dune/stuff/functions/global.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/l2.hh>
#include <functional>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/common/heterogenous.hh>
#include <dune/multiscale/common/df_io.hh>

using namespace Dune::GDT;

struct Projection : public GridAndSpaces {

  typedef DS::FunctionTypeGenerator<MsFEMTraits::LocalConstantFunctionType, DS::GlobalLambdaFunction>::type Lambda;

  LocalsolutionProxy::CorrectionsMapType fill_local_corrections(const Lambda& lambda,
                                                                const LocalGridList& localgrid_list) {
    LocalsolutionProxy::CorrectionsMapType local_corrections;
    for (const auto& coarse_entity : Dune::elements(coarseSpace.grid_view())) {
      LocalSolutionManager localSolManager(coarseSpace, coarse_entity, localgrid_list);
      auto& coarse_indexset = coarseSpace.grid_view().grid().leafIndexSet();
      const auto coarse_index = coarse_indexset.index(coarse_entity);
      local_corrections[coarse_index] =
          Dune::XT::Common::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(localSolManager.space(), " ");
      Dune::GDT::project(lambda, *local_corrections[coarse_index]);
    }
    return local_corrections;
  }

  void project() {
    const auto clearGuard = Dune::Multiscale::DiscreteFunctionIO::clear_guard();
    LocalGridList localgrid_list(*problem_, coarseSpace);
    const double constant(1);
    Lambda lambda([&](CommonTraits::DomainType /*x*/) { return constant; }, 0);
    auto local_corrections = fill_local_corrections(lambda, localgrid_list);

    LocalsolutionProxy proxy(std::move(local_corrections), coarseSpace, localgrid_list);

    CommonTraits::DiscreteFunctionType fine_scale_part(fineSpace);
    MsFEMProjection::project(proxy, fine_scale_part);

    const auto norm = std::sqrt(
        Dune::GDT::Products::L2<CommonTraits::GridViewType>(fineSpace.grid_view()).induced_norm(fine_scale_part));
    EXPECT_DOUBLE_EQ(norm, constant);
  }
};

struct Search : public GridAndSpaces {

  typedef DS::FunctionTypeGenerator<MsFEMTraits::LocalConstantFunctionType, DS::GlobalLambdaFunction>::type Lambda;

  void lg_search() {
    LocalGridList localgrid_list(*problem_, coarseSpace);
    LocalGridSearch lgs(coarseSpace, localgrid_list);

    for (auto&& i : Dune::XT::Common::value_range(coarseSpace.grid_view().size(0))) {
      auto&& lg = localgrid_list.getSubGrid(i);
      const auto lg_view = lg.leafGridView();
      for (auto&& lg_ent : Dune::elements(lg_view)) {
        auto center = lg_ent.geometry().center();
        lgs({center});
      }
    }
  }
};

TEST_F(Projection, Project) { this->project(); }
// TEST_P(Search, Project) {
//  this->lg_search();
//}
