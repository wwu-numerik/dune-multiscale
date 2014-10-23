#include <dune/multiscale/test/test_common.hxx>

#include <dune/stuff/functions/global.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/l2.hh>
#include <functional>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>

using namespace Dune::GDT;

struct Projection : public GridAndSpaces {
  typedef DS::FunctionTypeGenerator<MsFEMTraits::LocalConstantFunctionType, DS::GlobalLambdaFunction>::type Lambda;

  LocalsolutionProxy::CorrectionsMapType fill_local_corrections(const Lambda& lambda,
                                                                const LocalGridList& localgrid_list) {
    LocalsolutionProxy::CorrectionsMapType local_corrections;
    for (auto& coarse_entity : DSC::entityRange(coarseSpace.grid_view())) {
      LocalSolutionManager localSolManager(coarseSpace, coarse_entity, localgrid_list);
      auto& coarse_indexset = coarseSpace.grid_view().grid().leafIndexSet();
      const auto coarse_index = coarse_indexset.index(coarse_entity);
      local_corrections[coarse_index] =
          DSC::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(localSolManager.space(), " ");
      Operators::apply_projection(lambda, *local_corrections[coarse_index]);
    }
    return local_corrections;
  }

  void project() {
    LocalGridList localgrid_list(coarseSpace);
    const double constant(1);
    Lambda lambda([&](CommonTraits::DomainType /*x*/) { return constant;}, 0 );
    auto local_corrections = fill_local_corrections(lambda, localgrid_list);
    LocalGridSearch search(coarseSpace, localgrid_list);
    LocalsolutionProxy proxy(std::move(local_corrections), coarseSpace, localgrid_list);

    CommonTraits::DiscreteFunctionType fine_scale_part(fineSpace);
    DS::MsFEMProjection::project(proxy, fine_scale_part, search);

    const auto norm = std::sqrt(Dune::GDT::Products::L2< CommonTraits::GridViewType >(fineSpace.grid_view())
                                    .induced_norm(fine_scale_part));
    EXPECT_DOUBLE_EQ(norm, constant);
  }
};

TEST_P(Projection, Project) {
  this->project();
}


static const auto common_values = CommonTraits::world_dim < 3
                                  // Values need to have number of elements
                                  ? testing::Values(p_small/*, p_large, p_aniso, p_wover, p_fail*/)
                                  : testing::Values(p_small/*, p_minimal, p_minimal, p_minimal, p_minimal*/);

INSTANTIATE_TEST_CASE_P( TestNameB, Projection, common_values);

