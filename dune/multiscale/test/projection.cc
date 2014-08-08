#include <dune/multiscale/test/test_common.hh>

#include <dune/stuff/functions/global.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/l2.hh>
#include <functional>

using namespace Dune::GDT;

struct Projection : public GridAndSpaces {
  typedef DS::FunctionTypeGenerator<MsFEMTraits::LocalConstantFunctionType, DS::GlobalLambdaFunction>::type Lambda;

  LocalsolutionProxy::CorrectionsMapType fill_local_corrections(const Lambda& lambda) {
    LocalsolutionProxy::CorrectionsMapType local_corrections;
    for(auto& lc_pair : local_corrections) {
      auto& lc = *lc_pair.second;
      Operators::apply_projection(lambda, lc);
    }
    return local_corrections;
  }

  void project() {
    using namespace Dune::Multiscale::MsFEM;

    LocalGridList subgrid_list(coarseSpace);
    const double constant(1);
    Lambda lambda([&](CommonTraits::DomainType /*x*/) { return constant;}, 0 );
    auto local_corrections = fill_local_corrections(lambda);
    LocalGridSearch search(coarseSpace, subgrid_list);
    auto& coarse_indexset = coarseSpace.grid_view()->grid().leafIndexSet();
    LocalsolutionProxy proxy(local_corrections, coarse_indexset, search);

    CommonTraits::DiscreteFunctionType fine_scale_part(fineSpace);
    DS::MsFEMProjection::project(proxy, fine_scale_part, search);

    const auto norm = std::sqrt(Dune::GDT::Products::L2< CommonTraits::GridViewType >(*(fineSpace.grid_view()))
                                    .induced_norm(fine_scale_part));
    EXPECT_DOUBLE_EQ(norm, constant);
  }
};

TEST_P(Projection, Project) {
  this->project();
}

static const auto common_values = testing::Values(p_small, p_large, p_aniso, p_wover, /*p_huge,*/ p_fail);
//static const auto common_values = testing::Values(p_large);

INSTANTIATE_TEST_CASE_P( TestNameB, Projection, common_values);

int main(int argc, char** argv) {
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
