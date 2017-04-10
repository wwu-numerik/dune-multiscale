#include <dune/multiscale/test/test_common.hxx>

#include <unordered_set>
#include <dune/xt/common/float_cmp.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/multiscale/common/error_calc.hh>
#include <dune/multiscale/common/df_io.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/fem_solver.hh>

using namespace Dune;
using namespace Dune::Multiscale;

struct ErrorCheck : public GridAndSpaces
{

  void run_error_calc()
  {
    const auto clearGuard = Dune::Multiscale::DiscreteFunctionIO::clear_guard();
    LocalsolutionProxy::CorrectionsMapType local_corrections;
    auto& coarse_space = this->coarseSpace;
    LocalGridList localgrid_list(*problem_, coarse_space);
    const CommonTraits::InteriorGridViewType interior =
        coarse_space.grid_view().grid().leafGridView<InteriorBorder_Partition>();
    const auto& coarse_indexset = coarse_space.grid_view().grid().leafIndexSet();
    for (const auto& coarse_entity : Dune::elements(interior)) {
      LocalproblemSolutionManager localSolManager(coarse_space, coarse_entity, localgrid_list);

      const auto coarse_index = coarse_indexset.index(coarse_entity);
      local_corrections[coarse_index] = Dune::XT::Common::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(
          localSolManager.space(), "correction");
    }
    const auto msfem_solution =
        Dune::XT::Common::make_unique<LocalsolutionProxy>(std::move(local_corrections), coarse_space, localgrid_list);
    ErrorCalculator ec(*problem_, msfem_solution);
    auto errors = ec.print(DSC_LOG_INFO_0);
    return;
  }
};

TEST_F(ErrorCheck, LP)
{
  this->run_error_calc();
}
