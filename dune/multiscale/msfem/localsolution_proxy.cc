#include <config.h>

#include "localsolution_proxy.hh"

#include <dune/stuff/common/profiler.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/proxygridview.hh>
#include <dune/gdt/operators/prolongations.hh>

Dune::Multiscale::LocalsolutionProxy::LocalsolutionProxy(CorrectionsMapType&& corrections,
                                                         const CommonTraits::SpaceType& coarseSpace,
                                                         const LocalGridList& gridlist)
  : BaseType(*corrections.begin()->second)
  , corrections_(std::move(corrections))
  , view_(coarseSpace.grid_view())
  , index_set_(view_.grid().leafIndexSet())
  , search_(DSC::make_unique<LocalGridSearch>(coarseSpace, gridlist)) {
  assert(corrections_.size() == index_set_.size(0));
}

std::unique_ptr<Dune::Multiscale::LocalsolutionProxy::LocalFunctionType>
Dune::Multiscale::LocalsolutionProxy::local_function(const BaseType::EntityType& entity) const {
  const auto& coarse_cell = *search_->current_coarse_pointer();
  auto it = corrections_.find(index_set_.index(coarse_cell));
  if (it != corrections_.end())
    return it->second->local_function(entity);
  DUNE_THROW(InvalidStateException, "Coarse cell was not found!");
}

void
Dune::Multiscale::LocalsolutionProxy::add(const Dune::Multiscale::CommonTraits::DiscreteFunctionType& coarse_func) {
  DSC::ScopedTiming st("proxy.add");
  CorrectionsMapType targets;
  for (auto& cr : corrections_) {
    targets[cr.first] =
        DSC::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(cr.second->space(), "tmpcorrection");
  }
  for (auto& range_pr : targets) {
    auto id = range_pr.first;
    auto& range = *range_pr.second;
    const Dune::GDT::Operators::LagrangeProlongation<MsFEMTraits::LocalGridViewType> prolongation_operator(
        range.space().grid_view());
    prolongation_operator.apply(coarse_func, range);
    auto& correction = *corrections_[id];
    correction.vector() += range.vector();
  }
}

Dune::Multiscale::LocalGridSearch& Dune::Multiscale::LocalsolutionProxy::search() const {
  assert(search_);
  return *search_;
}
