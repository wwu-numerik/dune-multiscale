#include <config.h>

#include "localsolution_proxy.hh"

#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/multiscale/msfem/proxygridview.hh>

Dune::Multiscale::MsFEM::LocalsolutionProxy::LocalsolutionProxy(
    CorrectionsMapType&& corrections, const CommonTraits::DiscreteFunctionSpaceType& coarseSpace,
    const LocalGridList& gridlist)
  : BaseType(*corrections.begin()->second)
  , corrections_(std::move(corrections))
  , index_set_(coarseSpace.grid_view()->grid().leafIndexSet())
  , search_(DSC::make_unique<LocalGridSearch>(coarseSpace, gridlist))
  , gridlist_(gridlist) {
  assert(corrections_.size() == index_set_.size(0));
}

std::unique_ptr<Dune::Multiscale::MsFEM::LocalsolutionProxy::LocalFunctionType>
Dune::Multiscale::MsFEM::LocalsolutionProxy::local_function(const BaseType::EntityType& entity) const {
  const auto& coarse_cell = *search_->current_coarse_pointer();
  auto it = corrections_.find(index_set_.index(coarse_cell));
  if (it != corrections_.end())
    return it->second->local_function(entity);
  DUNE_THROW(InvalidStateException, "Coarse cell was not found!");
}

Dune::Multiscale::ProxyGridview Dune::Multiscale::MsFEM::LocalsolutionProxy::grid_view() const {
  return ProxyGridview(gridlist_);
}

Dune::Multiscale::MsFEM::LocalGridSearch& Dune::Multiscale::MsFEM::LocalsolutionProxy::search() const {
  assert(search_);
  return *search_;
}
