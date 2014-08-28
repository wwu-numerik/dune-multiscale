#include <config.h>

#include "localsolution_proxy.hh"

#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>

Dune::Multiscale::MsFEM::LocalsolutionProxy::LocalsolutionProxy(const CorrectionsMapType &corrections,
                                                                const LeafIndexSetType &index_set,
                                                                const LocalGridSearch & search)
  : BaseType(*corrections.begin()->second)
  , corrections_(corrections)
  , index_set_(index_set)
  , search_(search) {
  assert(corrections.size() == index_set.size(0));
}


std::unique_ptr<Dune::Multiscale::MsFEM::LocalsolutionProxy::LocalFunctionType>
Dune::Multiscale::MsFEM::LocalsolutionProxy::local_function(const BaseType::EntityType &entity) const {
  const auto& coarse_cell = *search_.current_coarse_pointer();
  auto it = corrections_.find(index_set_.index(coarse_cell));
  if (it != corrections_.end())
    return it->second->local_function(entity);
  DUNE_THROW(InvalidStateException, "Coarse cell was not found!");
}
