#include <config.h>

#include "localsolution_proxy.hh"

#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>

Dune::Multiscale::MsFEM::LocalsolutionProxy::LocalsolutionProxy(const Dune::Multiscale::MsFEM::LocalsolutionProxy::CorrectionsMapType &corrections, const Dune::Multiscale::MsFEM::LocalsolutionProxy::LeafIndexSetType &index_set, Dune::Multiscale::MsFEM::LocalGridSearch &search)
  : BaseType(*corrections.begin()->second)
  , corrections_(corrections)
  , index_set_(index_set)
  , search_(search) {}


std::unique_ptr<Dune::Multiscale::MsFEM::LocalsolutionProxy::LocalFunctionType> Dune::Multiscale::MsFEM::LocalsolutionProxy::local_function(const BaseType::EntityType &entity) const {
  const auto& coarse_cell = *search_.current_coarse_pointer();
  auto it = corrections_.find(index_set_.index(coarse_cell));
  if (it != corrections_.end())
    return it->second->local_function(entity);
  DUNE_THROW(InvalidStateException, "Coarse cell was not found!");
}
