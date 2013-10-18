#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/ranges.hh>

#include "cell_problem_numbering.hh"

namespace Dune {
namespace Multiscale {
namespace HMM {

CellProblemNumberingManager::CellProblemNumberingManager(const DiscreteFunctionSpaceType& discreteFunctionSpace) {
  std::size_t counter = 0;
  std::size_t number_of_entity = 0;
  for (const auto& entity : discreteFunctionSpace) {
    EntityPointerType ep(entity);
    cell_numbering_map_NL_.insert(std::make_pair(ep, number_of_entity));
    for (auto i : DSC::valueRange(discreteFunctionSpace.basisFunctionSet(entity).size())) {
      const KeyType idPair(ep, i);
      cell_numbering_map_.insert(std::make_pair(idPair, counter));
      counter++;
    }
    number_of_entity++;
  }
}

std::size_t CellProblemNumberingManager::get_number_of_cell_problem(const EntityPointerType& ent,
                                                                    const std::size_t& numOfBaseFunction) const {
  const auto it = cell_numbering_map_.find(std::make_pair(ent, numOfBaseFunction));
  if (it != cell_numbering_map_.end())
    return it->second;
  else
    DUNE_THROW(Dune::RangeError, "no number for entity");
}

std::size_t CellProblemNumberingManager::get_number_of_cell_problem(const EntityPointerType& ent) const {
  const auto it = cell_numbering_map_NL_.find(ent);
  if (it != cell_numbering_map_NL_.end())
    return it->second;
  else
    DUNE_THROW(Dune::RangeError, "no number for entity");
}

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {
