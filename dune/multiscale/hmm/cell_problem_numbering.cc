#include "cell_problem_numbering.hh"

namespace Dune {
namespace Multiscale {
namespace HMM {

CellProblemNumberingManager::CellProblemNumberingManager(const DiscreteFunctionSpaceType& discreteFunctionSpace)
{
  int counter = 0;
  int number_of_entity = 0;
  for (const auto& entity : discreteFunctionSpace)
  {
    EntityPointerType ep(entity);
    cell_numbering_map_NL_.insert( std::make_pair(ep, number_of_entity) );
    const int numBaseFunctions = discreteFunctionSpace.basisFunctionSet(entity).size();
    for (int i = 0; i < numBaseFunctions; ++i)
    {
      const std::pair< EntityPointerType, int > idPair(ep, i);
      cell_numbering_map_.insert( std::make_pair(idPair, counter) );
      counter++;
    }
    number_of_entity++;
  }
}

int CellProblemNumberingManager::get_number_of_cell_problem(const EntityPointerType& ent, const int& numOfBaseFunction) const {
  const typename CellNumMapType::key_type idPair(ent, numOfBaseFunction);
  const auto it = cell_numbering_map_.find(idPair);
  if (it != cell_numbering_map_.end() )
    return it->second;
  else
    DUNE_THROW(Dune::RangeError, "no number for entity");
}

int CellProblemNumberingManager::get_number_of_cell_problem(const EntityPointerType& ent) const {
  const auto it = cell_numbering_map_NL_.find(ent);
  if (it != cell_numbering_map_NL_.end() )
    return it->second;
  else
    DUNE_THROW(Dune::RangeError, "no number for entity");
}

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {
