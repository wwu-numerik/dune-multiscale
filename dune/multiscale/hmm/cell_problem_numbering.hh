// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNEMS_HMM_CELL_NUMBERING_HH
#define DUNEMS_HMM_CELL_NUMBERING_HH

#include <utility>
#include <map>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/hmm/entity_compare.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {

//! only for the combination entity + number of local base function on entity
class CellProblemNumberingManager
{
private:
  typedef Multiscale::CommonTraits::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType                  GridType;

  typedef typename GridType::Codim< 0 >::Entity         EntityType;
  typedef typename GridType::Codim< 0 >::EntityPointer EntityPointerType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::IteratorType        IteratorType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef classcomp< GridPartType, DomainType, EntityPointerType > CompClass;

  typedef std::map< std::pair< EntityPointerType, int >, int, CompClass > CellNumMapType;

  // for the comparison of two entities:
  typedef entity_compare< GridPartType, DomainType, EntityPointerType > CompEntityClass;

  typedef std::map< EntityPointerType, int, CompEntityClass > CellNumMapNLType;

  CellNumMapType cell_numbering_map_;
  CellNumMapNLType cell_numbering_map_NL_;

public:
  /** simpliefied: in general we need CellNumMapType for the cell problem numering in the linear setting (entity and
   * local number of base function) and in the nonlinear case we need CellNumMapNLType (NL stands for nonlinear).
   * CellNumMapType is also required in the nonlinear case if we use test function reconstruction (TFR)
   **/
  explicit CellProblemNumberingManager(const DiscreteFunctionSpaceType& discreteFunctionSpace);

  //! use 'cp_num_manager.get_number_of_cell_problem( it, i )'
  int get_number_of_cell_problem(const EntityPointerType& ent, const int& numOfBaseFunction) const;

  /** use 'cp_num_manager.get_number_of_cell_problem( it )'
   * \attention 'get_number_of_cell_problem( it )' is NOT equal to 'get_number_of_cell_problem( it , 0 )'!
   **/
  int get_number_of_cell_problem(const EntityPointerType& ent) const;
};

} //namespace Dune {

#endif // DUNEMS_HMM_CELL_NUMBERING_HH
