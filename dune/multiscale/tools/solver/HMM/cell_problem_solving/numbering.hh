#ifndef DUNEMS_HMM_CELL_NUMBERING_HH
#define DUNEMS_HMM_CELL_NUMBERING_HH

#include <utility>
#include <map>

#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {
//! comparison class for the CellProblemNumberingManager:
template< class GridPartType, class DomainType, class EntityPointerType >
struct classcomp
{
  bool operator()(const std::pair< EntityPointerType, int >& left_entity_pair,
                  const std::pair< EntityPointerType, int >& right_entity_pair) const {
    // compare the barycenteres of the entities with the lexicographic order, than compare the int's (number of local
    // base function)

    typedef Fem::CachingQuadrature< GridPartType, 0 > Quadrature;

    // ------ right element

    const typename EntityPointerType::Entity::Geometry& geometry_right = ( *(right_entity_pair.first) ).geometry();

    Quadrature quadrature_right( ( *(right_entity_pair.first) ), 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_right = quadrature_right.point(0);

    DomainType barycenter_right_entity = geometry_right.global(local_point_right);

    // ------ left element

    const typename EntityPointerType::Entity::Geometry& geometry_left = ( *(left_entity_pair.first) ).geometry();

    Quadrature quadrature_left( ( *(left_entity_pair.first) ), 0 );

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_left = quadrature_left.point(0);

    DomainType barycenter_left_entity = geometry_left.global(local_point_left);

    enum { dimension = GridPartType::GridType::dimension };

    int current_axis = dimension - 1;

    while (current_axis >= 0)
    {
      if (barycenter_left_entity[current_axis] < barycenter_right_entity[current_axis])
      { return true; } else if (barycenter_left_entity[current_axis] > barycenter_right_entity[current_axis])
      { return false; }

      current_axis -= 1;
    }

    if (left_entity_pair.second < right_entity_pair.second)
    {
      return true;
    } else
    { return false; }

    return true;
  } // ()
};

//! comparison class for the CellProblemNumberingManager (just comparison of two entities!)
template< class GridPartType, class DomainType, class EntityPointerType >
struct entity_compare
{
  bool operator()(EntityPointerType left_entity,
                  EntityPointerType right_entity) const {
    // compare the barycenteres of the entities with the lexicographic order

    typedef Fem::CachingQuadrature< GridPartType, 0 > Quadrature;

    // ------ right element

    const typename EntityPointerType::Entity::Geometry& geometry_right = (*right_entity).geometry();

    Quadrature quadrature_right(*right_entity, 0);

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_right = quadrature_right.point(0);

    DomainType barycenter_right_entity = geometry_right.global(local_point_right);

    // ------ left element

    const typename EntityPointerType::Entity::Geometry& geometry_left = (*left_entity).geometry();

    Quadrature quadrature_left(*left_entity, 0);

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_left = quadrature_left.point(0);

    DomainType barycenter_left_entity = geometry_left.global(local_point_left);

    enum { dimension = GridPartType::GridType::dimension };

    int current_axis = dimension - 1;

    while (current_axis >= 0)
    {
      if (barycenter_left_entity[current_axis] < barycenter_right_entity[current_axis])
      { return true; } else if (barycenter_left_entity[current_axis] > barycenter_right_entity[current_axis])
      { return false; }

      current_axis -= 1;
    }

    return false;
  } // ()
};

//! only for the combination entity + number of local base function on entity
template< class DiscreteFunctionSpaceType >
class CellProblemNumberingManager
{
public:
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType                  GridType;

  typedef typename GridType::template Codim< 0 >::Entity         EntityType;
  typedef typename GridType::template  Codim< 0 >::EntityPointer EntityPointerType;

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

  /** simpliefied: in general we need CellNumMapType for the cell problem numering in the linear setting (entity and
   * local number of base function) and in the nonlinear case we need CellNumMapNLType (NL stands for nonlinear).
   * CellNumMapType is also required in the nonlinear case if we use test function reconstruction (TFR)
   **/
  inline explicit CellProblemNumberingManager(DiscreteFunctionSpaceType& discreteFunctionSpace)
  {
    int counter = 0;
    int number_of_entity = 0;

    const IteratorType endit = discreteFunctionSpace.end();
    for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
    {
      cell_numbering_map_NL_.insert( std::make_pair(EntityPointerType(*it), number_of_entity) );

      const BaseFunctionSetType baseSet
        = discreteFunctionSpace.baseFunctionSet(*it);

      // number of base functions on entity
      const int numBaseFunctions = baseSet.size();

      for (int i = 0; i < numBaseFunctions; ++i)
      {
        std::pair< EntityPointerType, int > idPair(EntityPointerType(*it), i);
        cell_numbering_map_.insert( std::make_pair(idPair, counter) );
        counter++;
      }
      number_of_entity++;
    }
  }

  //! use 'cp_num_manager.get_number_of_cell_problem( it, i )'
  inline int get_number_of_cell_problem(const EntityPointerType& ent, const int& numOfBaseFunction) const {
    const typename CellNumMapType::key_type idPair(ent, numOfBaseFunction);
    auto it = cell_numbering_map_.find(idPair);
    if (it != cell_numbering_map_.end() )
      return it->second;
    else
      DUNE_THROW(Dune::RangeError, "no number for entity");
  }

  /** use 'cp_num_manager.get_number_of_cell_problem( it )'
   * \attention 'get_number_of_cell_problem( it )' is NOT equal to 'get_number_of_cell_problem( it , 0 )'!
   **/
  inline int get_number_of_cell_problem(const EntityPointerType& ent) const {
    auto it = cell_numbering_map_NL_.find(ent);
    if (it != cell_numbering_map_NL_.end() )
      return it->second;
    else
      DUNE_THROW(Dune::RangeError, "no number for entity");
  }
};
} //namespace Dune {
#endif // DUNEMS_HMM_CELL_NUMBERING_HH
