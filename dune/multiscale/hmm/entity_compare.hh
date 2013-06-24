#ifndef DUNEMS_ENTITY_COMPARE_HH
#define DUNEMS_ENTITY_COMPARE_HH

#include <utility>
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
    const auto& geometry_right = ( *(right_entity_pair.first) ).geometry();
    const Quadrature quadrature_right( ( *(right_entity_pair.first) ), 0 );

    // local barycenter (with respect to entity)
    const auto& local_point_right = quadrature_right.point(0);
    const DomainType barycenter_right_entity = geometry_right.global(local_point_right);

    // ------ left element
    const auto& geometry_left = ( *(left_entity_pair.first) ).geometry();
    const Quadrature quadrature_left( ( *(left_entity_pair.first) ), 0 );

    // local barycenter (with respect to entity)
    const auto& local_point_left = quadrature_left.point(0);
    const DomainType barycenter_left_entity = geometry_left.global(local_point_left);

    int current_axis = GridPartType::GridType::dimension - 1;

    while (current_axis >= 0)
    {
      if (barycenter_left_entity[current_axis] < barycenter_right_entity[current_axis])
      { return true; } else if (barycenter_left_entity[current_axis] > barycenter_right_entity[current_axis])
      { return false; }

      current_axis -= 1;
    }

    return left_entity_pair.second < right_entity_pair.second;
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
    const Quadrature quadrature_right(*right_entity, 0);

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_right = quadrature_right.point(0);
    const DomainType barycenter_right_entity = geometry_right.global(local_point_right);

    // ------ left element
    const typename EntityPointerType::Entity::Geometry& geometry_left = (*left_entity).geometry();
    const Quadrature quadrature_left(*left_entity, 0);

    // local barycenter (with respect to entity)
    const typename Quadrature::CoordinateType& local_point_left = quadrature_left.point(0);
    const DomainType barycenter_left_entity = geometry_left.global(local_point_left);

    int current_axis = GridPartType::GridType::dimension - 1;
    while (current_axis >= 0)
    {
      if (barycenter_left_entity[current_axis] < barycenter_right_entity[current_axis]) {
        return true;
      } else if (barycenter_left_entity[current_axis] > barycenter_right_entity[current_axis]) {
        return false;
      }
      current_axis--;
    }

    return false;
  } // ()
};

} // namespace Dune {

#endif // ENTITY_COMPARE_HH
