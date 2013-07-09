// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS__TOOLS_MISC_HH
#define DUNE_MS__TOOLS_MISC_HH

#include <dune/stuff/common/string.hh>
#include <dune/grid/common/entitypointer.hh>
#include <dune/grid/common/indexidset.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

/** ConstGridImp is supposed to be identical to const GridImp, alas different grid managers
 *  and different type deduction between gcc/clang make two different types necessary here
 */
template < class GridImp, class ConstGridImp, class IteratorImp, class IndexSetImp, class IndexTypeImp>
EntityPointer< GridImp, IteratorImp > make_father(const IndexSet< ConstGridImp, IndexSetImp, IndexTypeImp >& index_set,
                EntityPointer< GridImp, IteratorImp > entity )
{
  bool father_found = false;
  auto test_entity = entity;
  while (!father_found)
  {
    if (index_set.contains(*test_entity))
    {
      entity = test_entity;
    }

    if (!test_entity->hasFather())
    {
      father_found = true;
    } else {
      test_entity = test_entity->father();
    }
  }
  return entity;
}


template < class GridImp, class ConstGridImp, class IteratorImp, class IndexSetImp, class IndexTypeImp >
EntityPointer< GridImp, IteratorImp > make_father(const IndexSet< ConstGridImp, IndexSetImp, IndexTypeImp >& index_set,
                 EntityPointer< GridImp, IteratorImp > entity,
                 int level_difference)
{
    for (int lev = 0; lev < level_difference; ++lev)
          entity = entity->father();
    return make_father(index_set, entity);
}

template < int cd, int dim,
           template<int,int,class> class EntityImp,
           template<int,int,class> class EntityImpOther,
           class GridImp,
           class GridImpOther>
bool entities_identical(const Entity< cd, dim, GridImp, EntityImp >& entity,
                        const Entity< cd, dim, GridImpOther, EntityImpOther >& other)
{
    const int number_of_nodes = entity.template count< GridImp::dimension >();
    for (int k = 0; k < number_of_nodes; ++k )
    {
      if ( !( entity.geometry().corner(k) == other.geometry().corner(k) ) )
      {
        return false;
      }
    }
    return true;
}

} // namespace Grid
} // namespace Stuff
} // namespace Dune


#endif // DUNE_MS__TOOLS_MISC_HH
