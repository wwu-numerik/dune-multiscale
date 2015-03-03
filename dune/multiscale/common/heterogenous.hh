#ifndef DUNE_MULTISCALE_COMMON_HETEROGENOUS_HH
#define DUNE_MULTISCALE_COMMON_HETEROGENOUS_HH

#include <dune/common/deprecated.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/io/file/vtk/function.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/grid/search.hh>

#include <dune/gdt/spaces/cg/interface.hh>

namespace Dune {
namespace Multiscale {

template< class ImpTraits, size_t domainDim, size_t rangeDim >
std::vector<typename GDT::Spaces::CGInterface< ImpTraits, domainDim, rangeDim, 1 >::DomainType>
global_evaluation_points(const GDT::Spaces::CGInterface< ImpTraits, domainDim,rangeDim, 1 >& space,
         const typename GDT::Spaces::CGInterface< ImpTraits, domainDim, rangeDim, 1 >::EntityType& target_entity)
{
  const auto& target_lagrangepoint_set = space.lagrange_points(target_entity);
  const auto& target_geometry = target_entity.geometry();
  const auto quadNop = target_lagrangepoint_set.size();
  std::vector<typename GDT::Spaces::CGInterface< ImpTraits, domainDim, rangeDim, 1 >::DomainType> points(quadNop);
  for(size_t qP = 0; qP < quadNop ; ++qP) {
    points[qP] = target_geometry.global(target_lagrangepoint_set[qP]);
  }
  return points;
}

class MsFEMProjection {
public:
  //! signature for non-default SearchStrategy constructions
  template < class SourceSpaceImp, class TargetSpaceImp, class SourceVectorImp, class TargetVectorImp, class SearchStrategyImp >
  static void project(const GDT::ConstDiscreteFunction< SourceSpaceImp, SourceVectorImp >& source,
                      GDT::DiscreteFunction< TargetSpaceImp, TargetVectorImp >& target,
                      SearchStrategyImp& search)
  {
    constexpr size_t target_dimRange = TargetSpaceImp::dimRange;

    const auto& space =  target.space();

    preprocess(target);
    const auto interior = space.grid_view().grid().template leafGridView<Interior_Partition>();
    for(const auto& target_entity : DSC::entityRange(interior))
    {
      auto target_local_function = target.local_discrete_function(target_entity);
      const auto global_quads = global_evaluation_points(space, target_entity);
      const auto evaluation_entity_ptrs = search(global_quads);
      assert(evaluation_entity_ptrs.size() >= global_quads.size());

      size_t k = 0;
      typename GDT::DiscreteFunction< SourceSpaceImp, SourceVectorImp >::RangeType source_value;
      for(size_t qP = 0; qP < global_quads.size() ; ++qP)
      {
          const auto& source_entity_unique_ptr = evaluation_entity_ptrs[qP];
          if (source_entity_unique_ptr) {
            const auto& source_entity_ptr = (*source_entity_unique_ptr);
            const auto& source_geometry = source_entity_ptr->geometry();
            const auto& global_point = global_quads[qP];
            const auto& source_local_point = source_geometry.local(global_point);
            const auto& ent = *source_entity_ptr;
            const auto& source_local_function = source.local_function(ent);
            source_value = source_local_function->evaluate(source_local_point);
            for(size_t i = 0; i < target_dimRange; ++i, ++k) {
              target_local_function->vector().add(k, source_value[i]);
            }
          }
          else {
            DUNE_THROW(InvalidStateException, "Did not find the local lagrange point in the source mesh!");
          }
        }
      }

    postprocess(target);
  } // ... project(...)

protected:
  template<class TargetSpaceImp, class VectorImp>
  static void preprocess(GDT::DiscreteFunction< TargetSpaceImp, VectorImp >& func) {
    // set all DoFs to zero
    func.vector() *= 0;
  }

  template<class DofType, class SourceType >
  static void setDofValue(DofType& dof, const SourceType& value)
  {
    dof += value;
  }

  template<class TargetSpaceImp, class VectorImp>
  static void postprocess(GDT::DiscreteFunction< TargetSpaceImp, VectorImp >& func) {
    // compute node to entity relations
    constexpr size_t dimension = TargetSpaceImp::GridViewType::Grid::dimension;
    std::vector<int> nodeToEntity(func.space().grid_view().grid().size(dimension), 0);
    identifySharedNodes(func.space().grid_view(), nodeToEntity);

    auto factorsIt = nodeToEntity.begin();
    for (auto& dit : func.vector()) {
      assert(factorsIt!=nodeToEntity.end());
      assert(*factorsIt>0);
      dit /= *factorsIt;
      ++factorsIt;
    }
    return;
  }

  template<class GridPartType, class MapType>
  static void identifySharedNodes(const GridPartType& gridPart, MapType& map) {
    constexpr size_t dimension = GridPartType::Grid::dimension;
    const auto& indexSet = gridPart.indexSet();

    for (auto& entity : DSC::entityRange(gridPart.grid().leafGridView())) {
      const auto number_of_nodes_in_entity = entity.template count<dimension>();
      for (auto i : DSC::valueRange(number_of_nodes_in_entity)) {
        const auto node = entity.template subEntity<dimension>(i);
        const auto global_index_node = indexSet.index(*node);

        // make sure we don't access non-existing elements
        assert(map.size() > global_index_node);
        ++map[global_index_node];
      }
    }
  }

};

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_COMMON_HETEROGENOUS_HH
