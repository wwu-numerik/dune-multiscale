#include <config.h>
#include "heterogenous.hh"

#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>
#include <dune/xt/common/parallel/partitioner.hh>
#include <dune/xt/grid/parallel/partitioning/ranged.hh>

void Dune::Multiscale::MsFEMProjection::project(Dune::Multiscale::LocalsolutionProxy& source,
                                                Dune::Multiscale::CommonTraits::DiscreteFunctionType& target)
{
  constexpr size_t target_dimRange = CommonTraits::dimRange;
  static_assert(target_dimRange == 1, "");

  const auto& space = target.space();

  preprocess(target);
  const auto interior = space.grid_layer().grid().leafGridView();
  typedef std::remove_const<decltype(interior)>::type InteriorType;
  GDT::SystemAssembler<CommonTraits::SpaceType, InteriorType> walker(space, interior);

  Dune::XT::Grid::RangedPartitioning<InteriorType, 0, CommonTraits::InteriorBorderPartition> partitioning(
      interior, XT::Common::threadManager().max_threads());

  const std::function<void(const CommonTraits::EntityType&)> func = [&](const CommonTraits::EntityType& target_entity) {
    auto target_local_function = target.local_discrete_function(target_entity);
    const auto global_quads = global_evaluation_points(space, target_entity);
    auto& search = source.search();
    const auto evaluation_entity_ptrs = search(global_quads);
    assert(evaluation_entity_ptrs.size() >= global_quads.size());

    size_t k = 0;
    typename CommonTraits::DiscreteFunctionType::RangeType source_value;
    for (size_t qP = 0; qP < global_quads.size(); ++qP) {
      const auto& source_entity_unique_ptr = evaluation_entity_ptrs[qP];
      if (source_entity_unique_ptr) {
        const auto source_entity = (*source_entity_unique_ptr);
        const auto& source_geometry = source_entity.geometry();
        const auto& global_point = global_quads[qP];
        const auto& source_local_point = source_geometry.local(global_point);
        const auto& source_local_function = source.local_function(source_entity);
        source_value = source_local_function->evaluate(source_local_point);
        for (size_t i = 0; i < target_dimRange; ++i, ++k) {
          target_local_function->vector().set(k, source_value[i]);
        }
      } else {
        DUNE_THROW(InvalidStateException, "Did not find the local lagrange point in the source mesh!");
      }
    }
  };

  walker.add(func);
  walker.assemble(partitioning);

  postprocess(target);
}

void Dune::Multiscale::MsFEMProjection::preprocess(Dune::Multiscale::CommonTraits::DiscreteFunctionType& func)
{
  // set all DoFs to zero
  func.vector() *= 0;
}

void Dune::Multiscale::MsFEMProjection::postprocess(Dune::Multiscale::CommonTraits::DiscreteFunctionType& func)
{
}

void Dune::Multiscale::MsFEMProjection::identifySharedNodes(
    const Dune::Multiscale::CommonTraits::GridViewType& gridPart, std::vector<int>& map)
{
  const auto& indexSet = gridPart.indexSet();

  for (const auto& entity : Dune::elements(gridPart.grid().leafGridView())) {
    const auto number_of_nodes_in_entity = entity.subEntities(CommonTraits::world_dim);
    for (auto i : Dune::XT::Common::value_range(number_of_nodes_in_entity)) {
      const auto node = entity.template subEntity<CommonTraits::world_dim>(i);
      const auto global_index_node = indexSet.index(node);

      // make sure we don't access non-existing elements
      assert(map.size() > global_index_node);
      ++map[global_index_node];
    }
  }
}
