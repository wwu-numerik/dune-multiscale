#include <config.h>
#include "heterogenous.hh"

#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localsolution_proxy.hh>

void Dune::Multiscale::MsFEMProjection::project(const Dune::Multiscale::LocalsolutionProxy &source, Dune::Multiscale::CommonTraits::DiscreteFunctionType &target, Dune::Multiscale::LocalGridSearch &search) {
  constexpr size_t target_dimRange = CommonTraits::dimRange;

  const auto& space = target.space();

  preprocess(target);
  const auto interior = space.grid_view().grid().template leafGridView<Interior_Partition>();
  for (const auto& target_entity : DSC::entityRange(interior)) {
    auto target_local_function = target.local_discrete_function(target_entity);
    const auto global_quads = global_evaluation_points(space, target_entity);
    const auto evaluation_entity_ptrs = search(global_quads);
    assert(evaluation_entity_ptrs.size() >= global_quads.size());

    size_t k = 0;
    typename CommonTraits::DiscreteFunctionType::RangeType source_value;
    for (size_t qP = 0; qP < global_quads.size(); ++qP) {
      const auto& source_entity_unique_ptr = evaluation_entity_ptrs[qP];
      if (source_entity_unique_ptr) {
        const auto& source_entity_ptr = (*source_entity_unique_ptr);
        const auto& source_geometry = source_entity_ptr->geometry();
        const auto& global_point = global_quads[qP];
        const auto& source_local_point = source_geometry.local(global_point);
        const auto& ent = *source_entity_ptr;
        const auto& source_local_function = source.local_function(ent);
        source_value = source_local_function->evaluate(source_local_point);
        for (size_t i = 0; i < target_dimRange; ++i, ++k) {
          target_local_function->vector().add(k, source_value[i]);
        }
      } else {
        DUNE_THROW(InvalidStateException, "Did not find the local lagrange point in the source mesh!");
      }
    }
  }

  postprocess(target);
}

void Dune::Multiscale::MsFEMProjection::preprocess(Dune::Multiscale::CommonTraits::DiscreteFunctionType &func) {
  // set all DoFs to zero
  func.vector() *= 0;
}

void Dune::Multiscale::MsFEMProjection::postprocess(Dune::Multiscale::CommonTraits::DiscreteFunctionType &func) {
  // compute node to entity relations
  std::vector<int> nodeToEntity(func.space().grid_view().grid().size(CommonTraits::world_dim), 0);
  identifySharedNodes(func.space().grid_view(), nodeToEntity);

  auto factorsIt = nodeToEntity.begin();
  for (auto& dit : func.vector()) {
    assert(factorsIt != nodeToEntity.end());
    assert(*factorsIt > 0);
    dit /= *factorsIt;
    ++factorsIt;
  }
  return;
}

void Dune::Multiscale::MsFEMProjection::identifySharedNodes(const Dune::Multiscale::CommonTraits::GridViewType &gridPart, std::vector<int> &map) {
  const auto& indexSet = gridPart.indexSet();

  for (auto& entity : DSC::entityRange(gridPart.grid().leafGridView())) {
    const auto number_of_nodes_in_entity = entity.template count<CommonTraits::world_dim>();
    for (auto i : DSC::valueRange(number_of_nodes_in_entity)) {
      const auto node = entity.template subEntity<CommonTraits::world_dim>(i);
      const auto global_index_node = indexSet.index(*node);

      // make sure we don't access non-existing elements
      assert(map.size() > global_index_node);
      ++map[global_index_node];
    }
  }
}
