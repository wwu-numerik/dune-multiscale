#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>
#include <memory>
#include <algorithm>

#include "msfem_grid_specifier.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

MacroMicroGridSpecifier::MacroMicroGridSpecifier(DiscreteFunctionSpaceType& coarse_scale_space)
  : coarse_scale_space_(coarse_scale_space)
  , coarse_level_fine_level_difference_(std::numeric_limits<int>::max())
  , number_of_level_host_entities_(coarse_scale_space.gridPart().grid().size(0 /*codim*/))
  , coarseGridIsSimplex_(coarse_scale_space.gridPart().grid().leafIndexSet().geomTypes(0).size() == 1 &&
                         coarse_scale_space.gridPart().grid().leafIndexSet().geomTypes(0)[0].isSimplex()) {
  boundary_nodes_identified_ = false;
}
// get number of coarse grid entities
std::size_t MacroMicroGridSpecifier::getNumOfCoarseEntities() const { return number_of_level_host_entities_; }

//! Get the difference between coarse and fine level
int MacroMicroGridSpecifier::getLevelDifference() const { return coarse_level_fine_level_difference_; }

//! the coarse space
const MacroMicroGridSpecifier::DiscreteFunctionSpaceType& MacroMicroGridSpecifier::coarseSpace() const {
  return coarse_scale_space_;
}
//! the coarse space
MacroMicroGridSpecifier::DiscreteFunctionSpaceType& MacroMicroGridSpecifier::coarseSpace() {
  return coarse_scale_space_;
}

#ifdef ENBABLE_LOD_ONLY_CODE
void MacroMicroGridSpecifier::identify_coarse_boundary_nodes() {
  is_boundary_node_.resize(coarse_scale_space_.size());

  number_of_coarse_boundary_nodes_ = 0;

  const auto endit = coarse_scale_space_.end();
  for (auto it = coarse_scale_space_.begin(); it != endit; ++it) {

    std::vector<std::size_t> indices;
    coarse_scale_space_.mapper().map(*it, indices);

    auto intersection_it = coarse_scale_space_.gridPart().ibegin(*it);
    const auto endiit = coarse_scale_space_.gridPart().iend(*it);
    for (; intersection_it != endiit; ++intersection_it) {

      if (!intersection_it->boundary())
        continue;

      const auto& lagrangePointSet = coarse_scale_space_.lagrangePointSet(*it);

      const int face = (*intersection_it).indexInInside();
      for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face))
        is_boundary_node_[indices[lp]] = true;
    }
  }

  for (size_t i = 0; i < is_boundary_node_.size(); ++i) {
    if (is_boundary_node_[i])
      number_of_coarse_boundary_nodes_ += 1;
  }

  boundary_nodes_identified_ = true;
}

void MacroMicroGridSpecifier::identify_coarse_dirichlet_nodes() {
  is_dirichlet_node_.resize(coarse_scale_space_.size());

  number_of_coarse_dirichlet_nodes_ = 0;

  const auto endit = coarse_scale_space_.end();
  for (auto it = coarse_scale_space_.begin(); it != endit; ++it) {

    std::vector<std::size_t> indices;
    coarse_scale_space_.mapper().map(*it, indices);

    auto intersection_it = coarse_scale_space_.gridPart().ibegin(*it);
    const auto endiit = coarse_scale_space_.gridPart().iend(*it);
    for (; intersection_it != endiit; ++intersection_it) {

      if (!intersection_it->boundary())
        continue;

      if (intersection_it->boundary() && (intersection_it->boundaryId() != 1))
        continue;

      const auto& lagrangePointSet = coarse_scale_space_.lagrangePointSet(*it);

      const int face = (*intersection_it).indexInInside();
      for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face))
        is_dirichlet_node_[indices[lp]] = true;
    }
  }

  for (size_t i = 0; i < is_dirichlet_node_.size(); ++i) {
    if (is_dirichlet_node_[i])
      number_of_coarse_dirichlet_nodes_ += 1;
  }

  dirichlet_nodes_identified_ = true;
}

std::size_t MacroMicroGridSpecifier::get_number_of_coarse_boundary_nodes() const {
  assert(boundary_nodes_identified_);
  return number_of_coarse_boundary_nodes_;
}

std::size_t MacroMicroGridSpecifier::get_number_of_coarse_dirichlet_nodes() const {
  assert(dirichlet_nodes_identified_);
  return number_of_coarse_dirichlet_nodes_;
}

bool MacroMicroGridSpecifier::is_coarse_boundary_node(std::size_t global_index) const {
  assert(boundary_nodes_identified_);
  return is_boundary_node_[global_index];
}

bool MacroMicroGridSpecifier::is_coarse_dirichlet_node(std::size_t global_index) const {
  assert(dirichlet_nodes_identified_);
  return is_dirichlet_node_[global_index];
}
#endif // ENBABLE_LOD_ONLY_CODE

bool MacroMicroGridSpecifier::simplexCoarseGrid() const { return coarseGridIsSimplex_; }

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
