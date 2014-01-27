// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MSFEM_GRID_SPECIFIER_HH
#define MSFEM_GRID_SPECIFIER_HH

#include <dune/multiscale/common/traits.hh>
#include <cstddef>
#include <vector>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class MacroMicroGridSpecifier {
  //! \todo DiscreteFunctionSpaceType should be replaced be something like "CoarseDiscreteFunctionSpace" and
  // "FineDiscreteFunctionSpace"
  typedef typename CommonTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  static const int faceCodim = 1;

public:
  MacroMicroGridSpecifier(DiscreteFunctionSpaceType& coarse_scale_space);

  //! Get the difference between coarse and fine level
  int getLevelDifference() const;

  const DiscreteFunctionSpaceType& coarseSpace() const;
  DiscreteFunctionSpaceType& coarseSpace();

#if 0 // LOD Only
  void identify_coarse_boundary_nodes();
  void identify_coarse_dirichlet_nodes();

  std::size_t get_number_of_coarse_boundary_nodes() const;
  std::size_t get_number_of_coarse_dirichlet_nodes() const;
private:

  // number of coarse grid boundary nodes
  std::size_t number_of_coarse_boundary_nodes_;

  // number of coarse grid boundary nodes
  std::size_t number_of_coarse_dirichlet_nodes_;

  // have the Dirichlet boundary nodes been identified?
  bool dirichlet_nodes_identified_;
  // have the boundary nodes been identified?
  bool boundary_nodes_identified_;

  // is a given coarse node a boundary node of the coarse grid? true/false
  std::vector<bool> is_boundary_node_;

  // is a given coarse node a Dirichlet boundary node of the coarse grid? true/false
  std::vector<bool> is_dirichlet_node_;

public:
  bool is_coarse_boundary_node(std::size_t global_index) const;
  bool is_coarse_dirichlet_node(std::size_t global_index) const;
#endif //0 // LOD Only

  bool simplexCoarseGrid() const;

private:
  DiscreteFunctionSpaceType& coarse_scale_space_;

  // level difference between coarse grid level and fine grid level
  const int coarse_level_fine_level_difference_;



  const bool coarseGridIsSimplex_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif
