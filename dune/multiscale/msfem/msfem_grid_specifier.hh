// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MSFEM_GRID_SPECIFIER_HH
#define MSFEM_GRID_SPECIFIER_HH

#include <dune/multiscale/common/traits.hh>

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
  MacroMicroGridSpecifier(DiscreteFunctionSpaceType& coarse_scale_space, DiscreteFunctionSpaceType& fine_scale_space);

  // get number of coarse grid entities
  std::size_t getNumOfCoarseEntities() const;

  /** Set the number of overlay layers for each coarse element.
  *
  * @param[in] i The number of the coarse element
  * @param[in] number_of_layers_for_entity The number of overlay layers that shall be provided for the given coarse
  * element
  */
  void setNoOfLayers(std::size_t i, std::size_t number_of_layers_for_entity);

  /** Get the number of overlay layers for a given coarse element.
  *
  * @param[in] i The number of the coarse element
  * @return Returns the number of overlay layers for the given coarse element.
  */
  std::size_t getNoOfLayers(std::size_t i) const;

  /** Get the maximum number of overlay layers for the whole coarse grid.
  * @return Returns the maximum number of overlay layers for the whole coarse grid.
  */
  std::size_t maxNumberOverlayLayers() const;

  //! Get the difference between coarse and fine level
  int getLevelDifference() const;

  const DiscreteFunctionSpaceType& coarseSpace() const;
  DiscreteFunctionSpaceType& coarseSpace();
  const DiscreteFunctionSpaceType& fineSpace() const;
  DiscreteFunctionSpaceType& fineSpace();

  void setOversamplingStrategy(int oversampling_strategy);
  int getOversamplingStrategy() const;

  void initialize_local_error_manager();

  void set_loc_coarse_residual(std::size_t index, const RangeType& loc_coarse_residual);
  void set_loc_coarse_grid_jumps(std::size_t index, const RangeType& loc_coarse_grid_jumps);
  void set_loc_projection_error(std::size_t index, const RangeType& loc_projection_error);
  void set_loc_conservative_flux_jumps(std::size_t index, const RangeType& loc_conservative_flux_jumps);
  void set_loc_approximation_error(std::size_t index, const RangeType& loc_approximation_error);
  void set_loc_fine_grid_jumps(std::size_t index, const RangeType& loc_fine_grid_jumps);
  RangeType get_loc_coarse_residual(std::size_t index) const;
  RangeType get_loc_coarse_grid_jumps(std::size_t index) const;
  RangeType get_loc_projection_error(std::size_t index) const;
  RangeType get_loc_conservative_flux_jumps(std::size_t index) const;
  RangeType get_loc_approximation_error(std::size_t index) const;
  RangeType get_loc_fine_grid_jumps(std::size_t index) const;

  void identify_coarse_boundary_nodes();
  void identify_coarse_dirichlet_nodes();

  std::size_t get_number_of_coarse_boundary_nodes() const;
  std::size_t get_number_of_coarse_dirichlet_nodes() const;

  bool is_coarse_boundary_node(std::size_t global_index) const;
  bool is_coarse_dirichlet_node(std::size_t global_index) const;

  bool simplexCoarseGrid() const;

private:
  DiscreteFunctionSpaceType& coarse_scale_space_;
  DiscreteFunctionSpaceType& fine_scale_space_;

  // level difference between coarse grid level and fine grid level
  const int coarse_level_fine_level_difference_;

  // number of coarse grid entities
  const std::size_t number_of_level_host_entities_;

  // oversampling strategy - 1, 2 or 3. (1 and 2 for MsFEM and 3 for Rigorous MsFEM)
  int oversampling_strategy_;

  // layers for each coarse grid entity
  std::vector<std::size_t> number_of_layers;

  // have the boundary nodes been identified?
  bool boundary_nodes_identified_;

  // have the Dirichlet boundary nodes been identified?
  bool dirichlet_nodes_identified_;

  // is a given coarse node a boundary node of the coarse grid? true/false
  std::vector<bool> is_boundary_node_;

  // is a given coarse node a Dirichlet boundary node of the coarse grid? true/false
  std::vector<bool> is_dirichlet_node_;

  // number of coarse grid boundary nodes
  std::size_t number_of_coarse_boundary_nodes_;

  // number of coarse grid boundary nodes
  std::size_t number_of_coarse_dirichlet_nodes_;

  // ----- local error indicators (for each coarse grid element T) -------------

  // local coarse residual, i.e. H ||f||_{L^2(T)}
  typedef std::vector<RangeType> RangeTypeVector;
  RangeTypeVector loc_coarse_residual_;

  // local coarse grid jumps (contribute to the total coarse residual)
  RangeTypeVector loc_coarse_grid_jumps_;

  // local projection error (we project to get a globaly continous approximation)
  RangeTypeVector loc_projection_error_;

  // local jump in the conservative flux
  RangeTypeVector loc_conservative_flux_jumps_;

  // local approximation error
  RangeTypeVector loc_approximation_error_;

  // local sum over the fine grid jumps (for a fixed subgrid that cooresponds with a coarse entity T)
  RangeTypeVector loc_fine_grid_jumps_;

  const bool coarseGridIsSimplex_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif
