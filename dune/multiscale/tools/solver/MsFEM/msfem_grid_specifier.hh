#ifndef MSFEM_GRID_SPECIFIER_HH
#define MSFEM_GRID_SPECIFIER_HH

namespace Dune {
namespace Multiscale {
namespace MsFEM {

template< class DiscreteFunctionSpaceType >
class MacroMicroGridSpecifier
{
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  enum { faceCodim = 1 };
  
public:
  MacroMicroGridSpecifier(DiscreteFunctionSpaceType& coarse_scale_space,
                          DiscreteFunctionSpaceType& fine_scale_space)
    : coarse_scale_space_(coarse_scale_space)
    , fine_scale_space_(fine_scale_space)
    , coarse_level_fine_level_difference_( fine_scale_space.gridPart().grid().maxLevel()
                                           - coarse_scale_space.gridPart().grid().maxLevel() )
    , number_of_level_host_entities_( coarse_scale_space.gridPart().grid().size(0 /*codim*/) )
    , number_of_layers(number_of_level_host_entities_, 0)
  { boundary_nodes_identified_ = false; }

  // get number of coarse grid entities
  int getNumOfCoarseEntities() const {
    return number_of_level_host_entities_;
  }

  void setLayer(int i, int number_of_layers_for_entity) {
    if (i < number_of_level_host_entities_)
    { number_of_layers[i] = number_of_layers_for_entity; } else {
      DUNE_THROW(Dune::InvalidStateException,"Error. Assertion (i < number_of_level_host_entities_) not fulfilled.");
    }
  } // setLayer

  int getLayer(int i) const {
    if (i < number_of_level_host_entities_)
    { return number_of_layers[i]; } else {
      DUNE_THROW(Dune::InvalidStateException,"Error. Assertion (i < number_of_level_host_entities_) not fulfilled.");
    }
    return 0;
  } // getLayer

  // difference between coarse and fine level
  int getLevelDifference() const {
    return coarse_level_fine_level_difference_;
  }

  //! the coarse space
  const DiscreteFunctionSpaceType& coarseSpace() const {
    return coarse_scale_space_;
  }
  //! the coarse space
  DiscreteFunctionSpaceType& coarseSpace() {
    return coarse_scale_space_;
  }

  //! the fine space
  const DiscreteFunctionSpaceType& fineSpace() const {
    return fine_scale_space_;
  }
  //! the fine space
  DiscreteFunctionSpaceType& fineSpace() {
    return fine_scale_space_;
  }

  void setOversamplingStrategy( int oversampling_strategy )
  {
    switch ( oversampling_strategy )
    {
      case 1: break;
      case 2: break;
      case 3: break;
      default: DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1 or 2, or you must use 'rigorous MsFEM' (i.e. strategy 3).");
    }
    oversampling_strategy_ = oversampling_strategy;
  }

  int getOversamplingStrategy() const
  {
    switch ( oversampling_strategy_ )
    {
      case 1: break;
      case 2: break;
      case 3: break;
      default: DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1 or 2, or you must use 'rigorous MsFEM' (i.e. strategy 3).");
    }
    return oversampling_strategy_;
  }
  
  void initialize_local_error_manager() {
    //!TODO previous imp would append number_of_level_host_entities_ many zeroes every time it was called
    loc_coarse_residual_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_projection_error_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_coarse_grid_jumps_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_conservative_flux_jumps_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_approximation_error_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_fine_grid_jumps_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
  } // initialize_local_error_manager

  void set_loc_coarse_residual(int index, const RangeType& loc_coarse_residual) {
    loc_coarse_residual_[index] = loc_coarse_residual;
  }

  void set_loc_coarse_grid_jumps(int index, const RangeType& loc_coarse_grid_jumps) {
    loc_coarse_grid_jumps_[index] = loc_coarse_grid_jumps;
  }

  void set_loc_projection_error(int index, const RangeType& loc_projection_error) {
    loc_projection_error_[index] = loc_projection_error;
  }

  void set_loc_conservative_flux_jumps(int index, const RangeType& loc_conservative_flux_jumps) {
    loc_conservative_flux_jumps_[index] = loc_conservative_flux_jumps;
  }

  void set_loc_approximation_error(int index, const RangeType& loc_approximation_error) {
    loc_approximation_error_[index] = loc_approximation_error;
  }

  void set_loc_fine_grid_jumps(int index, const RangeType& loc_fine_grid_jumps) {
    loc_fine_grid_jumps_[index] = loc_fine_grid_jumps;
  }

  RangeType get_loc_coarse_residual(int index) const {
    if (loc_coarse_residual_.size() == 0)
    {
      DUNE_THROW(Dune::InvalidStateException,
                 "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_coarse_residual_[index];
  } // get_loc_coarse_residual

  RangeType get_loc_coarse_grid_jumps(int index) const {
    if (loc_coarse_grid_jumps_.size() == 0)
    {
      DUNE_THROW(Dune::InvalidStateException,
                 "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_coarse_grid_jumps_[index];
  } // get_loc_coarse_grid_jumps

  RangeType get_loc_projection_error(int index) const {
    if (loc_projection_error_.size() == 0)
    {
      DUNE_THROW(Dune::InvalidStateException,
                 "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_projection_error_[index];
  } // get_loc_projection_error

  RangeType get_loc_conservative_flux_jumps(int index) const {
    if (loc_conservative_flux_jumps_.size() == 0)
    {
      DUNE_THROW(Dune::InvalidStateException,
                 "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_conservative_flux_jumps_[index];
  } // get_loc_conservative_flux_jumps

  RangeType get_loc_approximation_error(int index) const {
    if (loc_approximation_error_.size() == 0)
    {
      DUNE_THROW(Dune::InvalidStateException,
                 "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_approximation_error_[index];
  } // get_loc_approximation_error

  RangeType get_loc_fine_grid_jumps(int index) const {
    if (loc_fine_grid_jumps_.size() == 0)
    {
      DUNE_THROW(Dune::InvalidStateException,
                 "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_fine_grid_jumps_[index];
  } // get_loc_fine_grid_jumps
  
  void identify_coarse_boundary_nodes()
  {
    is_boundary_node_.resize( coarse_scale_space_.size() );
    
    number_of_coarse_boundary_nodes_ = 0;
    
    const auto endit = coarse_scale_space_.end();
    for (auto it = coarse_scale_space_.begin(); it != endit; ++it)
    {

      auto intersection_it = coarse_scale_space_.gridPart().ibegin(*it);
      const auto endiit = coarse_scale_space_.gridPart().iend(*it);
      for ( ; intersection_it != endiit; ++intersection_it)
      {

        if ( !(*intersection_it).boundary() )
          continue;

        const auto& lagrangePointSet
          = coarse_scale_space_.lagrangePointSet(*it);

        const int face = (*intersection_it).indexInInside();

        auto faceIterator
          = lagrangePointSet.template beginSubEntity< faceCodim >(face);
        const auto faceEndIterator
          = lagrangePointSet.template endSubEntity< faceCodim >(face);
        for ( ; faceIterator != faceEndIterator; ++faceIterator)
          is_boundary_node_[coarse_scale_space_.mapToGlobal(*it, *faceIterator )] = true;

      }

    }

    for ( size_t i = 0; i < is_boundary_node_.size(); ++i )
    {
      if ( is_boundary_node_[i] )
        number_of_coarse_boundary_nodes_ += 1;
    }
    
    boundary_nodes_identified_ = true;

  }
  
  int get_number_of_coarse_boundary_nodes() const
  {
    assert( boundary_nodes_identified_ );
    return number_of_coarse_boundary_nodes_;
  }
  
  bool is_coarse_boundary_node( int global_index ) const
  {
    assert( boundary_nodes_identified_ );
    return is_boundary_node_[global_index];
  }
  
private:
  DiscreteFunctionSpaceType& coarse_scale_space_;
  DiscreteFunctionSpaceType& fine_scale_space_;

  // level difference bettween coarse grid level and fine grid level
  const int coarse_level_fine_level_difference_;

  // number of coarse grid entities
  const int number_of_level_host_entities_;

  // oversampling strategy - 1, 2 or 3. (1 and 2 for MsFEM and 3 for Rigorous MsFEM)
  int oversampling_strategy_;
  
  // layers for each coarse grid entity
  std::vector< int > number_of_layers;
  
  // have the boundary nodes been identified?
  bool boundary_nodes_identified_;
  
  // is a given coarse node a boundary node of the coarse grid? true/false
  std::vector< bool > is_boundary_node_;

  // number of coarse grid boundary nodes
  int number_of_coarse_boundary_nodes_;

  // ----- local error indicators (for each coarse grid element T) -------------

  // local coarse residual, i.e. H ||f||_{L^2(T)}
  typedef std::vector< RangeType > RangeTypeVector;
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
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif
