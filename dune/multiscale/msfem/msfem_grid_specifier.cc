#include "msfem_grid_specifier.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {


MacroMicroGridSpecifier::MacroMicroGridSpecifier(DiscreteFunctionSpaceType& coarse_scale_space,
                                                         DiscreteFunctionSpaceType& fine_scale_space)
    : coarse_scale_space_(coarse_scale_space)
    , fine_scale_space_(fine_scale_space)
    , coarse_level_fine_level_difference_( fine_scale_space.gridPart().grid().maxLevel()
                                           - coarse_scale_space.gridPart().grid().maxLevel() )
    , number_of_level_host_entities_( coarse_scale_space.gridPart().grid().size(0 /*codim*/) )
    , number_of_layers(number_of_level_host_entities_, 0)
{ boundary_nodes_identified_ = false; }

// get number of coarse grid entities
int MacroMicroGridSpecifier::getNumOfCoarseEntities() const {
    return number_of_level_host_entities_;
}

/** Set the number of overlay layers for each coarse element.
  *
  * @param[in] i The number of the coarse element
  * @param[in] number_of_layers_for_entity The number of overlay layers that shall be provided for the given coarse
  * element
  */
void MacroMicroGridSpecifier::setNoOfLayers(int i, int number_of_layers_for_entity) {
    if (i < number_of_level_host_entities_) {
        number_of_layers[i] = number_of_layers_for_entity;
    } else {
        DUNE_THROW(Dune::InvalidStateException,"Error. Assertion (i < number_of_level_host_entities_) not fulfilled.");
    }
} // setNoOfLayers

/** Get the number of overlay layers for a given coarse element.
  *
  * @param[in] i The number of the coarse element
  * @return Returns the number of overlay layers for the given coarse element.
  */
int MacroMicroGridSpecifier::getNoOfLayers(int i) const {
    if (i < number_of_level_host_entities_) {
        return number_of_layers[i];
    } else {
        DUNE_THROW(Dune::InvalidStateException,"Error. Assertion (i < number_of_level_host_entities_) not fulfilled.");
    }
    return 0;
} // getNoOfLayers

//! Get the difference between coarse and fine level
int MacroMicroGridSpecifier::getLevelDifference() const {
    return coarse_level_fine_level_difference_;
}

//! the coarse space
const MacroMicroGridSpecifier::DiscreteFunctionSpaceType& MacroMicroGridSpecifier::coarseSpace() const {
    return coarse_scale_space_;
}
//! the coarse space
MacroMicroGridSpecifier::DiscreteFunctionSpaceType& MacroMicroGridSpecifier::coarseSpace() {
    return coarse_scale_space_;
}

//! the fine space
const MacroMicroGridSpecifier::DiscreteFunctionSpaceType& MacroMicroGridSpecifier::fineSpace() const {
    return fine_scale_space_;
}
//! the fine space
MacroMicroGridSpecifier::DiscreteFunctionSpaceType& MacroMicroGridSpecifier::fineSpace() {
    return fine_scale_space_;
}

void MacroMicroGridSpecifier::setOversamplingStrategy( int oversampling_strategy )
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

int MacroMicroGridSpecifier::getOversamplingStrategy() const
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

void MacroMicroGridSpecifier::initialize_local_error_manager() {
    //!TODO previous imp would append number_of_level_host_entities_ many zeroes every time it was called
    loc_coarse_residual_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_projection_error_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_coarse_grid_jumps_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_conservative_flux_jumps_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_approximation_error_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
    loc_fine_grid_jumps_ = RangeTypeVector(number_of_level_host_entities_, 0.0);
} // initialize_local_error_manager

void MacroMicroGridSpecifier::set_loc_coarse_residual(int index, const RangeType& loc_coarse_residual) {
    loc_coarse_residual_[index] = loc_coarse_residual;
}

void MacroMicroGridSpecifier::set_loc_coarse_grid_jumps(int index, const RangeType& loc_coarse_grid_jumps) {
    loc_coarse_grid_jumps_[index] = loc_coarse_grid_jumps;
}

void MacroMicroGridSpecifier::set_loc_projection_error(int index, const RangeType& loc_projection_error) {
    loc_projection_error_[index] = loc_projection_error;
}

void MacroMicroGridSpecifier::set_loc_conservative_flux_jumps(int index, const RangeType& loc_conservative_flux_jumps) {
    loc_conservative_flux_jumps_[index] = loc_conservative_flux_jumps;
}

void MacroMicroGridSpecifier::set_loc_approximation_error(int index, const RangeType& loc_approximation_error) {
    loc_approximation_error_[index] = loc_approximation_error;
}

void MacroMicroGridSpecifier::set_loc_fine_grid_jumps(int index, const MacroMicroGridSpecifier::RangeType& loc_fine_grid_jumps) {
    loc_fine_grid_jumps_[index] = loc_fine_grid_jumps;
}

MacroMicroGridSpecifier::RangeType MacroMicroGridSpecifier::get_loc_coarse_residual(int index) const {
    if (loc_coarse_residual_.size() == 0)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_coarse_residual_[index];
} // get_loc_coarse_residual

MacroMicroGridSpecifier::RangeType MacroMicroGridSpecifier::get_loc_coarse_grid_jumps(int index) const {
    if (loc_coarse_grid_jumps_.size() == 0)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_coarse_grid_jumps_[index];
} // get_loc_coarse_grid_jumps

MacroMicroGridSpecifier::RangeType MacroMicroGridSpecifier::get_loc_projection_error(int index) const {
    if (loc_projection_error_.size() == 0)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_projection_error_[index];
} // get_loc_projection_error

MacroMicroGridSpecifier::RangeType MacroMicroGridSpecifier::get_loc_conservative_flux_jumps(int index) const {
    if (loc_conservative_flux_jumps_.size() == 0)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_conservative_flux_jumps_[index];
} // get_loc_conservative_flux_jumps

MacroMicroGridSpecifier::RangeType MacroMicroGridSpecifier::get_loc_approximation_error(int index) const {
    if (loc_approximation_error_.size() == 0)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_approximation_error_[index];
} // get_loc_approximation_error

MacroMicroGridSpecifier::RangeType MacroMicroGridSpecifier::get_loc_fine_grid_jumps(int index) const {
    if (loc_fine_grid_jumps_.size() == 0)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Error! Use: initialize_local_error_manager()-method for the grid specifier first!");
    }
    return loc_fine_grid_jumps_[index];
} // get_loc_fine_grid_jumps

void MacroMicroGridSpecifier::identify_coarse_boundary_nodes()
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
                    = lagrangePointSet.beginSubEntity< faceCodim >(face);
            const auto faceEndIterator
                    = lagrangePointSet.endSubEntity< faceCodim >(face);
            for ( ; faceIterator != faceEndIterator; ++faceIterator)
                is_boundary_node_[coarse_scale_space_.mapper().mapToGlobal(*it, *faceIterator )] = true;

        }

    }

    for ( size_t i = 0; i < is_boundary_node_.size(); ++i )
    {
        if ( is_boundary_node_[i] )
            number_of_coarse_boundary_nodes_ += 1;
    }

    boundary_nodes_identified_ = true;

}

int MacroMicroGridSpecifier::get_number_of_coarse_boundary_nodes() const
{
    assert( boundary_nodes_identified_ );
    return number_of_coarse_boundary_nodes_;
}

bool MacroMicroGridSpecifier::is_coarse_boundary_node( int global_index ) const
{
    assert( boundary_nodes_identified_ );
    return is_boundary_node_[global_index];
}

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

