#ifndef DUNE_MULTISCALE_ERROR_CONTAINER_HH
#define DUNE_MULTISCALE_ERROR_CONTAINER_HH

#include <dune/multiscale/common/traits.hh>
#include <vector>

namespace Dune {
namespace Multiscale {

struct ErrorContainer
{
    ErrorContainer(std::size_t max_loop_number)
        : loc_coarse_residual_(max_loop_number)
        , loc_coarse_grid_jumps_(max_loop_number)
        , loc_projection_error_(max_loop_number)
        , loc_conservative_flux_jumps_(max_loop_number)
        , loc_approximation_error_(max_loop_number)
        , loc_fine_grid_jumps_(max_loop_number)
        , total_coarse_residual_(max_loop_number)
        , total_projection_error_(max_loop_number)
        , total_coarse_grid_jumps_(max_loop_number)
        , total_conservative_flux_jumps_(max_loop_number)
        , total_approximation_error_(max_loop_number)
        , total_fine_grid_jumps_(max_loop_number)
        , total_estimated_H1_error_(max_loop_number)
        , locals({{ &loc_coarse_residual_, &loc_coarse_grid_jumps_,
                 &loc_projection_error_, &loc_conservative_flux_jumps_,
                 &loc_approximation_error_, &loc_fine_grid_jumps_}})
        , totals({{&total_coarse_residual_, &total_projection_error_,
                 &total_coarse_grid_jumps_, &total_conservative_flux_jumps_,
                 &total_approximation_error_, &total_fine_grid_jumps_ }})
    {}

    // local coarse residual, i.e. H ||f||_{L^2(T)}
    CommonTraits::RangeVectorVector loc_coarse_residual_;
    // local coarse grid jumps (contribute to the total coarse residual)
    CommonTraits::RangeVectorVector loc_coarse_grid_jumps_;
    // local projection error (we project to get a globaly continous approximation)
    CommonTraits::RangeVectorVector loc_projection_error_;
    // local jump in the conservative flux
    CommonTraits::RangeVectorVector loc_conservative_flux_jumps_;
    // local approximation error
    CommonTraits::RangeVectorVector loc_approximation_error_;
    // local sum over the fine grid jumps (for a fixed subgrid that cooresponds with a coarse entity T)
    CommonTraits::RangeVectorVector loc_fine_grid_jumps_;

    CommonTraits::RangeVector total_coarse_residual_;
    CommonTraits::RangeVector total_projection_error_;
    CommonTraits::RangeVector total_coarse_grid_jumps_;
    CommonTraits::RangeVector total_conservative_flux_jumps_;
    CommonTraits::RangeVector total_approximation_error_;
    CommonTraits::RangeVector total_fine_grid_jumps_;
    CommonTraits::RangeVector total_estimated_H1_error_;

    std::vector<CommonTraits::RangeVectorVector*> locals;
    std::vector<CommonTraits::RangeVector*> totals;
};

} //namespace Multiscale {
} //namespace Dune {


#endif // DUNE_MULTISCALE_ERROR_CONTAINER_HH
