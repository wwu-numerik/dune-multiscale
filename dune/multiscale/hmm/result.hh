// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_HMM_RESULT_HH
#define DUNE_MULTISCALE_HMM_RESULT_HH

#include <dune/multiscale/common/traits.hh>
#include <vector>

namespace Dune {
namespace Multiscale {
namespace HMM {

//! data container for all the different types of numerical errors in a HMM simulation
struct HMMResult {
    typename CommonTraits::RangeType estimated_source_error;
    typename CommonTraits::RangeType estimated_approximation_error;
    typename CommonTraits::RangeType estimated_residual_error;
    typename CommonTraits::RangeType estimated_residual_error_micro_jumps;
    typename CommonTraits::RangeType estimated_residual_error_macro_jumps;
    typename CommonTraits::RangeType estimated_tfr_error;
    std::vector<typename CommonTraits::RangeType> local_error_indicator;
    typename CommonTraits::RangeType minimal_loc_indicator;
    typename CommonTraits::RangeType maximal_loc_indicator;
    typename CommonTraits::RangeType average_loc_indicator;
    typename CommonTraits::RangeType estimated_error;
    double max_variation;
    double min_variation;
    HMMResult(std::size_t codim0_count)
        : estimated_source_error(0.0)
        , estimated_approximation_error(0.0)
        , estimated_residual_error(0.0)
        , estimated_residual_error_micro_jumps(0.0)
        , estimated_residual_error_macro_jumps(0.0)
        , estimated_tfr_error(0.0)
        , local_error_indicator(codim0_count, 0.0)
        , minimal_loc_indicator(10000.0)
        , maximal_loc_indicator(0.0)
        , average_loc_indicator(0.0)
        , estimated_error(0.0)
        , max_variation(0.0)
        , min_variation(0.0)
    {}
};

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {

#endif // RESULT_HH
