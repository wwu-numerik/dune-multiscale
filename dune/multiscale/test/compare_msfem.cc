
#include <dune/stuff/test/test_common.hh>
#include <gtest.h>

#include <string>
#include <array>
#include <initializer_list>
#include <vector>
#include <dune/stuff/common/fixed_map.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/ranges.hh>

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/multiscale/problems/elliptic/selector.hh>

using namespace Dune::Stuff::Common;

void set_param();

TEST(MSFEM, All) {
    using namespace Dune::Multiscale;
    using namespace Dune::Multiscale::MsFEM;

    int coarse_grid_level_ = 3;
    int number_of_layers_ = 4;
    int total_refinement_level_ = 5;
    const std::string macroGridName("compare_msfem.dgf") ;
    set_param();

    //! ---------------------- local error indicators --------------------------------
    // ----- local error indicators (for each coarse grid element T) -------------
    const int max_loop_number = 10;
    // local coarse residual, i.e. H ||f||_{L^2(T)}
    CommonTraits::RangeVectorVector loc_coarse_residual_(max_loop_number);
    // local coarse grid jumps (contribute to the total coarse residual)
    CommonTraits::RangeVectorVector loc_coarse_grid_jumps_(max_loop_number);
    // local projection error (we project to get a globaly continous approximation)
    CommonTraits::RangeVectorVector loc_projection_error_(max_loop_number);
    // local jump in the conservative flux
    CommonTraits::RangeVectorVector loc_conservative_flux_jumps_(max_loop_number);
    // local approximation error
    CommonTraits::RangeVectorVector loc_approximation_error_(max_loop_number);
    // local sum over the fine grid jumps (for a fixed subgrid that cooresponds with a coarse entity T)
    CommonTraits::RangeVectorVector loc_fine_grid_jumps_(max_loop_number);

    CommonTraits::RangeVector total_coarse_residual_(max_loop_number);
    CommonTraits::RangeVector total_projection_error_(max_loop_number);
    CommonTraits::RangeVector total_coarse_grid_jumps_(max_loop_number);
    CommonTraits::RangeVector total_conservative_flux_jumps_(max_loop_number);
    CommonTraits::RangeVector total_approximation_error_(max_loop_number);
    CommonTraits::RangeVector total_fine_grid_jumps_(max_loop_number);
    CommonTraits::RangeVector total_estimated_H1_error_(max_loop_number);

    //! TODO put these into something like a named tuple/class
    std::vector<CommonTraits::RangeVectorVector*> locals = {{ &loc_coarse_residual_, &loc_coarse_grid_jumps_,
                                                             &loc_projection_error_, &loc_conservative_flux_jumps_,
                                                             &loc_approximation_error_, &loc_fine_grid_jumps_}};
    std::vector<CommonTraits::RangeVector*> totals = {{&total_coarse_residual_, &total_projection_error_,
                                                      &total_coarse_grid_jumps_, &total_conservative_flux_jumps_,
                                                      &total_approximation_error_, &total_fine_grid_jumps_ }};

    unsigned int loop_number = 0;
    while (algorithm(macroGridName, loop_number++, total_refinement_level_, coarse_grid_level_,
                     number_of_layers_, locals, totals, total_estimated_H1_error_))
    {}

    for(const auto total : totals)
    {
        for(auto error : *total)
        {
            // fails iff 1e-2 <= error
            EXPECT_GT(double(1e-2),error);
        }
    }
}

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}


void set_param()
{

}
