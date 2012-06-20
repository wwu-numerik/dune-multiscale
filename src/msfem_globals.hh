#ifndef MSFEM_GLOBALS_HH
#define MSFEM_GLOBALS_HH

// ! ---------------------- important variables ---------------------------------
enum { dimension = GridType::dimension };
// name of the error file in which the data will be saved
std::string filename_;
std::string path_;
int total_refinement_level_;
int coarse_grid_level_;
// only for uniform computations / use the same number of layers for every coarse grid entity
int number_of_layers_;
// ! -----------------------------------------------------------------------------

// ! ---------------------- local error indicators --------------------------------
// ----- local error indicators (for each coarse grid element T) -------------
int max_loop_number = 10;
bool local_indicators_available_ = false;
// local coarse residual, i.e. H ||f||_{L^2(T)}
std::vector< std::vector< RangeType > > loc_coarse_residual_(max_loop_number);
// local coarse grid jumps (contribute to the total coarse residual)
std::vector< std::vector< RangeType > > loc_coarse_grid_jumps_(max_loop_number);
// local projection error (we project to get a globaly continous approximation)
std::vector< std::vector< RangeType > > loc_projection_error_(max_loop_number);
// local jump in the conservative flux
std::vector< std::vector< RangeType > > loc_conservative_flux_jumps_(max_loop_number);
// local approximation error
std::vector< std::vector< RangeType > > loc_approximation_error_(max_loop_number);
// local sum over the fine grid jumps (for a fixed subgrid that cooresponds with a coarse entity T)
std::vector< std::vector< RangeType > > loc_fine_grid_jumps_(max_loop_number);

std::vector< RangeType > total_coarse_residual_(max_loop_number);
std::vector< RangeType > total_projection_error_(max_loop_number);
std::vector< RangeType > total_coarse_grid_jumps_(max_loop_number);
std::vector< RangeType > total_conservative_flux_jumps_(max_loop_number);
std::vector< RangeType > total_approximation_error_(max_loop_number);
std::vector< RangeType > total_fine_grid_jumps_(max_loop_number);
std::vector< RangeType > total_estimated_H1_error_(max_loop_number);

bool repeat_algorithm_ = true;
int loop_number_ = 0;

#ifdef ADAPTIVE
double error_tolerance_;
#endif

#endif // MSFEM_GLOBALS_HH
