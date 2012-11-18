#include "common.hh"

// #define ADAPTIVE

#ifndef ADAPTIVE
 #define UNIFORM
#endif

#if HAVE_GRAPE
 #include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#define PGF

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>


// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/stuff/common/filesystem.hh>
//! local (dune-multiscale) includes
#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/multiscale/tools/solver/FEM/fem_solver.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_solver.hh>
#include <dune/multiscale/tools/misc/h1error.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/tools/improved_l2error.hh>
#include <dune/multiscale/tools/errorestimation/MsFEM/msfem_elliptic_error_estimator.hh>

//!-----------------------------------------------------------------------------

#include <dune/multiscale/msfem/msfem_traits.hh>

#include "msfem_globals.hh"
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/walk.hh>
#include <dune/stuff/grid/walk_functors.hh>
#include <dune/stuff/aliases.hh>

template<class DFSpace >
struct LGFunctor
{
  LGFunctor(const DFSpace& df_space)
    : df_space_(df_space)
  {}

  template< class Entity >
  void operator()(const Entity& ent, const int /*ent_idx*/) {
    auto lg_points = df_space_.lagrangePointSet(ent);
    for(const auto i : DSC::valueRange(lg_points.size())) {
      auto key = lg_points.localKey(i);
      if (key.codim() == Entity::dimension) {
        auto geo = ent.geometry();
        auto kk = lg_points.point(i);
        auto ll = geo.global(kk);
        DSC_LOG_DEBUG << kk << " (" << ll << ") ";
      }
    }
    DSC_LOG_DEBUG << "\n" << std::endl;
  }
  const DFSpace& df_space_;
};

//!---------------------------------------------------------------------------------------
//! algorithm
void algorithm(const std::string& macroGridName) {
  using namespace Dune;
  DSC_LOG_INFO << "loading dgf: " << macroGridName << std::endl;
  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
  // for the parameters:
  // create a grid pointer for the DGF file belongig to the macro grid:
  MsfemTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine(coarse_grid_level_);
  //! ---- tools ----
  // model problem data
  Problem::ModelProblemData problem_info;
  L2Error< MsfemTraits::DiscreteFunctionType > l2error;
  // expensive hack to deal with discrete functions, defined on different grids
  Dune::ImprovedL2Error< MsfemTraits::DiscreteFunctionType > DUNE_UNUSED(impL2error);

  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for MsFEM-macro-problem
  MsfemTraits::GridPartType gridPart(*macro_grid_pointer);
  MsfemTraits::GridType& grid = gridPart.grid();
  //! --------------------------------------------------------------------------------------

  // coarse grid
  MsfemTraits::GridPointerType macro_grid_pointer_coarse(macroGridName);
  macro_grid_pointer_coarse->globalRefine(coarse_grid_level_);
  MsfemTraits::GridPartType gridPart_coarse(*macro_grid_pointer_coarse);
  MsfemTraits::GridType& grid_coarse = gridPart_coarse.grid();

  grid.globalRefine(total_refinement_level_ - coarse_grid_level_);

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  MsfemTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  MsfemTraits::DiscreteFunctionSpaceType discreteFunctionSpace_coarse(gridPart_coarse);

  //! ---------------------- solve MsFEM problem ---------------------------
  //! solution vector
  // solution of the standard finite element method
  MsfemTraits::DiscreteFunctionType msfem_solution(filename_ + " MsFEM Solution", discreteFunctionSpace);
  msfem_solution.clear();

  MsfemTraits::DiscreteFunctionType coarse_part_msfem_solution(filename_ + " Coarse Part MsFEM Solution", discreteFunctionSpace);
  coarse_part_msfem_solution.clear();

  MsfemTraits::DiscreteFunctionType fine_part_msfem_solution(filename_ + " Fine Part MsFEM Solution", discreteFunctionSpace);
  fine_part_msfem_solution.clear();

  const int number_of_level_host_entities = grid_coarse.size(0 /*codim*/);
  const int coarse_level_fine_level_difference = grid.maxLevel() - grid_coarse.maxLevel();

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MsfemTraits::MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
  for (int i = 0; i < number_of_level_host_entities; i += 1)
  {
    specifier.setLayer(i, number_of_layers_);
  }

  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const MsfemTraits::DiffusionType diffusion_op;
  // define (first) source term:
  const MsfemTraits::FirstSourceType f; // standard source f

  //! create subgrids:
  const bool silence = false;
  {
    MsfemTraits::SubGridListType subgrid_list(specifier, silence);

    auto walk = DSG::make_gridwalk(grid.leafView());
    //! GridWalk functor that refines all entitites above given volume
    LGFunctor<MsfemTraits::DiscreteFunctionSpaceType> fun(discreteFunctionSpace);
    walk(fun);
  }
  //! -------------------------- writing data output FEM Solution ----------


} // function algorithm

int main(int argc, char** argv) {
  try {
    init(argc, argv);
    namespace DSC = Dune::Stuff::Common;
    //!TODO include base in config
    DSC_PROFILER.startTiming("msfem_all");

    const std::string datadir = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(datadir);

    DSC_LOG_INFO << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log")
                            % datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");

    // syntax: info_from_par_file / default  / validation of the value

    // coarse_grid_level denotes the (starting) grid refinement level for the global coarse scale problem, i.e. it describes 'H'
    coarse_grid_level_ = DSC_CONFIG_GETV( "msfem.coarse_grid_level", 4, DSC::ValidateLess< int >( -1 ) );

    // syntax: info_from_par_file / default
    number_of_layers_ = DSC_CONFIG_GET("msfem.oversampling_layers", 4);

    #ifdef ADAPTIVE
    error_tolerance_ = DSC_CONFIG_GET("msfem.error_tolerance", 1e-6);
    #endif   // ifdef ADAPTIVE

    // data for the model problem; the information manager
    // (see 'problem_specification.hh' for details)
    const Problem::ModelProblemData info;

    // total_refinement_level denotes the (starting) grid refinement level for the global fine scale problem, i.e. it describes 'h'
    total_refinement_level_
      = DSC_CONFIG_GETV( "msfem.fine_grid_level", 4, DSC::ValidateLess< int >(coarse_grid_level_-1) );

    // name of the grid file that describes the macro-grid:
    const std::string macroGridName = info.getMacroGridFile();

    algorithm(macroGridName);

    const auto cpu_time = DSC_PROFILER.stopTiming("msfem_all");
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "ms" << std::endl;
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
