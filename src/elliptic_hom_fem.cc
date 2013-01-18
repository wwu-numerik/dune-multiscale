#include "common.hh"

// The following FEM code requires an access to the 'ModelProblemData' class,
// which provides us with information about f, A, \Omega, etc.


#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG


//! local (dune-multiscale) includes
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/multiscale/fem/algorithm.hh>

//! (very restrictive) homogenizer
#include <dune/multiscale/tools/homogenizer/elliptic_analytical_homogenizer.hh>
#include <dune/multiscale/tools/homogenizer/elliptic_homogenizer.hh>

//! the main FEM computation
template < class FEMTraits >
void algorithm_hom_fem(typename FEMTraits::GridPointerType& macro_grid_pointer,   // grid pointer that belongs to the macro grid
                       const std::string filename) {
  typedef FEMTraits FEM;
  using namespace Dune;

  const typename FEM::ModelProblemDataType problem_data;
  print_info(problem_data, DSC_LOG_INFO);
  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for the finite element problem
  typename FEM::GridPartType gridPart(*macro_grid_pointer);
  //! --------------------------------------------------------------------------------------

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  typename FEM::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  //! --------------------------------------------------------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const typename FEM::DiffusionType diffusion_op;

  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  Dune::RightHandSideAssembler< typename FEM::DiscreteFunctionType > rhsassembler;
  const typename FEM::FirstSourceType f;   // standard source f

  //! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
  const typename FEM::EllipticOperatorType discrete_elliptic_op(discreteFunctionSpace, diffusion_op);

  // unit cube grid for the computations of cell problems
  const std::string unit_cell_location = "../dune/multiscale/grids/cell_grids/unit_cube.dgf";
  // descretized homogenizer:

  typedef Dune::Homogenizer< typename FEM::GridType, typename FEM::DiffusionType > HomogenizerType;
  
  // to create an empty diffusion matrix that can be filled with constant values
  typedef Problem::ConstantDiffusionMatrix< typename FEM::FunctionSpaceType, typename HomogenizerType::HomTensorType >
     HomDiffusionType;
    
  const HomogenizerType disc_homogenizer(unit_cell_location);
  const typename HomogenizerType::HomTensorType A_hom = disc_homogenizer.getHomTensor(diffusion_op);
  const HomDiffusionType hom_diffusion_op(A_hom);
  
  //!TODO check: hatte nur 2 tmp parameter, Masse hinzugefUGT
  typedef DiscreteEllipticOperator< typename FEM::DiscreteFunctionType,
                                    HomDiffusionType, typename FEM::MassTermType > HomEllipticOperatorType;

  HomEllipticOperatorType hom_discrete_elliptic_op( discreteFunctionSpace, hom_diffusion_op); 
  
  typename FEM::FEMMatrix hom_stiff_matrix("homogenized stiffness matrix", discreteFunctionSpace, discreteFunctionSpace);
  
  typename FEM::DiscreteFunctionType hom_rhs("homogenized rhs", discreteFunctionSpace);
  hom_rhs.clear();
  
  //! solution vector
  // - By solution, we denote the (discrete) homogenized solution determined with FEM on the coarse scale and FEM for the cell problems
  typename FEM::DiscreteFunctionType homogenized_solution(filename + " Homogenized Solution", discreteFunctionSpace);
  homogenized_solution.clear();
  hom_discrete_elliptic_op.assemble_matrix(hom_stiff_matrix);
  
  constexpr int hmm_polorder = 2* FEM::DiscreteFunctionSpaceType::polynomialOrder + 2;
  rhsassembler.template assemble < hmm_polorder >(f, hom_rhs);

  // set Dirichlet Boundary to zero
  boundaryTreatment(hom_rhs);

  const typename FEM::InverseFEMMatrix hom_biCGStab(hom_stiff_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("global.cgsolver_verbose", false));
  hom_biCGStab(hom_rhs, homogenized_solution);
  
  // write FEM solution to a file and produce a VTK output
  // ---------------------------------------------------------------------------------
  
  // write the final (discrete) solution to a file
  std::string solution_file = (boost::format("/homogenized_solution_macro_refLevel_%d")
                                % DSC_CONFIG_GET("fem.grid_level", 4) ).str();
  DiscreteFunctionWriter(solution_file).append(homogenized_solution);

  // writing paraview data output
  // general output parameters
  Dune::myDataOutputParameters outputparam;

  // create and initialize output class
  typename FEM::IOTupleType hom_fem_solution_series(&homogenized_solution);
  outputparam.set_prefix((boost::format("/homogenized_solution")).str());
  typename FEM::DataOutputType homfemsol_dataoutput(homogenized_solution.space().gridPart().grid(),
                                                    hom_fem_solution_series, outputparam);
  homfemsol_dataoutput.writeData( 1.0 /*dummy*/, "homogenized-solution" );
  
  // ---------------------------------------------------------------------------------
  
}


int main(int argc, char** argv) {
  try {
    init(argc, argv);

    if ( !DSC_CONFIG_GET("problem.linear", 1 ) )
     DUNE_THROW(Dune::InvalidStateException, "Problem is declared to be nonlinear. Homogenizers are only available for linear homogenization problems.");
     
    namespace DSC = Dune::Stuff::Common;

    DSC_PROFILER.startTiming("total_cpu");

    const std::string path = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(path);

    // name of the error file in which the data will be saved
    std::string filename_;
    const Problem::ModelProblemData info;

    if ( !info.problemIsPeriodic() )
     DUNE_THROW(Dune::InvalidStateException, "Problem is declared to be non-periodic. Homogenizers are only available for periodic homogenization problems.");

    if ( DSC_CONFIG_GET("problem.stochastic_pertubation", false) )
     DUNE_THROW(Dune::InvalidStateException, "Homogenizers are only available for non-stochastically perturbed problems. Please, switch off the key 'problem.stochastic_pertubation'.");
    
    const std::string save_filename = std::string(path + "/logdata/ms.log.log");
    DSC_LOG_INFO << "LOG FILE " << std::endl << std::endl;
    std::cout << "Data will be saved under: " << save_filename << std::endl;
    DSC_LOG_INFO << "Solving the homogenized problem with a standard Finite Element method." << std::endl << std::endl;

    // refinement_level denotes the grid refinement level for the global problem, i.e. it describes 'H'
    const int refinement_level = DSC_CONFIG_GET("fem.grid_level", 4);

    // name of the grid file that describes the macro-grid:
    const std::string gridName = info.getMacroGridFile();
    DSC_LOG_INFO << "loading dgf: " << gridName << std::endl;

    // create a grid pointer for the DGF file belongig to the macro grid:
    FEMTraits::GridPointerType grid_pointer(gridName);

    // refine the grid 'starting_refinement_level' times:
    grid_pointer->globalRefine(refinement_level);

    algorithm_hom_fem<FEMTraits> (grid_pointer, filename_);

    const auto cpu_time = DSC_PROFILER.stopTiming("total_cpu") / 1000.f;
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "ms" << std::endl;
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
