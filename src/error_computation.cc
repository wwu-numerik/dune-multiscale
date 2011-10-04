
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


// polynomial order of discrete space
#define POLORDER 1

// grid type
#define GRIDTYPE ALBERTAGRID
//#define GRIDTYPE UGGRID

// computational domain is a subset of \R^{GRIDDIM}
#define GRIDDIM 2
#define WORLDDIM GRIDDIM


#ifndef USE_GRAPE
#define USE_GRAPE HAVE_GRAPE
#endif

#define USE_TWISTFREE_MAPPER
#define VERBOSE false

#include <iostream>
#include <sstream>

// for creation of directories
#include <sys/types.h>
#include <sys/stat.h>
#define DIRMODUS ,0711

#include <stdio.h>
#include <stdlib.h>
//-----------------------------

//! is an exact solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro EXACTSOLUTION_AVAILABLE

//! is the homogenized solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro HOMOGENIZEDSOL_AVAILABLE
// (if HOMOGENIZEDSOL_AVAILABLE == true, it means that it is already computed.)
#define HOMOGENIZEDSOL_AVAILABLE


//! Do we have/want a fine-scale reference solution?
#define FINE_SCALE_REFERENCE


//! Do we have a HMM reference solution? (precomputed detailed HMM simulation)
// we might use a detailed HMM computation as a reference! (if it is available)
//#define HMM_REFERENCE


#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

// for yasp grid
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
// for ug grid
#include <dune/grid/io/file/dgfparser/dgfug.hh>
// for alu grid:
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
// for alberta grid:
#include <dune/grid/albertagrid/dgfparser.hh>


#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif


// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>



#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>




//! local (dune-multiscale) includes

#include <dune/multiscale/operators/disc_func_writer/discretefunctionwriter.hh>

#include <dune/multiscale/operators/meanvalue.hh>



using namespace Dune;

//! check for GridType:
#if GRIDTYPE==ALBERTAGRID
 typedef AlbertaGrid< GRIDDIM, WORLDDIM > GridType;
#elif GRIDTYPE==UGGRID
 // the Grid Type ( UG-Grid )
 typedef UGGrid< GRIDDIM > GridType;
#elif GRIDTYPE==YASPGRID
  // the Grid Type ( Yasp-Grid )
 typedef YaspGrid< GRIDDIM > GridType;
#elif GRIDTYPE==ALUGRID
    #ifdef CUBEGRID
        typedef ALUCubeGrid<GRIDDIM,GRIDDIM> GridType;
    #else
        typedef ALUSimplexGrid<GRIDDIM,GRIDDIM> GridType;
    #endif
#endif


//! --------- typedefs for the macro grid and the corresponding discrete space -------------

//Dune::InteriorBorder_Partition or Dune::All_Partition >?
//see: http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
typedef AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;

typedef GridPtr< GridType > GridPointerType;

typedef FunctionSpace < double , double , WORLDDIM , 1 > FunctionSpaceType;

//!-----------------------------------------------------------------------------------------




//! ---------  typedefs for the standard discrete function space (macroscopic) -------------

typedef FunctionSpaceType::DomainType DomainType; 

//! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
typedef FunctionSpaceType::RangeType RangeType;

//! defines the function space to which the numerical solution belongs to 
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace < FunctionSpaceType, GridPartType, 1 > //1=POLORDER
   DiscreteFunctionSpaceType;


typedef DiscreteFunctionSpaceType :: DomainFieldType TimeType;

typedef DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

typedef GridType :: Codim<0> :: Entity EntityType; 
typedef GridType :: Codim<0> :: EntityPointer EntityPointerType; 
typedef GridType :: Codim<0> :: Geometry EntityGeometryType; 
typedef GridType :: Codim<1> :: Geometry FaceGeometryType; 

typedef DiscreteFunctionSpaceType     :: BaseFunctionSetType      BaseFunctionSetType;

typedef CachingQuadrature < GridPartType , 0 > EntityQuadratureType;
typedef CachingQuadrature < GridPartType , 1 > FaceQuadratureType;

typedef DiscreteFunctionSpaceType :: DomainFieldType TimeType;
typedef DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;


typedef DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
typedef DiscreteFunctionType :: DofIteratorType DofIteratorType;

//!-----------------------------------------------------------------------------------------



//! ------------ cell problem solver and numbering manager -----------------------------------------

typedef LocalProblemNumberingManager< DiscreteFunctionSpaceType > LocProbNumberingManagerType;

//!-----------------------------------------------------------------------------------------








//! -------------------------------- important variables ---------------------------------

enum { dimension = GridType :: dimension};

//name of the error file in which the data will be saved
std :: string filename_;

int refinement_level_macrogrid_;
int refinement_level_referenceprob_;

//std :: ofstream data_file_; // file where we save the data

//! -----------------------------------------------------------------------------








void algorithm ( std :: string &RefElementName,
                 GridPointerType &macro_grid_pointer, // grid pointer that belongs to the macro grid
                 GridPointerType &fine_macro_grid_pointer, // grid pointer that belongs to the fine macro grid (for reference computations)
                 GridPointerType &ref_simplex_grid_pointer, // grid pointer that belongs to the grid for the refernce element
                 int refinement_difference, //refinement difference for the macro grid (problem-to-solve vs. reference problem)
                 std :: ofstream &data_file )
{


  //! ---- tools ----

  // model problem data
  Problem::ModelProblemData problem_info;

  L2Error< DiscreteFunctionType > l2error;

  // expensive hack to deal with discrete functions, defined on different grids
  ImprovedL2Error< DiscreteFunctionType > impL2error;

  //! ---------------------------- grid parts ----------------------------------------------

  // grid part for the global function space, required for MsFEM-macro-problem
  GridPartType gridPart( *macro_grid_pointer);

  // grid part for the function space for T_0, required for the local HsFEM-problems
  GridPartType refSimplexGridPart ( *ref_simplex_grid_pointer );

  // grid part for the global function space, required for the detailed fine-scale computation (very high resolution)
  GridPartType gridPartFine( *fine_macro_grid_pointer );

  GridType &grid = gridPart.grid();
  GridType &gridFine = gridPartFine.grid();

  //! --------------------------------------------------------------------------------------



  //! ------------------------- discrete function spaces -----------------------------------

  // the global-problem function space:
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  // the global-problem function space for the reference computation:
  DiscreteFunctionSpaceType finerDiscreteFunctionSpace( gridPartFine );

  // the local-problem function space:
  DiscreteFunctionSpaceType refSimplexDiscreteFunctionSpace( refSimplexGridPart );

  //! --------------------------------------------------------------------------------------





#ifdef HOMOGENIZEDSOL_AVAILABLE

  DiscreteFunctionType homogenized_solution( filename_ + " Homogenized Solution", finerDiscreteFunctionSpace );
  homogenized_solution.clear();

#endif


#ifdef FINE_SCALE_REFERENCE

  DiscreteFunctionType fem_solution( filename_ + " Reference FEM Solution", finerDiscreteFunctionSpace );
  fem_solution.clear();

  char modeprob[50];
  sprintf( modeprob, "/Model_Problem_%d", problem_info.get_Number_of_Model_Problem() );
  std::string modeprob_s(modeprob);

  char reference_solution_directory[50];
  sprintf( reference_solution_directory, "/reference_solution_ref_%d", refinement_level_referenceprob_ );
  std::string reference_solution_directory_s(reference_solution_directory);

  char reference_solution_name[50];
  sprintf( reference_solution_name, "/finescale_solution_discFunc_refLevel_%d", refinement_level_referenceprob_ );
  std::string reference_solution_name_s(reference_solution_name);

  std :: string location_fine_scale_ref = "data/MsFEM/" + modeprob_s + reference_solution_directory_s + reference_solution_name_s;

  bool reader_is_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_ref( (location_fine_scale_ref).c_str() );
  discrete_function_reader_ref.open();

  discrete_function_reader_ref.read( 0, fem_solution );
  std :: cout << "fine scale reference read." << std :: endl;

// end: FINE_SCALE_REFERENCE defined or not defined
#endif 


//noch per Hand die Daten eingetragen:
#ifdef HMM_REFERENCE

  int gridLevel_refHMM = 10; // Macro_'gridLevel'

  std :: string macroGridName_refHMM;
  problem_info.getMacroGridFile( macroGridName_refHMM );
  std :: cout << "loading dgf: " << macroGridName_refHMM << std :: endl;

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer_refHMM( macroGridName_refHMM );
  // refine the grid 'gridLevel_refHMM' times:
  macro_grid_pointer_refHMM->globalRefine( gridLevel_refHMM );

  GridPartType gridPart_refHMM( *macro_grid_pointer_refHMM);
  GridType &grid_refHMM = gridPart_refHMM.grid();

  DiscreteFunctionSpaceType discreteFunctionSpace_refHMM( gridPart_refHMM );

  DiscreteFunctionType hmm_reference_solution( filename_ + " Reference (HMM) Solution", discreteFunctionSpace_refHMM );
  hmm_reference_solution.clear();

#if 0
  char modeprob_name[50];
  sprintf( modeprob_name, "/Model_Problem_%d", problem_info.get_Number_of_Model_Problem() );
  std::string modeprob_name_s(modeprob_name);

  char reference_msfem_solution_directory[50];
  sprintf( reference_msfem_solution_directory, ".......", 10 );
  std::string reference_solution_directory_s(reference_solution_directory);

  char reference_solution_name[50];
  sprintf( reference_solution_name, "....", refinement_level_referenceprob_ );
  std::string reference_solution_name_s(reference_solution_name);

  std :: string location_hmm_ref = "data/MsFEM/" + modeprob_s + reference_solution_directory_s + reference_solution_name_s;
#endif

  std :: string location_hmm_ref = "data/MsFEM/Model_Problem_1/Macro_10_Micro_8/msfem_solution_discFunc_refLevel_10";

  bool hmm_ref_reader_is_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_hmm_ref( (location_hmm_ref).c_str() );
  discrete_function_reader_hmm_ref.open();

  discrete_function_reader_hmm_ref.read( 0, hmm_reference_solution );
  std :: cout << "HMM reference read." << std :: endl;

#endif

  // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain local problem
  LocProbNumberingManagerType lp_num_manager( discreteFunctionSpace, refSimplexDiscreteFunctionSpace, filename_);

  //! solution vector
  // solution of the heterogeneous multiscale finite element method, where we used the Newton method to solve the non-linear system of equations
  DiscreteFunctionType msfem_solution( filename_ + " MsFEM Solution", discreteFunctionSpace );
  msfem_solution.clear();



  std :: cout << std :: endl << "The L2 errors:" << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "The L2 errors:" << std :: endl << std :: endl; }

  //! ----------------- compute L2-errors -------------------


#ifdef ERROR_COMPUTATION

#ifdef FINE_SCALE_REFERENCE
  long double timeadapt = clock();

  RangeType msfem_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(msfem_solution, fem_solution);

  std :: cout << "|| u_msfem - u_fine_scale ||_L2 =  " << msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
         { data_file << "|| u_msfem - u_fine_scale ||_L2 =  " << msfem_error << std :: endl; }

  timeadapt = clock() - timeadapt;
  timeadapt = timeadapt / CLOCKS_PER_SEC;

  // if it took longer then 1 minute to compute the error:
  if ( timeadapt > 60 )
   {
     std :: cout << "WARNING! EXPENSIVE! Error assembled in " << timeadapt << "s." << std :: endl << std :: endl;

     if (data_file.is_open())
         { std :: cout << "WARNING! EXPENSIVE! Error assembled in " << timeadapt << "s." << std :: endl << std :: endl; }
   }
#endif

//! STILL A TEST:
#if 0

  RangeType msfem_finescale_error = impL2error.error_L2_with_corrector
       < LocProbNumberingManagerType , 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(fem_solution, msfem_solution, lp_num_manager);

  std :: cout << "|| u_msfem_finescale - u_fine_scale ||_L2 = " << msfem_finescale_error << std :: endl << std :: endl;


#endif





#ifdef HMM_REFERENCE
  long double timeadapthmmref = clock();

  RangeType msfem_vs_hmm_ref_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(msfem_solution,hmm_reference_solution);

  std :: cout << "|| u_msfem - u_hmm_ref ||_L2 =  " << msfem_vs_hmm_ref_error << std :: endl << std :: endl;
  if (data_file.is_open())
         { data_file << "|| u_msfem - u_hmm_ref ||_L2 =  " << msfem_vs_hmm_ref_error << std :: endl; }

  timeadapthmmref = clock() - timeadapthmmref;
  timeadapthmmref = timeadapthmmref / CLOCKS_PER_SEC;

  // if it took longer then 1 minute to compute the error:
  if ( timeadapthmmref > 60 )
   {
     std :: cout << "WARNING! EXPENSIVE! Error assembled in " << timeadapthmmref << "s." << std :: endl << std :: endl;

     if (data_file.is_open())
         { std :: cout << "WARNING! EXPENSIVE! Error assembled in " << timeadapthmmref << "s." << std :: endl << std :: endl; }
   }
#endif


#ifdef HOMOGENIZEDSOL_AVAILABLE
// not yet modified according to a generalized L2-error, here, homogenized_solution and fem_solution still need to be defined on the same grid!
  #ifdef FINE_SCALE_REFERENCE
  // RangeType hom_error = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( homogenized_solution, fem_solution );
  RangeType hom_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( homogenized_solution, fem_solution );

  std :: cout << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std :: endl; }
  #endif

  //RangeType hom_msfem_error = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( msfem_solution, homogenized_solution );
  RangeType hom_msfem_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( msfem_solution, homogenized_solution );

  std :: cout << "|| u_hom - u_msfem ||_L2 =  " << hom_msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_hom - u_msfem ||_L2 =  " << hom_msfem_error << std :: endl; }
#endif


#ifdef EXACTSOLUTION_AVAILABLE
  RangeType exact_msfem_error = l2error.norm< ExactSolutionType >( u, msfem_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_msfem - u_exact ||_L2 =  " << exact_msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_msfem - u_exact ||_L2 =  " << exact_msfem_error << std :: endl; }

 #ifdef FINE_SCALE_REFERENCE
  RangeType fem_error = l2error.norm< ExactSolutionType >( u, fem_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl; }
 #endif

#endif

// endif for macro ERROR_COMPUTATION
#endif


//! -------------------------------------------------------

}

int main(int argc, char** argv)
{

  if(argc != 2)
  {
    fprintf(stderr,"usage: %s <starting_level_for_grid_refinement> \n",argv[0]);
    exit(1);
  }

 #ifndef LINEAR_PROBLEM
  std :: cout << "Nonlinear case not implemented, please define LINEAR_PROBLEM." <<std :: endl;
 #endif

  Dune::MPIManager::initialize(argc, argv);

  // name of the file in which you want to save the data:
  std :: cout << "Enter name for data directory: ";
  std :: cin >> filename_;

  // generate directories for data output
  if (mkdir(("data/MsFEM/" + filename_).c_str() DIRMODUS) == -1)
   {
    std::cout << "Directory already exists! Overwrite? y/n: ";
    char answer;
    std :: cin >> answer;
    if (!(answer=='y'))
     {std :: abort();}
   }
  else
   {
     mkdir(("data/MsFEM/" + filename_).c_str() DIRMODUS);
   }


  std :: string save_filename = "data/MsFEM/" + filename_ + "/computed-errors.txt";
  std :: cout << "Data will be saved under: " << save_filename << std :: endl;

  // data for the model problem; the information manager
  // (see 'problem_specification.hh' for details)
  Problem::ModelProblemData info( filename_ );

  //epsilon is specified in ModelProblemData, which is specified in problem_specification.hh
  epsilon_ = info.getEpsilon();

  //estimated epsilon (specified in ModelProblemData)
  epsilon_est_ = info.getEpsilonEstimated();

  //edge length of the cells in the cells, belonging to the cell problems
  //note that (delta/epsilon_est) needs to be a positive integer!
  delta_ = info.getDelta();

  // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
  refinement_level_macrogrid_ = atoi( argv[ 1 ] );

  // grid refinement level for solving the cell problems, i.e. it describes 'h':
  int refinement_level_grid_T0;
  std :: cout << "Enter refinement level for the local problems: ";
  std :: cin >> refinement_level_grid_T0;
  std :: cout << std :: endl;


  // (starting) grid refinement level for solving the reference problem
  refinement_level_referenceprob_ = info.getRefinementLevelReferenceProblem();
  // in general: for the homogenized case = 11 and for the high resolution case = 14
  // Note that this depends on the model problem!
#ifndef FINE_SCALE_REFERENCE
  refinement_level_referenceprob_ = 0;
#endif


  // how many times finer do we solve the reference problem (it is either a homogenized problem or the exact problem with a fine-scale resolution)
  int refinement_difference_for_referenceproblem = refinement_level_referenceprob_ - refinement_level_macrogrid_;


  //name of the grid file that describes the macro-grid:
  std :: string macroGridName;
  info.getMacroGridFile( macroGridName );
  std :: cout << "loading dgf: " << macroGridName << std :: endl;

  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values for the parameters:

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer( macroGridName );
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine( refinement_level_macrogrid_ );

  // create a finer GridPart for either the homogenized or the fine-scale problem.
  // this shall be used to compute an approximation of the exact solution.
  GridPointerType fine_macro_grid_pointer( macroGridName );
  // refine the grid 'starting_refinement_level_reference' times:
  fine_macro_grid_pointer->globalRefine( refinement_level_referenceprob_ );


  // after transformation, the cell problems are problems on the 0-centered unit cube [-½,½]²:
  std :: string RefElementName( "../dune/multiscale/grids/cell_grids/ref_simplex_2d.dgf" ); // --> the 0-centered unit cube, i.e. [-1/2,1/2]^2
 // std :: string RefElementName( "../dune/multiscale/grids/cell_grids/unit_cube_0_centered.dgf" );


  // to solve the local MsFEM problems, we always need a gridPart for T_0.
  // Here it is always the refernece simplex that needs to be used (after transformation, local problems are always formulated on such a grid )
  GridPtr< GridType > ref_simplex_grid_pointer( RefElementName );
  ref_simplex_grid_pointer->globalRefine( refinement_level_grid_T0 );


  // to save all information in a file
  std :: ofstream data_file( (save_filename).c_str() );
  if (data_file.is_open())
            {
               data_file << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << "." << std :: endl << std :: endl;
#ifdef LINEAR_PROBLEM
               data_file << "Problem is declared as being LINEAR." << std :: endl;
#else
               data_file << "Problem is declared as being NONLINEAR." << std :: endl;
#endif
#ifdef EXACTSOLUTION_AVAILABLE
               data_file << "Exact solution is available." << std :: endl << std :: endl;
#else
               data_file << "Exact solution is not available." << std :: endl << std :: endl;
#endif
               data_file << "Computations were made for:" << std :: endl << std :: endl;
               data_file << "Refinement Level for (uniform) Macro Grid = " << refinement_level_macrogrid_ << std :: endl;
               data_file << "Refinement Level for Micro Grid (grid on macro entity) = " << refinement_level_grid_T0 << std :: endl << std :: endl;
#ifdef PGF
               data_file << "We use MsFEM in Petrov-Galerkin formulation." << std :: endl;
#else
               data_file << "We use MsFEM in its standard formulation." << std :: endl;
#endif
               data_file << "Cell problems are solved and saved (in a pre-process)." << std :: endl << std :: endl;
               data_file << "Epsilon = " << epsilon_ << std :: endl;
               data_file << "Estimated Epsilon = " << epsilon_est_ << std :: endl;
               data_file << "Delta (edge length of cell-cube) = " << delta_ << std :: endl;
#ifdef STOCHASTIC_PERTURBATION
               data_file << std :: endl << "Stochastic perturbation added. Variance = " << VARIANCE << std :: endl;
#endif
               data_file << std :: endl << std :: endl;
            }

     algorithm( RefElementName,
                macro_grid_pointer,
                fine_macro_grid_pointer,
                ref_simplex_grid_pointer,
                refinement_difference_for_referenceproblem,
                data_file );
    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the normal macro problem


   long double cpu_time = clock();
   cpu_time = cpu_time / CLOCKS_PER_SEC;
   std :: cout << "Total runtime of the program: " << cpu_time << "s" << std :: endl;

   if (data_file.is_open())
    {
      data_file << "Total runtime of the program: " << cpu_time << "s" << std :: endl;
    }

  data_file.close();

  return 0;

}
