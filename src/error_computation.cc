
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

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>


//! local (dune-multiscale) includes

#include <dune/multiscale/operators/disc_func_writer/discretefunctionwriter.hh>

#include <dune/multiscale/problems/elliptic_problems/model_problem_6/problem_specification.hh>
#include <dune/multiscale/operators/msfem_localproblems/localproblemsolver.hh>
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

std :: string location_fine_scale_reference_;
std :: string location_msfem_solution_;
std :: string location_homogenized_solution_;

int refinement_level_macrogrid_;
int refinement_level_finescale_referenceprob_;
int refinement_level_grid_T0_;
int refinement_level_homogenized_prob_;

//std :: ofstream data_file_; // file where we save the data

//! -----------------------------------------------------------------------------








void algorithm ( GridPointerType &macro_grid_pointer_msfem_sol, // grid pointer that belongs to the macro grid
                 GridPointerType &fine_macro_grid_pointer, // grid pointer that belongs to the fine macro grid (for reference computations)
                 GridPointerType &homog_macro_grid_pointer, // grid pointer that belongs to the macro grid for the homog. problem
                 GridPointerType &ref_simplex_grid_pointer, // grid pointer that belongs to the grid for the reference element T_0
                 std :: ofstream &data_file )
{


  //! ---- tools ----

  // model problem data
  Problem::ModelProblemData problem_info;
  std :: cout << "Loaded Problem Data." << std :: endl;

  L2Error< DiscreteFunctionType > l2error;

  // expensive hack to deal with discrete functions, defined on different grids
  ImprovedL2Error< DiscreteFunctionType > impL2error;

  //! ---------------------------- grid parts ----------------------------------------------


  // grid part for the global function space, required for MsFEM-macro-problem
  GridPartType gridPart( *macro_grid_pointer_msfem_sol);

  std :: cout << "Created Grid Partition for the MsFEM Macro Grid." << std :: endl;


  // grid part for the function space for T_0, required for the local HsFEM-problems
  GridPartType refSimplexGridPart ( *ref_simplex_grid_pointer );

  std :: cout << "Created Grid Partition for the Reference Simplex." << std :: endl;


  // grid part for the global function space, required for the detailed fine-scale computation (very high resolution)
  GridPartType gridPartFine( *fine_macro_grid_pointer );

  std :: cout << "Created Grid Partition for the Fine Macro Grid of the fine-scale reference Problem." << std :: endl;


  // grid part for the discrete function space that belongs to the available homogenized solution
  GridPartType gridPartHomog( *homog_macro_grid_pointer );

  std :: cout << "Created Grid Partition for the Grid of the homogenized Problem." << std :: endl;


  GridType &grid = gridPart.grid();
  GridType &gridFine = gridPartFine.grid();
  GridType &gridHomog = gridPartHomog.grid();

  std :: cout << "Created Grids." << std :: endl;

  //! --------------------------------------------------------------------------------------



  //! ------------------------- discrete function spaces -----------------------------------

  // the global-problem function space:
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  // the global-problem function space for the reference computation:
  DiscreteFunctionSpaceType finerDiscreteFunctionSpace( gridPartFine );

  // the global-problem function space for the solution of the homog. problem:
  DiscreteFunctionSpaceType homogDiscreteFunctionSpace( gridPartHomog );

  // the local-problem function space:
  DiscreteFunctionSpaceType refSimplexDiscreteFunctionSpace( refSimplexGridPart );

  std :: cout << "Created Discrete Function Spaces." << std :: endl;

  //! --------------------------------------------------------------------------------------



  //! solution vector
  // solution of the MsFEM
  DiscreteFunctionType msfem_solution( filename_ + " MsFEM Solution", discreteFunctionSpace );
  msfem_solution.clear();

  DiscreteFunctionReader discrete_function_reader_msfem( (location_msfem_solution_).c_str() );
  discrete_function_reader_msfem.open();

  discrete_function_reader_msfem.read( 0, msfem_solution );
  std :: cout << "MsFEM solution read." << std :: endl;




#ifdef HOMOGENIZEDSOL_AVAILABLE

  DiscreteFunctionType homogenized_solution( filename_ + " Homogenized Solution", homogDiscreteFunctionSpace );
  homogenized_solution.clear();

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_hom( (location_homogenized_solution_).c_str() );
  discrete_function_reader_hom.open();

  discrete_function_reader_hom.read( 0, homogenized_solution );
  std :: cout << "homogenized solution read." << std :: endl;

#endif


#ifdef FINE_SCALE_REFERENCE

  DiscreteFunctionType fine_scale_fem_solution( filename_ + " Reference FEM Solution", finerDiscreteFunctionSpace );
  fine_scale_fem_solution.clear();

  bool reader_is_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_ref( (location_fine_scale_reference_).c_str() );
  discrete_function_reader_ref.open();

  discrete_function_reader_ref.read( 0, fine_scale_fem_solution );
  std :: cout << "fine scale reference read." << std :: endl;

// end: FINE_SCALE_REFERENCE defined or not defined
#endif 


//noch per Hand die Daten eingetragen:
#ifdef HMM_REFERENCE
#if 0

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


  std :: string location_hmm_ref = "data/MsFEM/Model_Problem_1/Macro_10_Micro_8/msfem_solution_discFunc_refLevel_10";

  bool hmm_ref_reader_is_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_hmm_ref( (location_hmm_ref).c_str() );
  discrete_function_reader_hmm_ref.open();

  discrete_function_reader_hmm_ref.read( 0, hmm_reference_solution );
  std :: cout << "HMM reference read." << std :: endl;

#endif
#endif

  // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain local problem
  LocProbNumberingManagerType lp_num_manager( discreteFunctionSpace, refSimplexDiscreteFunctionSpace, filename_);


  std :: cout << std :: endl << "The L2 errors:" << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "The L2 errors:" << std :: endl << std :: endl; }

  //! ----------------- compute L2-errors -------------------


#ifdef FINE_SCALE_REFERENCE

  RangeType msfem_error;
  if ( refinement_level_finescale_referenceprob_ == refinement_level_homogenized_prob_ )
    { msfem_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(msfem_solution, fine_scale_fem_solution); }
  else
    { msfem_error = impL2error.norm_L2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(msfem_solution, fine_scale_fem_solution); }


  std :: cout << "|| u_msfem - u_fine_scale ||_L2 =  " << msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
         { data_file << "|| u_msfem - u_fine_scale ||_L2 =  " << msfem_error << std :: endl; }

#endif

//! STILL A TEST:
#if 0

  RangeType msfem_finescale_error = impL2error.error_L2_with_corrector
       < LocProbNumberingManagerType , 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(fine_scale_fem_solution, msfem_solution, lp_num_manager);

  std :: cout << "|| u_msfem_finescale - u_fine_scale ||_L2 = " << msfem_finescale_error << std :: endl << std :: endl;


#endif





#ifdef HOMOGENIZEDSOL_AVAILABLE

// not yet modified according to a generalized L2-error, here, homogenized_solution and fine_scale_fem_solution still need to be defined on the same grid!
  #ifdef FINE_SCALE_REFERENCE
  // RangeType hom_error = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( homogenized_solution, fine_scale_fem_solution );
  RangeType hom_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( homogenized_solution, fine_scale_fem_solution );

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




//! -------------------------------------------------------

}

int main(int argc, char** argv)
{

  Dune::MPIManager::initialize(argc, argv);

  // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
  refinement_level_macrogrid_ = 6;

  // grid refinement level for solving the cell problems, i.e. it describes 'h':
  refinement_level_grid_T0_ = 14;

  // grid refinement level of the solution of the fine-scale reference problem
  refinement_level_finescale_referenceprob_ = 14;
  //refinement_level_finescale_referenceprob_2_ = 4; // in case of two ref. problems.

  // grid refinement level of the solution of the homogenized problem
  refinement_level_homogenized_prob_ = 12;





  //! location of the fine scale reference solution (CHECKEN! - NOCH FEHLER!!!)

  char modeprob[50];
  //sprintf( modeprob, "Model_Problem_%d", problem_info.get_Number_of_Model_Problem() );
  std::string modeprob_s(modeprob);

  char reference_solution_directory[50];
  sprintf( reference_solution_directory, "/fem_fine_scale_%d", refinement_level_finescale_referenceprob_ );
  std::string reference_solution_directory_s(reference_solution_directory);

  char reference_solution_name[50];
  sprintf( reference_solution_name, "/finescale_solution_discFunc_refLevel_%d", refinement_level_finescale_referenceprob_ );
  std::string reference_solution_name_s(reference_solution_name);

  location_fine_scale_reference_ = "data/MsFEM/" + modeprob_s + reference_solution_directory_s + reference_solution_name_s;

  //! --------------------



  //! location of the homogenized solution

  //char modeprob[50];
  //sprintf( modeprob, "Model_Problem_%d", problem_info.get_Number_of_Model_Problem() );
  //std::string modeprob_s(modeprob);

  char homogenized_solution_directory[50];
  sprintf( homogenized_solution_directory, "/fem_homogenized_solution_%d", refinement_level_homogenized_prob_ );
  std::string homogenized_solution_directory_s(homogenized_solution_directory);

  char homogenized_solution_name[50];
  sprintf( homogenized_solution_name, "/homogenized_solution_discFunc_refLevel_%d", refinement_level_homogenized_prob_ );
  std::string homogenized_solution_name_s(homogenized_solution_name);

  location_homogenized_solution_ = "data/MsFEM/" + modeprob_s + homogenized_solution_directory_s + homogenized_solution_name_s;

  //! --------------------




  //! location of the msfem solution 

  //char modeprob[50];
  //sprintf( modeprob, "Model_Problem_%d", problem_info.get_Number_of_Model_Problem() );
  //std::string modeprob_s(modeprob);

  char msfem_solution_directory[50];
  sprintf( msfem_solution_directory, "macro_%d_micro_%d", refinement_level_macrogrid_, refinement_level_grid_T0_ );
  std::string msfem_solution_directory_s(msfem_solution_directory);
  filename_ = msfem_solution_directory_s;
  msfem_solution_directory_s = "/" + msfem_solution_directory_s;

  char msfem_solution_name[50];
  sprintf( msfem_solution_name, "/msfem_solution_discFunc_refLevel_%d", refinement_level_macrogrid_ );
  std::string msfem_solution_name_s( msfem_solution_name );

  location_msfem_solution_ = "data/MsFEM/" + modeprob_s + msfem_solution_directory_s + msfem_solution_name_s;

  //! --------------------





  // check directories for data output
  if (mkdir(("data/MsFEM/" + filename_).c_str() DIRMODUS) == -1)
   {
    std::cout << "Directory '" << "data/MsFEM/" << filename_ << "' exists. Loading data.";
   }
  else
   {
    std::cout << "Directory '" << "data/MsFEM/" << filename_ << "' does not exist!";
    abort();
   }




  std :: string save_filename = "data/MsFEM/" + filename_ + "/computed-errors.txt";
  std :: cout << "Data will be saved under: " << save_filename << std :: endl;

  // data for the model problem; the information manager
  // (see 'problem_specification.hh' for details)
  Problem::ModelProblemData info( filename_ );





  //name of the grid file that describes the macro-grid:
  std :: string macroGridName;
  info.getMacroGridFile( macroGridName );
  std :: cout << "loading dgf: " << macroGridName << std :: endl;

  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values for the parameters:

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer_msfem_sol( macroGridName );
  // refine the grid 'refinement_level_macrogrid_' times:
  macro_grid_pointer_msfem_sol->globalRefine( refinement_level_macrogrid_ );

  // create a finer GridPart for the fine-scale problem.
  // this shall be used as an approximation of the exact solution.
  GridPointerType fine_macro_grid_pointer( macroGridName );
  fine_macro_grid_pointer->globalRefine( refinement_level_finescale_referenceprob_ );

  // create a GridPart for the solution of the homogenized problem. 
  GridPointerType homog_macro_grid_pointer( macroGridName );
  homog_macro_grid_pointer->globalRefine( refinement_level_homogenized_prob_ );

    
  // after transformation, the cell problems are problems on the 0-centered unit cube [-½,½]²:
  std :: string RefElementName( "../dune/multiscale/grids/cell_grids/ref_simplex_2d.dgf" ); // --> the 0-centered unit cube, i.e. [-1/2,1/2]^2
 // std :: string RefElementName( "../dune/multiscale/grids/cell_grids/unit_cube_0_centered.dgf" );


  // to solve the local MsFEM problems, we always need a gridPart for T_0.
  // Here it is always the refernece simplex that needs to be used (after transformation, local problems are always formulated on such a grid )
  GridPtr< GridType > ref_simplex_grid_pointer( RefElementName );
  ref_simplex_grid_pointer->globalRefine( refinement_level_grid_T0_ );



  // to save all information in a file
  std :: ofstream data_file( (save_filename).c_str() );
  if (data_file.is_open())
            {
               data_file << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << "." << std :: endl << std :: endl;
               data_file << "For details on problem and method, see 'problem-info.txt'." << std :: endl << std :: endl;
               data_file << "Refinement level of finescale reference problem: " << refinement_level_finescale_referenceprob_ << std :: endl << std :: endl;
               data_file << "Refinement level of homogenized problem: " << refinement_level_homogenized_prob_ << std :: endl << std :: endl;
               data_file << "This file just contains the associated L2 and H1 errors." << std :: endl << std :: endl;
               data_file << std :: endl << std :: endl;
            }

  algorithm( macro_grid_pointer_msfem_sol,
             fine_macro_grid_pointer,
             homog_macro_grid_pointer,
             ref_simplex_grid_pointer,
             data_file );

  data_file.close();


  return 0;

}
