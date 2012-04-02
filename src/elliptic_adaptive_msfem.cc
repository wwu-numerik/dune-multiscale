
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
#include <vector>

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
#define EXACTSOLUTION_AVAILABLE


#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

// for yasp grid
// #include <dune/grid/io/file/dgfparser/dgfyasp.hh>
// for ug grid
// #include <dune/grid/io/file/dgfparser/dgfug.hh>
// for alu grid:
// #include <dune/grid/io/file/dgfparser/dgfalu.hh>
// for alberta grid:
#include <dune/grid/albertagrid/dgfparser.hh>


#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#define PGF

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>


// dune-fem includes:

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>


//! local (dune-multiscale) includes
//#include <dune/multiscale/problems/elliptic_problems/model_problem_toy/problem_specification.hh>
//#include <dune/multiscale/problems/elliptic_problems/model_problem_easy/problem_specification.hh>
#include <dune/multiscale/problems/elliptic_problems/model_problem_9/problem_specification.hh>

#include <dune/multiscale/tools/solver/FEM/fem_solver.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_solver.hh>
#include <dune/multiscale/tools/misc/h1error.hh>

#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>

//! Bei Gelegenheit LOESCHEN:
#include <dune/multiscale/tools/meanvalue.hh>


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

//! ----- typedefs for the macro grid and the corresponding discrete space -----

//Dune::InteriorBorder_Partition or Dune::All_Partition >?
//see: http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
typedef AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;

typedef GridPtr< GridType > GridPointerType;

typedef FunctionSpace < double , double , WORLDDIM , 1 > FunctionSpaceType;

//!-----------------------------------------------------------------------------



//! --------- typedefs for the coefficient and data functions ------------------

// type of first source term (right hand side of differential equation or type of 'f')
typedef Problem::FirstSource< FunctionSpaceType > FirstSourceType;

// type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
typedef Problem::Diffusion< FunctionSpaceType > DiffusionType;

// default type for any missing coefficient function (e.g. advection,...)
typedef Problem::DefaultDummyFunction< FunctionSpaceType > DefaultDummyFunctionType;

#ifdef EXACTSOLUTION_AVAILABLE
// type of exact solution (in general unknown)
typedef Problem::ExactSolution< FunctionSpaceType > ExactSolutionType;
typedef DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
  DiscreteExactSolutionType; //for data output with paraview or grape
#endif

//!-----------------------------------------------------------------------------




//! ----  typedefs for the standard discrete function space (macroscopic) -----

typedef FunctionSpaceType::DomainType DomainType; 

//! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
typedef FunctionSpaceType::RangeType RangeType;

//! defines the function space to which the numerical solution belongs to 
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace < FunctionSpaceType, GridPartType, 1 > //1=POLORDER
   DiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

//!-----------------------------------------------------------------------------









//! ---------------------- important variables ---------------------------------

enum { dimension = GridType :: dimension};

//name of the error file in which the data will be saved
std :: string filename_;
std :: string path_;

int total_refinement_level_;
int coarse_grid_level_;

//! -----------------------------------------------------------------------------



//! ------------------ typedefs and classes for data output ---------------------

typedef Tuple<DiscreteFunctionType*> IOTupleType;
typedef DataOutput< GridType, IOTupleType> DataOutputType;



#ifdef EXACTSOLUTION_AVAILABLE
// just for the discretized exact solution (in case it is available)
typedef Tuple<DiscreteExactSolutionType*> ExSolIOTupleType;
// just for the discretized exact solution (in case it is available)
typedef DataOutput<GridType, ExSolIOTupleType> ExSolDataOutputType;
#endif


// define output traits
struct myDataOutputParameters : public DataOutputParameters {

public:

  std::string my_prefix_;
  std::string my_path_;

  void set_prefix( std::string my_prefix )
    {
      my_prefix_ = my_prefix;
      // std :: cout << "Set prefix. my_prefix_ = " << my_prefix_ << std :: endl;
    }

  void set_path( std::string my_path )
    {
      my_path_ = my_path;
    }

  // base of file name for data file
  std::string prefix() const
    {
      if (my_prefix_ == "")
        return "solutions";
      else
        return my_prefix_;
    }

  // path where the data is stored
  std::string path() const 
    {
      if (my_path_ == "")
        return "data_output_msfem";
      else
        return my_path_;

    }


  // format of output:
  int outputformat() const
    {
      //return 0; // GRAPE (lossless format)
      return 1; // VTK
      //return 2; // VTK vertex data
      //return 3; // gnuplot
    }


};

//!---------------------------------------------------------------------------------------



void algorithm ( GridPointerType &macro_grid_pointer, // grid pointer that belongs to the macro grid
                 std :: ofstream &data_file )
{

  //! ---- tools ----

  // model problem data
  Problem::ModelProblemData problem_info;

  L2Error< DiscreteFunctionType > l2error;
  H1Error< DiscreteFunctionType > h1error;

  // expensive hack to deal with discrete functions, defined on different grids
  ImprovedL2Error< DiscreteFunctionType > impL2error;

  //! ---------------------------- grid parts ----------------------------------------------

  // grid part for the global function space, required for MsFEM-macro-problem
  GridPartType gridPart( *macro_grid_pointer);

  GridType &grid = gridPart.grid();

  //! --------------------------------------------------------------------------------------


  //! ------------------------- discrete function spaces -----------------------------------

  // the global-problem function space:
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  //! --------------------------------------------------------------------------------------



  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  DiffusionType diffusion_op;

  // define (first) source term:
  FirstSourceType f; // standard source f

  // exact solution unknown?
#ifdef EXACTSOLUTION_AVAILABLE
  ExactSolutionType u;
  DiscreteExactSolutionType discrete_exact_solution( "discrete exact solution ", u, gridPart );
#endif

  
  //! ---------------------------- general output parameters ------------------------------
  
  // general output parameters
  myDataOutputParameters outputparam;
  outputparam.set_path( path_ );

  // sequence stamp
  std::stringstream outstring;
  
  //! --------------------------------------------------------------------------------------


  
  //! -------------------------- writing data output Exact Solution ------------------------

#ifdef EXACTSOLUTION_AVAILABLE
  // --------- data output discrete exact solution --------------

  // create and initialize output class
  ExSolIOTupleType exact_solution_series( &discrete_exact_solution );
  outputparam.set_prefix("exact_solution");
  ExSolDataOutputType exactsol_dataoutput( gridPart.grid(), exact_solution_series, outputparam );

  // write data
  outstring << "exact-solution";
  exactsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());
  // -------------------------------------------------------
#endif
  
  //! --------------------------------------------------------------------------------------
  
  
  
  
  //! ---------------------- solve MsFEM problem ---------------------------

  // coarse space
#if 1
  Problem::ModelProblemData info( filename_ );
  std :: string macroGridName;
  info.getMacroGridFile( macroGridName );

  GridPointerType macro_grid_pointer_coarse( macroGridName );
  macro_grid_pointer_coarse->globalRefine( coarse_grid_level_ );
  GridPartType gridPart_coarse( *macro_grid_pointer_coarse);

  GridType &grid_coarse = gridPart_coarse.grid();
  DiscreteFunctionSpaceType discreteFunctionSpace_coarse( gridPart_coarse );
#endif

  //! solution vector
  // solution of the standard finite element method
  DiscreteFunctionType msfem_solution( filename_ + " MsFEM Solution", discreteFunctionSpace );
  msfem_solution.clear();

  DiscreteFunctionType coarse_part_msfem_solution( filename_ + " Coarse Part MsFEM Solution", discreteFunctionSpace );
  coarse_part_msfem_solution.clear();

  DiscreteFunctionType fine_part_msfem_solution( filename_ + " Fine Part MsFEM Solution", discreteFunctionSpace );
  fine_part_msfem_solution.clear();

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  int number_of_level_host_entities = grid.size( coarse_grid_level_, 0 /*codim*/ );  
  int coarse_level_fine_level_difference = grid.maxLevel() - grid_coarse.maxLevel();
  
  MacroMicroGridSpecifier specifier( number_of_level_host_entities , coarse_level_fine_level_difference );
  for ( int i = 0; i < number_of_level_host_entities; i+=1 )
    { specifier.setLayer( i , 30 ); }

#if 1
  // just for Dirichlet zero-boundary condition
  Elliptic_MsFEM_Solver< DiscreteFunctionType > msfem_solver( discreteFunctionSpace, data_file, path_ );
  msfem_solver.solve_dirichlet_zero( diffusion_op, f, discreteFunctionSpace_coarse, specifier,
                                     coarse_part_msfem_solution, fine_part_msfem_solution, msfem_solution );
#endif

  //! ----------------------------------------------------------------------

  std :: cout << "Data output for MsFEM Solution." << std :: endl;
  
  //! ----------------- writing data output MsFEM Solution -----------------
  
  // --------- VTK data output for MsFEM solution --------------------------

  // create and initialize output class
  IOTupleType msfem_solution_series( &msfem_solution );
  outputparam.set_prefix("msfem_solution");
  DataOutputType msfem_dataoutput( gridPart.grid(), msfem_solution_series, outputparam );

  // write data
  outstring << "msfem_solution";
  msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());


  // create and initialize output class
  IOTupleType coarse_msfem_solution_series( &coarse_part_msfem_solution );
  outputparam.set_prefix("coarse_part_msfem_solution");
  DataOutputType coarse_msfem_dataoutput( gridPart.grid(), coarse_msfem_solution_series, outputparam );

  // write data
  outstring << "coarse_msfem_solution";
  coarse_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());



  // create and initialize output class
  IOTupleType fine_msfem_solution_series( &fine_part_msfem_solution );
  outputparam.set_prefix("fine_part_msfem_solution");
  DataOutputType fine_msfem_dataoutput( gridPart.grid(), fine_msfem_solution_series, outputparam );

  // write data
  outstring << "fine_msfem_solution";
  fine_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());


  // ---------------------------------------------------------------------- 
  
  // ---------------------- write discrete msfem solution to file ---------
  bool writer_is_open = false;

  char fname[50];
  sprintf( fname, "/msfem_solution_discFunc_refLevel_%d_%d", total_refinement_level_, coarse_grid_level_);
  std :: string fname_s( fname );

  std :: string location = path_ + fname_s;
  DiscreteFunctionWriter dfw( (location).c_str() );
  writer_is_open = dfw.open();
  if ( writer_is_open )
    dfw.append( msfem_solution );

  //! --------------------------------------------------------------------
    


  //! ---------------------- solve FEM problem ---------------------------

  //! solution vector
  // solution of the standard finite element method
  DiscreteFunctionType fem_solution( filename_ + " FEM Solution", discreteFunctionSpace );
  fem_solution.clear();

  // just for Dirichlet zero-boundary condition
  Elliptic_FEM_Solver< DiscreteFunctionType > fem_solver( discreteFunctionSpace,  data_file );
  fem_solver.solve_dirichlet_zero( diffusion_op, f, fem_solution );

  //! ----------------------------------------------------------------------

   std :: cout << "Data output for FEM Solution." << std :: endl;
 
  //! -------------------------- writing data output FEM Solution ----------

  // ------------- VTK data output for FEM solution --------------

  // create and initialize output class
  IOTupleType fem_solution_series( &fem_solution );
  outputparam.set_prefix("fem_solution");
  DataOutputType fem_dataoutput( gridPart.grid(), fem_solution_series, outputparam );

  // write data
  outstring << "fem_solution";
  fem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------------
 
   
   
  // ------------- write discrete fem solution to file -----------

  writer_is_open = false;

  char fem_fname[50];
  sprintf( fem_fname, "/fem_solution_discFunc_refLevel_%d", total_refinement_level_ );
  std :: string fem_fname_s( fem_fname );

  std :: string fine_location = path_ + fem_fname_s;
  DiscreteFunctionWriter fem_dfw( (fine_location).c_str() );
  writer_is_open = fem_dfw.open();
  if ( writer_is_open )
    fem_dfw.append( fem_solution );
  
  // -------------------------------------------------------------

  //! ------------------------------------------------------------
    
  
    
  std :: cout << std :: endl << "The L2 errors:" << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "The L2 errors:" << std :: endl << std :: endl; }

  //! ----------------- compute L2- and H1- errors -------------------

#ifdef EXACTSOLUTION_AVAILABLE

  RangeType fem_error = l2error.norm< ExactSolutionType >( u, fem_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl; }

  RangeType h1_fem_error(0.0);
  h1_fem_error = h1error.semi_norm< ExactSolutionType >( u, fem_solution );
  h1_fem_error += fem_error;

  std :: cout << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std :: endl << std :: endl; }

#endif

#ifdef EXACTSOLUTION_AVAILABLE

  RangeType msfem_error = l2error.norm< ExactSolutionType >( u, msfem_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std :: endl; }

  RangeType h1_msfem_error(0.0);
  h1_msfem_error = h1error.semi_norm< ExactSolutionType >( u, msfem_solution );
  h1_msfem_error += msfem_error;

  std :: cout << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std :: endl; }

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

  Dune::MPIManager::initialize(argc, argv);

  // name of the file in which you want to save the data:
  std :: cout << "Enter name for data directory: ";
  std :: cin >> filename_;
  
  // generate directories for data output
  if (mkdir((path_).c_str() DIRMODUS) == -1)
   {
    std::cout << "Directory already exists! Overwrite? y/n: ";
    char answer;
    std :: cin >> answer;
    if (!(answer=='y'))
     {std :: abort();}
   }
  else
   {
     mkdir((path_).c_str() DIRMODUS);
   }

  path_ = "data/MsFEM/" + filename_;
  
  std :: cout << "Enter coarse grid level: ";
  std :: cin >> coarse_grid_level_;
  if ( coarse_grid_level_ > atoi( argv[ 1 ] ))
   {
     std :: cout << "Coarse grid level must be smaller than the starting refinment level." << std :: endl;
     abort();
   }
   
  std :: string save_filename = path_ + "/problem-info.txt";
  std :: cout << "Data will be saved under: " << save_filename << std :: endl;

  // data for the model problem; the information manager
  // (see 'problem_specification.hh' for details)
  Problem::ModelProblemData info( filename_ );

  // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
  total_refinement_level_ = atoi( argv[ 1 ] );

  //name of the grid file that describes the macro-grid:
  std :: string macroGridName;
  info.getMacroGridFile( macroGridName );
  std :: cout << "loading dgf: " << macroGridName << std :: endl;

  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values for the parameters:

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer( macroGridName );
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine( total_refinement_level_ );


  // to save all information in a file
  std :: ofstream data_file( (save_filename).c_str() );
  if (data_file.is_open())
            {
               data_file << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << "." << std :: endl << std :: endl;
               data_file << "Computations were made for:" << std :: endl << std :: endl;
               data_file << "Refinement Level for (uniform) Macro Grid = " << total_refinement_level_ << std :: endl;
               data_file << std :: endl << std :: endl;
            }

    algorithm( macro_grid_pointer, data_file );
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
