#ifdef HAVE_CMAKE_CONFIG
#  include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#  include <config.h>
#endif

// polynomial order of discrete space
#define POLORDER 1

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


#define ADAPTIVE

#ifndef ADAPTIVE
#define UNIFORM
#endif


//! is an exact solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro EXACTSOLUTION_AVAILABLE
// #define EXACTSOLUTION_AVAILABLE


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
#include <dune/fem/misc/h1norm.hh>


//! local (dune-multiscale) includes
//#include <dune/multiscale/problems/elliptic_problems/model_problem_toy/problem_specification.hh>
//#include <dune/multiscale/problems/elliptic_problems/model_problem_easy/problem_specification.hh>
#include <dune/multiscale/problems/elliptic_problems/model_problem_10/problem_specification.hh>

#include <dune/multiscale/tools/solver/FEM/fem_solver.hh>

#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_solver.hh>
#include <dune/multiscale/tools/misc/h1error.hh>

#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>

//! Bei Gelegenheit LOESCHEN:
#include <dune/multiscale/tools/meanvalue.hh>

#include <dune/multiscale/tools/errorestimation/MsFEM/elliptic_error_estimator.hh>

using namespace Dune;

//! ----- typedefs for the macro grid and the corresponding discrete space -----
typedef GridSelector::GridType
    GridType;
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
typedef LagrangeDiscreteFunctionSpace < FunctionSpaceType, GridPartType, 1 > // 1 = POLORDER
   DiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

//!-----------------------------------------------------------------------------




//!------------------------- for adaptive grid refinement ---------------------------------

//! For adaption:

//! type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunctionType >
  RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< GridType, RestrictProlongOperatorType >
  AdaptationManagerType;

//!---------------------------------------------------------------------------------------


typedef MacroMicroGridSpecifier< DiscreteFunctionSpaceType > MacroMicroGridSpecifierType;
typedef SubGrid< GridType :: dimension , GridType > SubGridType; 
typedef SubGridList< DiscreteFunctionType, SubGridType, MacroMicroGridSpecifierType > SubGridListType;


//! -------------------------- MsFEM error estimator ----------------------------
typedef MsFEMErrorEstimator< DiscreteFunctionType,
                             DiffusionType,
                             FirstSourceType,
                             MacroMicroGridSpecifierType,
                             SubGridListType > MsFEMErrorEstimatorType;
//! -----------------------------------------------------------------------------




//! ---------------------- important variables ---------------------------------

enum { dimension = GridType :: dimension};

//name of the error file in which the data will be saved
std :: string filename_;
std :: string path_;

int total_refinement_level_;
int coarse_grid_level_;

// only for uniform computations / use the same number of layers for every coarse grid entity
int number_of_layers_;

//! -----------------------------------------------------------------------------



//! ---------------------- local error indicators --------------------------------

// ----- local error indicators (for each coarse grid element T) -------------

int max_loop_number = 10;

bool local_indicators_available_ = false;

// local coarse residual, i.e. H ||f||_{L^2(T)}
std :: vector < std :: vector < RangeType > > loc_coarse_residual_( max_loop_number );

// local coarse grid jumps (contribute to the total coarse residual)
std :: vector < std :: vector < RangeType > > loc_coarse_grid_jumps_( max_loop_number );

// local projection error (we project to get a globaly continous approximation)
std :: vector < std :: vector < RangeType > > loc_projection_error_( max_loop_number );

// local jump in the conservative flux
std :: vector < std :: vector < RangeType > > loc_conservative_flux_jumps_( max_loop_number );

// local approximation error
std :: vector < std :: vector < RangeType > > loc_approximation_error_( max_loop_number );

// local sum over the fine grid jumps (for a fixed subgrid that cooresponds with a coarse entity T)
std :: vector < std :: vector < RangeType > > loc_fine_grid_jumps_( max_loop_number );


std :: vector < RangeType > total_coarse_residual_( max_loop_number );
std :: vector < RangeType > total_projection_error_( max_loop_number );
std :: vector < RangeType > total_coarse_grid_jumps_( max_loop_number );
std :: vector < RangeType > total_conservative_flux_jumps_( max_loop_number );
std :: vector < RangeType > total_approximation_error_( max_loop_number );
std :: vector < RangeType > total_fine_grid_jumps_( max_loop_number );

std :: vector < RangeType > total_estimated_H1_error_( max_loop_number );

bool repeat_algorithm_ = true;
int loop_number_ = 0;

#ifdef ADAPTIVE
double error_tolerance_;
#endif

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
      return 1;   // VTK
      //return 2; // VTK vertex data
      //return 3; // gnuplot
    }


};

//!---------------------------------------------------------------------------------------



void algorithm ( std :: string& macroGridName,
                 std :: ofstream& data_file )
{

  std :: cout << "loading dgf: " << macroGridName << std :: endl;

  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values for the parameters:

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer( macroGridName );
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine( coarse_grid_level_ );

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

  // coarse grid

  GridPointerType macro_grid_pointer_coarse( macroGridName );
  macro_grid_pointer_coarse->globalRefine( coarse_grid_level_ );

  GridPartType gridPart_coarse( *macro_grid_pointer_coarse);
  GridType &grid_coarse = gridPart_coarse.grid();


// strategy for adaptivity:
#ifdef ADAPTIVE

  typedef GridType :: LeafGridView GridView;
  typedef GridView :: Codim < 0 >:: Iterator ElementLeafIterator ;

  typedef GridType :: Traits :: LeafIndexSet GridLeafIndexSet;

  if ( local_indicators_available_ == true )
   {

      bool coarse_scale_error_dominant = false;
      bool fine_scale_error_dominant = false; // wird noch nicht benoetigt, dass wir diese Verfeinerung uniform regeln
      bool oversampling_error_dominant = false; // wird noch nicht benoetigt, dass wir diese Adadption uniform regeln

      // identify the dominant contribution:
      double average_est_error = total_estimated_H1_error_[ loop_number_ - 1 ] / 6.0; // 6 contributions


      if ( (total_approximation_error_[ loop_number_ - 1 ] >= average_est_error ) || 
           (total_fine_grid_jumps_[ loop_number_ - 1 ] >= average_est_error ) )
       {
         fine_scale_error_dominant = true; /* increase fine level resolution by 2 levels */
         total_refinement_level_ += 2; // 'the fine grid level'
         data_file << "Fine scale error identified as being dominant. Decrease the number of global refinements by 2." << std :: endl;
         data_file << "NEW: Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std :: endl;
         data_file << "Note that this means: the fine grid is " << total_refinement_level_ - coarse_grid_level_
                   << " refinement levels finer than the coarse grid." << std :: endl;
      }


      if ( (total_projection_error_[ loop_number_ - 1 ] >= average_est_error ) || 
           (total_conservative_flux_jumps_[ loop_number_ - 1 ] >= average_est_error ) )
       {
         oversampling_error_dominant = true; /* increase number of layers by 1 */
         number_of_layers_ += 1;
         data_file << "Oversampling error identified as being dominant. Increase the number of layers for each subgrid by 5." << std :: endl;
         data_file << "NEW: Number of layers = " << number_of_layers_ << std :: endl; 
       }


      if ( (total_coarse_residual_[ loop_number_ - 1 ] >= average_est_error ) || 
           (total_coarse_grid_jumps_[ loop_number_ - 1 ] >= average_est_error ) )
       {
         data_file << "Coarse scale error identified as being dominant. Start adaptive coarse grid refinment." << std :: endl;
         coarse_scale_error_dominant = true; /* mark elementwise for 2 refinments */
       }

      std :: vector < RangeType > average_coarse_error_indicator( loop_number_ ); //arithmetic average
      for ( int l = 0; l < loop_number_; ++l )
        {
          assert( loc_coarse_residual_[ l ].size() == loc_coarse_grid_jumps_[ l ].size() );
          average_coarse_error_indicator[ l ] = 0.0;
          for ( int i = 0; i < loc_coarse_grid_jumps_[ l ].size(); ++i )
            {
              average_coarse_error_indicator[ l ] += loc_coarse_residual_[ l ][ i ]  + loc_coarse_grid_jumps_[ l ][ i ];
            }
          average_coarse_error_indicator[ l ] = average_coarse_error_indicator[ l ] / loc_coarse_residual_[ l ].size();
        }

      // allgemeineren Algorithmus vorgestellt, aber noch nicht implementiert

      if ( coarse_scale_error_dominant == true ) {

      int number_of_refinements = 4;

      // allowed varianve from average ( in percent )
      double variance = 1.1; // = 110 % of the average

      for ( int l = 0; l < loop_number_; ++l )
        {

          const GridLeafIndexSet& gridLeafIndexSet = grid.leafIndexSet();
          GridView gridView = grid.leafView();

          int total_number_of_entities = 0;
          int number_of_marked_entities = 0;
          for ( ElementLeafIterator it = gridView.begin<0>();
                it != gridView.end<0>(); ++it )
            {
              RangeType loc_coarse_error_indicator = loc_coarse_grid_jumps_[ l ][ gridLeafIndexSet.index( *it ) ] + loc_coarse_residual_[ l ][ gridLeafIndexSet.index( *it ) ];
              total_number_of_entities += 1;

              if ( loc_coarse_error_indicator >= variance * average_coarse_error_indicator[ l ] )
                {
                  grid.mark( number_of_refinements , *it );
                  number_of_marked_entities += 1;
                }
            }

          if ( l == (loop_number_-1) )
           {
             data_file << number_of_marked_entities << " of " << total_number_of_entities << " coarse grid entities marked for mesh refinement." << std :: endl;
           }

          grid.preAdapt();
          grid.adapt();
          grid.postAdapt();

          // coarse grid

          GridView gridView_coarse = grid_coarse.leafView();
          const GridLeafIndexSet& gridLeafIndexSet_coarse = grid_coarse.leafIndexSet();

          for ( ElementLeafIterator it = gridView_coarse.begin<0>();
                it != gridView_coarse.end<0>(); ++it )
            {
              RangeType loc_coarse_error_indicator = loc_coarse_grid_jumps_[ l ][ gridLeafIndexSet_coarse.index( *it ) ] + loc_coarse_residual_[ l ][ gridLeafIndexSet_coarse.index( *it ) ];

              if ( loc_coarse_error_indicator >= variance * average_coarse_error_indicator[ l ] )
               { grid_coarse.mark( number_of_refinements , *it ); }
            }

          grid_coarse.preAdapt();
          grid_coarse.adapt();
          grid_coarse.postAdapt();
	}

      }

   }

#endif

  grid.globalRefine( total_refinement_level_ - coarse_grid_level_ );

  //! ------------------------- discrete function spaces -----------------------------------

  // the global-problem function space:
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  DiscreteFunctionSpaceType discreteFunctionSpace_coarse( gridPart_coarse );

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



  //! --------------- writing data output for the coarse grid visualization ------------------

  DiscreteFunctionType coarse_grid_visualization( filename_ + " Visualization of the coarse grid", discreteFunctionSpace_coarse );
  coarse_grid_visualization.clear();

  // -------------------------- data output -------------------------

  // create and initialize output class
  IOTupleType coarse_grid_series( &coarse_grid_visualization );



  char coarse_grid_fname[50];
  sprintf( coarse_grid_fname, "coarse_grid_visualization_%d_", loop_number_ );
  std :: string coarse_grid_fname_s( coarse_grid_fname );
  outputparam.set_prefix( coarse_grid_fname_s );
  DataOutputType coarse_grid_dataoutput( gridPart_coarse.grid(), coarse_grid_series, outputparam );
  // write data
  outstring << coarse_grid_fname_s;

  coarse_grid_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------

  //! --------------------------------------------------------------------------------------



  //! ---------------------- solve MsFEM problem ---------------------------  

  //! solution vector
  // solution of the standard finite element method
  DiscreteFunctionType msfem_solution( filename_ + " MsFEM Solution", discreteFunctionSpace );
  msfem_solution.clear();

  DiscreteFunctionType coarse_part_msfem_solution( filename_ + " Coarse Part MsFEM Solution", discreteFunctionSpace );
  coarse_part_msfem_solution.clear();

  DiscreteFunctionType fine_part_msfem_solution( filename_ + " Fine Part MsFEM Solution", discreteFunctionSpace );
  fine_part_msfem_solution.clear();

  int number_of_level_host_entities = grid_coarse.size( 0 /*codim*/ );  
  int coarse_level_fine_level_difference = grid.maxLevel() - grid_coarse.maxLevel();
  
  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MacroMicroGridSpecifierType specifier( discreteFunctionSpace_coarse, discreteFunctionSpace );
  for ( int i = 0; i < number_of_level_host_entities; i+=1 )
    { specifier.setLayer( i , number_of_layers_ ); }

  //! create subgrids:
  bool silence = false;
  SubGridListType subgrid_list( specifier, silence );

#if 1
  // just for Dirichlet zero-boundary condition
  Elliptic_MsFEM_Solver< DiscreteFunctionType > msfem_solver( discreteFunctionSpace, data_file, path_ );
  msfem_solver.solve_dirichlet_zero( diffusion_op, f, specifier, subgrid_list,
                                     coarse_part_msfem_solution, fine_part_msfem_solution, msfem_solution );
#endif




  std :: cout << "Data output for MsFEM Solution." << std :: endl;

  //! ----------------- writing data output MsFEM Solution -----------------

  // --------- VTK data output for MsFEM solution --------------------------

  // create and initialize output class
  IOTupleType msfem_solution_series( &msfem_solution );

#ifdef ADAPTIVE
  char msfem_fname[50];
  sprintf( msfem_fname, "msfem_solution_%d_", loop_number_ );
  std :: string msfem_fname_s( msfem_fname );
  outputparam.set_prefix( msfem_fname_s );
  DataOutputType msfem_dataoutput( gridPart.grid(), msfem_solution_series, outputparam );
  // write data
  outstring << msfem_fname_s;
#else
  outputparam.set_prefix("msfem_solution");
  DataOutputType msfem_dataoutput( gridPart.grid(), msfem_solution_series, outputparam );
  // write data
  outstring << "msfem_solution";
#endif
  
  msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());


  // create and initialize output class
  IOTupleType coarse_msfem_solution_series( &coarse_part_msfem_solution );
  
#ifdef ADAPTIVE
  char coarse_msfem_fname[50];
  sprintf( coarse_msfem_fname, "coarse_part_msfem_solution_%d_", loop_number_ );
  std :: string coarse_msfem_fname_s( coarse_msfem_fname );
  outputparam.set_prefix( coarse_msfem_fname_s );
  DataOutputType coarse_msfem_dataoutput( gridPart.grid(), coarse_msfem_solution_series, outputparam );
  // write data
  outstring << coarse_msfem_fname_s;
#else
  outputparam.set_prefix("coarse_part_msfem_solution");
  DataOutputType coarse_msfem_dataoutput( gridPart.grid(), coarse_msfem_solution_series, outputparam );
  // write data
  outstring << "coarse_part_msfem_solution";
#endif

  coarse_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());



  // create and initialize output class
  IOTupleType fine_msfem_solution_series( &fine_part_msfem_solution );
  
#ifdef ADAPTIVE
  char fine_msfem_fname[50];
  sprintf( fine_msfem_fname, "fine_part_msfem_solution_%d_", loop_number_ );
  std :: string fine_msfem_fname_s( fine_msfem_fname );
  outputparam.set_prefix( fine_msfem_fname_s );
  DataOutputType fine_msfem_dataoutput( gridPart.grid(), fine_msfem_solution_series, outputparam );
  // write data
  outstring << fine_msfem_fname_s;
#else
  outputparam.set_prefix("fine_part_msfem_solution");
  DataOutputType fine_msfem_dataoutput( gridPart.grid(), fine_msfem_solution_series, outputparam );
  // write data
  outstring << "fine_msfem_solution";
#endif
  
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




  //! ------------------------------------------------------------

  std :: cout << std :: endl << "The L2 errors:" << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "The L2 errors:" << std :: endl << std :: endl; }

  //! ----------------- compute L2- and H1- errors -------------------

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

  //! ----------------------------------------------------------------------



  //! ----------------------------------------------------------------------

// error estimation
#if 1

  RangeType total_estimated_H1_error( 0.0 );

  MsFEMErrorEstimatorType estimator( discreteFunctionSpace, specifier, subgrid_list, diffusion_op, f, path_ );
  total_estimated_H1_error = estimator.adaptive_refinement( grid_coarse, msfem_solution, coarse_part_msfem_solution, fine_part_msfem_solution, data_file );

  loc_coarse_residual_[ loop_number_ ].clear();
  loc_coarse_grid_jumps_[ loop_number_ ].clear();
  loc_projection_error_[ loop_number_ ].clear();
  loc_conservative_flux_jumps_[ loop_number_ ].clear();
  loc_approximation_error_[ loop_number_ ].clear();
  loc_fine_grid_jumps_[ loop_number_ ].clear();

  for ( int m = 0; m < specifier.getNumOfCoarseEntities(); ++m )
    {
      loc_coarse_residual_[ loop_number_ ].push_back( 0.0 );
      loc_coarse_grid_jumps_[ loop_number_ ].push_back( 0.0 );
      loc_projection_error_[ loop_number_ ].push_back( 0.0 );
      loc_conservative_flux_jumps_[ loop_number_ ].push_back( 0.0 );
      loc_approximation_error_[ loop_number_ ].push_back( 0.0 );
      loc_fine_grid_jumps_[ loop_number_ ].push_back( 0.0 );
    }

  total_coarse_residual_[ loop_number_ ] = 0.0;
  total_projection_error_[ loop_number_ ] = 0.0;
  total_coarse_grid_jumps_[ loop_number_ ] = 0.0;
  total_conservative_flux_jumps_[ loop_number_ ] = 0.0;
  total_approximation_error_[ loop_number_ ] = 0.0;
  total_fine_grid_jumps_[ loop_number_ ] = 0.0;

  for ( int m = 0; m < specifier.getNumOfCoarseEntities(); ++m )
    {
      loc_coarse_residual_[ loop_number_ ][ m ] = specifier.get_loc_coarse_residual( m );
      loc_coarse_grid_jumps_[ loop_number_ ][ m ] = specifier.get_loc_coarse_grid_jumps( m );
      loc_projection_error_[ loop_number_ ][ m ] = specifier.get_loc_projection_error( m );
      loc_conservative_flux_jumps_[ loop_number_ ][ m ] = specifier.get_loc_conservative_flux_jumps( m );
      loc_approximation_error_[ loop_number_ ][ m ] = specifier.get_loc_approximation_error( m );
      loc_fine_grid_jumps_[ loop_number_ ][ m ] = specifier.get_loc_fine_grid_jumps( m );

      total_coarse_residual_[ loop_number_ ] += pow( loc_coarse_residual_[ loop_number_ ][ m ] , 2.0 );
      total_projection_error_[ loop_number_ ] += pow( loc_projection_error_[ loop_number_ ][ m ] , 2.0 );
      total_coarse_grid_jumps_[ loop_number_ ] += pow( loc_coarse_grid_jumps_[ loop_number_ ][ m ] , 2.0 );
      total_conservative_flux_jumps_[ loop_number_ ] += pow( loc_conservative_flux_jumps_[ loop_number_ ][ m ] , 2.0 );
      total_approximation_error_[ loop_number_ ] += pow( loc_approximation_error_[ loop_number_ ][ m ] , 2.0 );
      total_fine_grid_jumps_[ loop_number_ ] += pow( loc_fine_grid_jumps_[ loop_number_ ][ m ] , 2.0 );
    }

  total_coarse_residual_[ loop_number_ ] = sqrt( total_coarse_residual_[ loop_number_ ] );
  total_projection_error_[ loop_number_ ] = sqrt( total_projection_error_[ loop_number_ ] );
  total_coarse_grid_jumps_[ loop_number_ ] = sqrt( total_coarse_grid_jumps_[ loop_number_ ] );
  total_conservative_flux_jumps_[ loop_number_ ] = sqrt( total_conservative_flux_jumps_[ loop_number_ ] );
  total_approximation_error_[ loop_number_ ] = sqrt( total_approximation_error_[ loop_number_ ] );
  total_fine_grid_jumps_[ loop_number_ ] = sqrt( total_fine_grid_jumps_[ loop_number_ ] );

  total_estimated_H1_error_[ loop_number_ ] = total_coarse_residual_[ loop_number_ ];
  total_estimated_H1_error_[ loop_number_ ] += total_projection_error_[ loop_number_ ];
  total_estimated_H1_error_[ loop_number_ ] += total_coarse_grid_jumps_[ loop_number_ ];
  total_estimated_H1_error_[ loop_number_ ] += total_conservative_flux_jumps_[ loop_number_ ];
  total_estimated_H1_error_[ loop_number_ ] += total_approximation_error_[ loop_number_ ];
  total_estimated_H1_error_[ loop_number_ ] += total_fine_grid_jumps_[ loop_number_ ];

  local_indicators_available_ = true;

  #ifdef ADAPTIVE
  if ( total_estimated_H1_error <= error_tolerance_ )
    repeat_algorithm_ = false;
  #endif

  #ifndef ADAPTIVE
  repeat_algorithm_ = false;
  #endif

#endif


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


#ifndef EXACTSOLUTION_AVAILABLE

  std :: cout << "Exact solution not available. Use fine-scale FEM-approximation as a reference solution." << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "Exact solution not available. Use fine-scale FEM-approximation as a reference solution." << std :: endl << std :: endl; }


  RangeType approx_msfem_error = l2error.norm2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( fem_solution, msfem_solution );

  std :: cout << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std :: endl; }


  H1Norm< GridPartType > h1norm( gridPart );
  RangeType h1_approx_msfem_error = h1norm.distance( fem_solution, msfem_solution );

  std :: cout << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std :: endl; }

#endif

//! -------------------------------------------------------




//! irgendwann loeschen!
// dirty hack for oversampling error visualization!
#if 0

     DiscreteFunctionType oversampling_error_visualization( filename_ + " Visualization of the oversampling error", discreteFunctionSpace_coarse );
     oversampling_error_visualization.clear();

     typedef DiscreteFunctionSpaceType :: IteratorType :: Entity EntityType;

     typedef EntityType :: EntityPointer EntityPointerType;
     typedef EntityType :: Codim< 2 > :: EntityPointer NodePointerType;
     
     typedef DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
     
     int number_of_coarse_nodes = gridPart_coarse.grid().size( 2 /*codim*/ );
     std :: vector< std :: vector < EntityPointerType > > entities_sharing_same_node( number_of_coarse_nodes );

     GridView gridView_coarse_new = grid_coarse.leafView();
     const GridLeafIndexSet& gridLeafIndexSet_coarse_new = grid_coarse.leafIndexSet();

     for ( ElementLeafIterator it = gridView_coarse_new.begin<0>();
           it != gridView_coarse_new.end<0>(); ++it )
        {
          int number_of_nodes_in_entity = (*it).count<2>();
          for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const NodePointerType node = (*it).subEntity<2>(i);
	      int global_index_node = gridPart_coarse.indexSet().index( *node );

	      entities_sharing_same_node[ global_index_node ].push_back( it );
	    }
        }

     RangeType average_conservative_flux_jumps = 0.0;

     std :: cout << "total_conservative_flux_jumps_[ loop_number_ ] = " << total_conservative_flux_jumps_[ loop_number_ ] << std :: endl;
     std :: cout << "average_conservative_flux_jumps = " << average_conservative_flux_jumps << std :: endl << std :: endl;


     double cum = 0;
     for ( ElementLeafIterator it = gridView_coarse_new.begin<0>();
           it != gridView_coarse_new.end<0>(); ++it )
        {
           average_conservative_flux_jumps += loc_conservative_flux_jumps_[ loop_number_ ][ gridLeafIndexSet_coarse_new.index( *it ) ];
           cum += 1.0;
        }
     average_conservative_flux_jumps /= cum;
        

     for ( ElementLeafIterator it = gridView_coarse_new.begin<0>();
           it != gridView_coarse_new.end<0>(); ++it )
        {

           int number_of_nodes = (*it).count<2>();

           LocalFunctionType loc_func = oversampling_error_visualization.localFunction( *it );

           RangeType loc_oversampling_indicator = loc_conservative_flux_jumps_[ loop_number_ ][ gridLeafIndexSet_coarse_new.index( *it ) ];
           std :: cout << "loc_oversampling_indicator = " << loc_oversampling_indicator << std :: endl;

           for ( int i = 0; i < number_of_nodes; i += 1 )
             {
                const NodePointerType node = (*it).subEntity<2>(i);
                int global_index_node = gridPart_coarse.indexSet().index( *node );

                if ( loc_oversampling_indicator <= 0.2 * average_conservative_flux_jumps )
                 { loc_func[ i ] += ( 0.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 0.2 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 0.4 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 1.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 0.4 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 0.6 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 2.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 0.6 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 0.8 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 3.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 0.8 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 1.0 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 4.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 1.0 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 1.2 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 5.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 1.2 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 1.4 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 6.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 1.4 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 1.6 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 7.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 1.6 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 1.8 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 8.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

              if ( ( loc_oversampling_indicator > 1.8 * average_conservative_flux_jumps) &&
                   ( loc_oversampling_indicator <= 200.0 * average_conservative_flux_jumps) )
                 { loc_func[ i ] += ( 9.0 + number_of_layers_ ) / entities_sharing_same_node[ global_index_node ].size(); }

             }

        }

  // -------------------------- data output -------------------------

  // create and initialize output class
  IOTupleType oversampling_coarse_grid_series( &oversampling_error_visualization );

  char oversampling_coarse_grid_fname[50];
  sprintf( oversampling_coarse_grid_fname, "oversampling_error_visualization_%d_", loop_number_ );
  std :: string oversampling_coarse_grid_fname_s( oversampling_coarse_grid_fname );
  outputparam.set_prefix( oversampling_coarse_grid_fname_s );
  DataOutputType oversampling_coarse_grid_dataoutput( gridPart_coarse.grid(), oversampling_coarse_grid_series, outputparam );
  // write data
  outstring << oversampling_coarse_grid_fname_s;

  oversampling_coarse_grid_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------

#endif

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

  path_ = "data/MsFEM/" + filename_;

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

  std :: string save_filename = path_ + "/problem-info.txt";
  std :: cout << "Data will be saved under: " << save_filename << std :: endl;

  std :: cout << "Enter coarse grid level: ";
  std :: cin >> coarse_grid_level_;
  if ( coarse_grid_level_ > atoi( argv[ 1 ] ))
   {
     std :: cout << "Coarse grid level must be smaller than the starting refinment level." << std :: endl;
     abort();
   }

#ifdef UNIFORM
  std :: cout << "Enter number of layers for oversampling: ";
  std :: cin >> number_of_layers_;
#else
  std :: cout << "Enter initial number of layers for oversampling: ";
  std :: cin >> number_of_layers_;  
#endif

#ifdef ADAPTIVE
  std :: cout << "Enter error tolerance: ";
  std :: cin >> error_tolerance_;
#endif

  // data for the model problem; the information manager
  // (see 'problem_specification.hh' for details)
  Problem::ModelProblemData info( filename_ );

  // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
  total_refinement_level_ = atoi( argv[ 1 ] );

  //name of the grid file that describes the macro-grid:
  std :: string macroGridName;
  info.getMacroGridFile( macroGridName );

  // to save all information in a file
  std :: ofstream data_file( (save_filename).c_str() );
  if (data_file.is_open())
            {
               data_file << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << " with epsilon = " << info.getEpsilon() << "." << std :: endl << std :: endl;
#ifdef UNIFORM
               data_file << "Use MsFEM with an uniform computation, i.e.:" << std :: endl;
               data_file << "Uniformly refined coarse and fine mesh and" << std :: endl;
               data_file << "the same number of layers for each (oversampled) local grid computation." << std :: endl << std :: endl;
               data_file << "Computations were made for:" << std :: endl << std :: endl;
               data_file << "Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std :: endl;
               data_file << "Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std :: endl;
               data_file << "Number of layers for oversampling = " << number_of_layers_ << std :: endl;
               data_file << std :: endl << std :: endl;
#else
               data_file << "Use MsFEM with an adaptive computation, i.e.:" << std :: endl;
               data_file << "Starting with a uniformly refined coarse and fine mesh and" << std :: endl;
               data_file << "the same number of layers for each (oversampled) local grid computation." << std :: endl << std :: endl;
               data_file << "Error tolerance = " <<  error_tolerance_ << std :: endl << std :: endl;
               data_file << "Computations were made for:" << std :: endl << std :: endl;
               data_file << "(Starting) Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std :: endl;
               data_file << "(Starting) Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std :: endl;
               data_file << "(Starting) Number of layers for oversampling = " << number_of_layers_ << std :: endl;
               data_file << std :: endl << std :: endl;
#endif
            }

    loop_number_ = 0;
    while ( repeat_algorithm_ == true )
     {
#ifdef ADAPTIVE
       data_file << "------------------ run " << loop_number_ + 1<< " --------------------" << std :: endl << std :: endl;
#endif
       algorithm( macroGridName , data_file );
#ifdef ADAPTIVE
       data_file << std :: endl << std :: endl;
       data_file << "---------------------------------------------" << std :: endl;
       loop_number_ += 1;
#endif
     }
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
