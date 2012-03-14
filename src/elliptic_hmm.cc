
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


//! do we have a linear elliptic problem?
// if yes, #define LINEAR_PROBLEM
// #define LINEAR_PROBLEM

//! TFR-HMM or simple HMM?
//#define TFR

//! is an exact solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro EXACTSOLUTION_AVAILABLE

//! is the homogenized solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro HOMOGENIZEDSOL_AVAILABLE
// (if HOMOGENIZEDSOL_AVAILABLE == true, it means that it can be computed. It still needs to be determined by using a homogenizer )
// #define HOMOGENIZEDSOL_AVAILABLE

//! Do we solve the cell problems ad hoc or in a pre-process?
// the second possibility requires a data file where the solutions are saved (file becomes large)
//#define AD_HOC_COMPUTATION

//! Do we have/want a fine-scale reference solution?
//#define FINE_SCALE_REFERENCE
#ifdef FINE_SCALE_REFERENCE

  // load the precomputed fine scale reference from a file
  #define FSR_LOAD

  #ifndef FSR_LOAD
  // compute the fine scale reference (on the fly)
    #define FSR_COMPUTE

    #ifdef FSR_COMPUTE
     // Do we write the discrete fine-scale solution to a file? (for later usage)
      #define WRITE_FINESCALE_SOL_TO_FILE
    #endif

  #endif

#endif


//! Do we have a HMM reference solution? (precomputed detailed HMM simulation)
// we might use a detailed HMM computation as a reference! (if it is available)
//#define HMM_REFERENCE

//! Do we write the discrete HMM solution to a file? (for later usage)
#define WRITE_HMM_SOL_TO_FILE

//! Do we want to use error estimation (a-posteriori estimate and adaptivity)?
// Not possible for ad-hoc computations! (in this case, error estimation is far too expensive)
#ifndef AD_HOC_COMPUTATION
  //
  //#define ERRORESTIMATION
  // only possible if we use error estimation:
  #ifdef ERRORESTIMATION
   // Do you want to allow adaptive mesh refinement?
   //#define ADAPTIVE
  #endif

#endif

//! Do we want to add a stochastic perturbation on the data?
//#define STOCHASTIC_PERTURBATION
#ifdef STOCHASTIC_PERTURBATION
 // size of variance:
 #define VARIANCE 0.01
 //! Do we want to force the algorithm to come to an end?
 // (was auch immer der Grund war, dass das Programm zuvor endlos lange weiter gelaufen ist. z.B. Tolerenzen nicht erreicht etc.)
 #define FORCE
#endif

//! if a computation was broken (after a certain HMM Newton step), we might want to resume to this computation,
//! loading the solution of the last step that was succesfully carried out (it has to be saved somewhere!)
//  (this only works for non-adaptive computations!)
//#define RESUME_TO_BROKEN_COMPUTATION

#ifdef RESUME_TO_BROKEN_COMPUTATION
 // last HMM Newton step that was succesfully carried out, saving the iterate afterwards
 #define HMM_NEWTON_ITERATION_STEP 2
#else
 // default: we need a full computation. start with step 1:
 #define HMM_NEWTON_ITERATION_STEP 0 
#endif


//!NOTE: All the multiscale code requires an access to the 'ModelProblemData' class (typically defined in problem_specification.hh), which provides us with information about epsilon, delta, etc.
// HMM Assembler, Error Estimator, ... they all hark back to 'ModelProblemData'. Probably there is a better solution, but for me, it works perfectly.

namespace Multiscale
{


  // parameters for the current realization of the HMM
  class HMMParameters
  {


  public:
    // Constructor for ModelProblemData
    inline explicit HMMParameters ( )
    {
    }

    // do you want to save the solutions of the cell problems for a later usage?
    // (e.g. for error estimation and adaptivity)
    inline bool save_cell_problems ( ) const
    {
      return true;
    }

  };

}



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
#include <dune/fem/gridpart/periodicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>

#if 1
//hmfemmain:
//!---------
#include <dune/fem/solver/oemsolver/oemsolver.hh>


#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
//#include <dune/fem/io/visual/grape/datadisp/errordisplay.hh>
//!-----------
#endif
#if 0
//poisson
//!--------

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/visual/grape/datadisp/errordisplay.hh>
//!---------
#endif
#include <dune/fem/solver/inverseoperators.hh>




//! local (dune-multiscale) includes
#include <dune/multiscale/problems/elliptic_problems/model_problem_8/problem_specification.hh>


#include <dune/multiscale/operators/righthandside_assembler.hh>

#include <dune/multiscale/operators/disc_func_writer/discretefunctionwriter.hh>

#include <dune/multiscale/operators/cell_problem_solving/cellproblemsolver.hh>

#include <dune/multiscale/operators/reconstruction_manager/elliptic/reconstructionintegrator.hh>

// we only use error estimation, if the solutions of the cell problems have been determined in a pre-process. Otherwise it is far too expensive!
#ifndef AD_HOC_COMPUTATION
 #include <dune/multiscale/operators/errorestimation/elliptic_error_estimator.hh>
#endif

#include <dune/multiscale/operators/meanvalue.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/multiscale/operators/matrix_assembler/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/operators/matrix_assembler/elliptic_hmm_matrix_assembler.hh>

//! (very restrictive) homogenizer
#ifdef LINEAR_PROBLEM
#include <dune/multiscale/operators/homogenizer/elliptic_analytical_homogenizer.hh>
#include <dune/multiscale/operators/homogenizer/elliptic_homogenizer.hh>
#else
// dummy (does not work, since identical to HMM assembler)
#include <dune/multiscale/operators/homogenizer/nonlinear_elliptic_homogenizer.hh>
#endif

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





//! --------- typedefs for the coefficient and data functions ------------------------------

// type of first source term (right hand side of differential equation or type of 'f')
typedef Problem::FirstSource< FunctionSpaceType > FirstSourceType;

// type of second source term 'G' (second right hand side of differential equation 'div G')
typedef Problem::SecondSource< FunctionSpaceType > SecondSourceType;

// type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
typedef Problem::Diffusion< FunctionSpaceType > DiffusionType;

// type of mass (or reaction) term (i.e. 'm' or 'c')
typedef Problem::MassTerm< FunctionSpaceType > MassTermType;

// default type for any missing coefficient function (e.g. advection,...)
typedef Problem::DefaultDummyFunction< FunctionSpaceType > DefaultDummyFunctionType;

#ifdef EXACTSOLUTION_AVAILABLE
// type of exact solution (in general unknown)
typedef Problem::ExactSolution< FunctionSpaceType > ExactSolutionType;
typedef DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
  DiscreteExactSolutionType; //for data output with paraview or grape
#endif

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





//! --------- typedefs for the periodic micro grid and the corresponding discrete space ----

typedef PeriodicLeafGridPart< GridType > PeriodicGridPartType;

typedef LagrangeDiscreteFunctionSpace
        < FunctionSpaceType, PeriodicGridPartType, 1 > // 1 =POLORDER
  PeriodicDiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction < PeriodicDiscreteFunctionSpaceType > PeriodicDiscreteFunctionType;

//!-----------------------------------------------------------------------------------------



//! ------------ cell problem solver and numbering manager -----------------------------------------

typedef CellProblemNumberingManager< DiscreteFunctionSpaceType > CellProblemNumberingManagerType;

//!-----------------------------------------------------------------------------------------


//! --------------------- the standard matrix traits -------------------------------------

struct MatrixTraits
{
  typedef DiscreteFunctionSpaceType RowSpaceType;
  typedef DiscreteFunctionSpaceType ColumnSpaceType;
  typedef LagrangeMatrixSetup< false > StencilType;
  typedef ParallelScalarProduct< DiscreteFunctionSpaceType > ParallelScalarProductType;

  template< class M >
  struct Adapter
  {
    typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
  };
};

//! --------------------------------------------------------------------------------------


//! --------------------- type of fem stiffness matrix -----------------------------------

typedef SparseRowMatrixOperator< DiscreteFunctionType, DiscreteFunctionType, MatrixTraits > FEMMatrix;

//! --------------------------------------------------------------------------------------


//! --------------- solver for the linear system of equations ----------------------------

// use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric ones. GMRES seems to be more stable, but is extremely slow!
typedef OEMBICGSQOp/*OEMBICGSTABOp*/< DiscreteFunctionType, FEMMatrix > InverseFEMMatrix;

//! --------------------------------------------------------------------------------------


//! --------------- the discrete operators (standard FEM and HMM) ------------------------

// discrete elliptic operator (corresponds with FEM Matrix)
typedef DiscreteEllipticOperator< DiscreteFunctionType, DiffusionType, MassTermType > EllipticOperatorType;

// discrete elliptic HMM operator (corresponds with HMM (or HMFEM) Matrix)
typedef DiscreteEllipticHMMOperator< DiscreteFunctionType, PeriodicDiscreteFunctionType, DiffusionType, CellProblemNumberingManagerType > EllipticHMMOperatorType;

//! --------------------------------------------------------------------------------------




//! --------------- ERROR ESTIMATOR NOT YET IMPLEMENTED ------------------------
typedef ErrorEstimator< PeriodicDiscreteFunctionType, 
                        DiscreteFunctionType,
                        DiffusionType > ErrorEstimatorType;
//! -----------------------------------------------------------------------------





//! -------------------------------- important variables ---------------------------------

enum { dimension = GridType :: dimension};

//name of the error file in which the data will be saved
std :: string filename_;

double epsilon_; // 'epsilon' in for instance A^{epsilon}(t,x) = A(t,x/epsilon)
double epsilon_est_; // estimated epsilon in case epsilon is unknown
double delta_; //edge length of the cells in the cell proplems,

int refinement_level_macrogrid_;
int refinement_level_referenceprob_;

double error_tolerance_;

//std :: ofstream data_file_; // file where we save the data

//! -----------------------------------------------------------------------------



//! --------- typedefs and classes for data output -----------------------------------------

typedef Tuple<DiscreteFunctionType*> IOTupleType;
typedef DataOutput<GridType, IOTupleType> DataOutputType;

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
        return "data_output_hmm";
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






//!------------------------- for adaptive grid refinement ---------------------------------

//! For adaption:

//! type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunctionType >
  RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< GridType, RestrictProlongOperatorType >
  AdaptationManagerType;

//!---------------------------------------------------------------------------------------



//! set the dirichlet points to zero
template< class EntityType, class DiscreteFunctionType >
void boundaryTreatment( const EntityType &entity, DiscreteFunctionType &rhs )
{
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
    LagrangePointSetType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  
  enum { faceCodim = 1 };

  typedef typename GridPartType :: IntersectionIteratorType
    IntersectionIteratorType;

  typedef typename LagrangePointSetType :: template Codim< faceCodim > 
                                        :: SubEntityIteratorType
    FaceDofIteratorType;

  const DiscreteFunctionSpaceType &discreteFunctionSpace = rhs.space();

  const GridPartType &gridPart = discreteFunctionSpace.gridPart();

  IntersectionIteratorType it = gridPart.ibegin( entity );
  const IntersectionIteratorType endit = gridPart.iend( entity );
  for( ; it != endit; ++it ) {
    if( !(*it).boundary() )
      continue;

    LocalFunctionType rhsLocal = rhs.localFunction( entity );
    const LagrangePointSetType &lagrangePointSet
      = discreteFunctionSpace.lagrangePointSet( entity );

    const int face = (*it).indexInInside();
    FaceDofIteratorType faceIterator
      = lagrangePointSet.template beginSubEntity< faceCodim >( face );
    const FaceDofIteratorType faceEndIterator
      = lagrangePointSet.template endSubEntity< faceCodim >( face );
    for( ; faceIterator != faceEndIterator; ++faceIterator )
      rhsLocal[ *faceIterator ] = 0;
  }
}


template < class Stream, class DiscFunc >
void oneLinePrint( Stream& stream, const DiscFunc& func )
{
    typedef typename DiscFunc::ConstDofIteratorType
        DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it )
        stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
}


RangeType get_size_of_domain( DiscreteFunctionSpaceType& discreteFunctionSpace )
{

   typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
   IteratorType endit = discreteFunctionSpace.end();
   RangeType size_of_domain = 0.0;

   for(IteratorType it = discreteFunctionSpace.begin(); it != endit ; ++it)
        {
          CachingQuadrature <GridPartType , 0 > entityQuadrature (*it, 0 );

          // get geoemetry of entity
          const  GridType :: Codim< 0 > :: Geometry& geometry = it->geometry();

          const double volumeEntity = entityQuadrature.weight( 0 ) * 
               geometry.integrationElement(entityQuadrature.point( 0 ));

          size_of_domain += volumeEntity;
        }

   return size_of_domain;

}


void algorithm ( std :: string &UnitCubeName,
                 GridPointerType &macro_grid_pointer, // grid pointer that belongs to the macro grid
                 GridPointerType &fine_macro_grid_pointer, // grid pointer that belongs to the fine macro grid (for reference computations)
                 GridPointerType &periodic_grid_pointer, // grid pointer that belongs to the periodic micro grid
                 int refinement_difference, //refinement difference for the macro grid (problem-to-solve vs. reference problem)
                 std :: ofstream &data_file )
{


  //! ---- tools ----

  // model problem data
  Problem::ModelProblemData problem_info;

  // set of hmm parameters/information
  Multiscale::HMMParameters method_info;

  L2Error< DiscreteFunctionType > l2error;

  // expensive hack to deal with discrete functions, defined on different grids
  ImprovedL2Error< DiscreteFunctionType > impL2error;

  //! ---------------------------- grid parts ----------------------------------------------

  // grid part for the global function space, required for HMM-macro-problem
  GridPartType gridPart( *macro_grid_pointer);

  // grid part for the periodic function space, required for HMM-cell-problems
  PeriodicGridPartType periodicGridPart ( *periodic_grid_pointer );

  // auxiliary grid part for the periodic function space, required for HMM-cell-problems
  GridPartType         auxiliaryGridPart( *periodic_grid_pointer );
  // auxiliaryGridPart for the error estimator (the auxiliaryGridPart yields an intersection iterator, which can not be done by the periodicGridPart)

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

  // the local-problem function space (containing periodic functions):
  PeriodicDiscreteFunctionSpaceType periodicDiscreteFunctionSpace( periodicGridPart );
  // and the corresponding auxiliary one:
  DiscreteFunctionSpaceType auxiliaryDiscreteFunctionSpace( auxiliaryGridPart );

  //! --------------------------------------------------------------------------------------



 //! --------------------------------------------------------------------------------------
//sollte bald in eigens Programm ausgelagert werden:
// (hier werden bereits berechnete diskrete HMM solutions eingelesen und die L^2-Differenz berechnet
// Leben tun alle diese Funktionen auf dem Makrogitter mit 10 Verfeinerungsleveln, wenn sie auf einem groeberen Gitter bestimmt worden sind, dann wurden sie spaeter darauf projeziert
#if 0

  //name of the grid file that describes the macro-grid:
  std :: string macroGridName;
  problem_info.getMacroGridFile( macroGridName );
  std :: cout << "loading dgf: " << macroGridName << std :: endl;


  std :: string discFunc_location_1;
  std :: string discFunc_location_2;


//    discFunc_location_1 = "data/HMM/Model_Problem_1/Macro_6_Micro_2/hmm_solution_discFunc_refLevel_6";
//  discFunc_location_1 = "data/HMM/Model_Problem_1/Macro_4_Micro_4/hmm_solution_discFunc_refLevel_4";
//  discFunc_location_1 = "data/HMM/Model_Problem_1/Macro_10_Micro_8_tolerance3.5e-06/hmm_solution_discFunc_refLevel_10";
//  discFunc_location_1 = "data/HMM/Model_Problem_1/Macro_6_Micro_6/hmm_solution_discFunc_refLevel_6";
//  discFunc_location_1 = "data/HMM/Model_Problem_1/Macro_8_Micro_8/hmm_solution_discFunc_refLevel_8";
//  discFunc_location_1 = "data/HMM/Model_Problem_2/reference_solution_ref_16/finescale_solution_discFunc_refLevel_16";
//  discFunc_location_1 = "data/HMM/Model_Problem_1/Macro_10_Micro_8/hmm_solution_discFunc_refLevel_10";
//  discFunc_location_1 = "data/HMM/Model_Problem_2/Macro_8_Micro_10_OVERSAMPLING/hmm_solution_discFunc_refLevel_8";
  discFunc_location_1 = "data/HMM/Model_Problem_2/Macro_6_Micro_6_tol_1e-08/hmm_solution_discFunc_refLevel_6";
//  discFunc_location_1 = "data/HMM/Model_Problem_2/Macro_4_Micro_10_STRANGE_OVERSAMPLING_TFR/hmm_solution_discFunc_refLevel_4";


//  discFunc_location_2 = "data/HMM/Model_Problem_1/Macro_8_Micro_8/hmm_solution_discFunc_refLevel_8";
//  discFunc_location_2 = "data/HMM/Model_Problem_1/Macro_10_Micro_8/hmm_solution_discFunc_refLevel_10";
//  discFunc_location_2 = "data/HMM/Model_Problem_1/reference_solution_ref_18/finescale_solution_discFunc_refLevel_18";
//  discFunc_location_2 = "data/HMM/Model_Problem_2/Macro_4_Micro_6_OVERSAMPLING/hmm_solution_discFunc_refLevel_4";
  discFunc_location_2 = "data/HMM/Model_Problem_2/reference_solution_ref_18/finescale_solution_discFunc_refLevel_18";
//  discFunc_location_2 = "data/HMM/Model_Problem_2/Macro_2_Micro_8_STRANGE_OVERSAMPLING/hmm_solution_discFunc_refLevel_2";
//  discFunc_location_2 = "data/HMM/Model_Problem_2/Macro_4_Micro_4/hmm_solution_discFunc_refLevel_4";
//  discFunc_location_2 = "data/HMMModel_Problem_2/zzz_inProgress/done/DELTA_0.1_EPSILON_0.05/Macro_8_Micro_10/hmm_solution_discFunc_refLevel_8";

  int gridLevel_1 = 6; // Macro_'gridLevel_1'...
  int gridLevel_2 = 18; // Macro_'gridLevel_2'...
//! Note: gridLevel_2 >= gridLevel_1

  std :: cout << "gridLevel_1 = " << gridLevel_1 << std :: endl;
  std :: cout << "gridLevel_2 = " << gridLevel_2 << std :: endl << std :: endl;

  std :: cout << "discFunc_location_1 = " << discFunc_location_1 << std :: endl;
  std :: cout << "discFunc_location_2 = " << discFunc_location_2 << std :: endl << std :: endl;

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer_1( macroGridName );
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer_1->globalRefine( gridLevel_1 );

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer_2( macroGridName );
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer_2->globalRefine( gridLevel_2 );

  GridPartType gridPart_1( *macro_grid_pointer_1);
  GridPartType gridPart_2( *macro_grid_pointer_2);

  GridType &grid_1 = gridPart_1.grid();
  GridType &grid_2 = gridPart_2.grid();

  std :: cout << "Grid 1 size vorher = " << grid_1.size(0) << std :: endl;
  std :: cout << "Grid 2 size vorher = " << grid_2.size(0) << std :: endl;

  DiscreteFunctionSpaceType discreteFunctionSpace_1( gridPart_1 );
  DiscreteFunctionSpaceType discreteFunctionSpace_2( gridPart_2 );

  DiscreteFunctionType discrete_function_1( " discrete_function_1 ", discreteFunctionSpace_1 );
  discrete_function_1.clear();
  DiscreteFunctionType discrete_function_2( " discrete_function_2 ", discreteFunctionSpace_2 );
  discrete_function_2.clear();

  DiscreteFunctionType zero_function_1( " zero_function_1 ", discreteFunctionSpace_1 );
  zero_function_1.clear();
  DiscreteFunctionType zero_function_2( " zero_function_2 ", discreteFunctionSpace_1 );
  zero_function_2.clear();
//L2 norm = 0.642166
//L2 norm = 0.623333
//L2 norm = 0.72218
//L2 norm = 0.72218
//L2 norm = 0.719281
//L2 norm = 0.792442


  bool reader_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_1( (discFunc_location_1).c_str() );
  discrete_function_reader_1.open();

  DiscreteFunctionReader discrete_function_reader_2( (discFunc_location_2).c_str() );
  discrete_function_reader_2.open();

  discrete_function_reader_1.read( 0, discrete_function_1 );
  std :: cout << "discrete_function_1 read." << std :: endl;


  discrete_function_reader_2.read( 0, discrete_function_2 );
  std :: cout << "discrete_function_2 read." << std :: endl;

#if 1

  //!warum wird das gebraucht, um das richtige Ergebnis zu bekommen???????
  #if 1
  L2Norm < GridPartType > norm_L2(gridPart_1);
  RangeType value = 0;
  value = norm_L2.norm(zero_function_1);
  std :: cout << "Value = " << value << std :: endl;
  #endif

#if 0
  // general output parameters
  myDataOutputParameters outputparam_1_vorher;
  outputparam_1_vorher.set_path( "data/HMM/" );

  // sequence stamp
  std::stringstream outstring_1_vorher;

  IOTupleType discrete_function_1_series_vorher( &zero_function_1/*discrete_function_1*/ );
  outputparam_1_vorher.set_prefix("discrete_function_1_vorher");
  DataOutputType discrete_function_1_dataoutput_vorher( grid_1, discrete_function_1_series_vorher, outputparam_1_vorher );

  // write data
  outstring_1_vorher << "discrete_function_1_vorher";
  discrete_function_1_dataoutput_vorher.writeData( 1.0 /*dummy*/, outstring_1_vorher.str() );
  // clear the std::stringstream:
  outstring_1_vorher.str(std::string());

  // -------------------------------------------------------
#endif

#if 0
  std :: cout << "Starting adaption 1...";

  // one for the discreteFunctionSpace
  RestrictProlongOperatorType rp_1( discrete_function_1 );
  AdaptationManagerType adaptationManager_1( grid_1, rp_1 );

  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType endit_1 = discreteFunctionSpace_1.end();
  for( IteratorType it = discreteFunctionSpace_1.begin(); it != endit_1; ++it )
    { grid_1.mark( (gridLevel_2-gridLevel_1) , *it ); }

  adaptationManager_1.adapt();

  std :: cout << " done." << std :: endl;
#endif

#if 0
  std :: cout << "Starting adaption 2...";

  // one for the discreteFunctionSpace
  RestrictProlongOperatorType rp_2( discrete_function_2 );
  AdaptationManagerType adaptationManager_2( gridPart_2.grid(), rp_2 );

  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType endit_2 = discreteFunctionSpace_2.end();
  for( IteratorType it = discreteFunctionSpace_2.begin(); it != endit_2; ++it )
    { gridPart_2.grid().mark( 0 , *it ); }
  adaptationManager_2.adapt();

  std :: cout << " done." << std :: endl;
#endif

  // --------- data output --------------

#if 0
  // general output parameters
  myDataOutputParameters outputparam_1;
  outputparam_1.set_path( "data/HMM/" );

  // sequence stamp
  std::stringstream outstring_1;

  IOTupleType discrete_function_1_series( &discrete_function_1 );
  outputparam_1.set_prefix("discrete_function_1");
  DataOutputType discrete_function_1_dataoutput( grid_1, discrete_function_1_series, outputparam_1 );

  // write data
  outstring_1 << "discrete_function_1";
  discrete_function_1_dataoutput.writeData( 1.0 /*dummy*/, outstring_1.str() );
  // clear the std::stringstream:
  outstring_1.str(std::string());




  myDataOutputParameters outputparam_2;
  outputparam_2.set_path( "data/HMM/" );

  // sequence stamp
  std::stringstream outstring_2;

  IOTupleType discrete_function_2_series( &discrete_function_2 );
  outputparam_2.set_prefix("discrete_function_2");
  DataOutputType discrete_function_2_dataoutput( grid_2, discrete_function_2_series, outputparam_2 );

  // write data
  outstring_2 << "discrete_function_2";
  discrete_function_2_dataoutput.writeData( 1.0 /*dummy*/, outstring_2.str() );
  // clear the std::stringstream:
  outstring_2.str(std::string());

  // -------------------------------------------------------
#endif

#if 1
  std :: cout << "Grid 1 size = " << grid_1.size(0) << std :: endl;
  std :: cout << "Grid 2 size = " << grid_2.size(0) << std :: endl;

  ImprovedL2Error< DiscreteFunctionType > improved_l2error;
  L2Error< DiscreteFunctionType > l2error_test;

  RangeType difference_L2_alternative = 0.0; //!improved_l2error.norm_L2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( discrete_function_2, discrete_function_1 );
  RangeType norm_L2_function_2 = l2error_test.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( zero_function_2, discrete_function_2 );
  RangeType difference_L2_test = 0.0; //!improved_l2error.norm_L2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( discrete_function_1, discrete_function_2 );

  RangeType difference_L2 = improved_l2error.norm_adaptive_grids_2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( discrete_function_1, discrete_function_2 );


  std :: cout << "L2 difference  = " << difference_L2 << std :: endl;
//!  std :: cout << "L2 difference (check) = " << difference_L2_test << std :: endl;
  std :: cout << "L2 norm function_2 = " << norm_L2_function_2 << std :: endl;
//!  std :: cout << "L2 difference (alternative check) = " << difference_L2_alternative << std :: endl;

#endif
#endif


  std::abort();

#endif
 //! --------------------------------------------------------------------------------------




  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  DiffusionType diffusion_op;

  // define (first) source term:
  FirstSourceType f; // standard source f

  // if we have some additional source term (-div G), define:
  SecondSourceType G;
  // - div ( A^{\epsilon} \nabla u^{\epsilon} ) = f - div G

  //! Ueberdenken, ob wir das nicht rausschmeisen und nur im Hintergrund fuer die Zellprobleme verwenden:
  // define mass (just for cell problems \lambda w - \div A \nabla w = rhs)
  MassTermType mass;


  // dummy coefficient (mass, advection, etc.)
  DefaultDummyFunctionType dummy_coeff;

  // exact solution unknown?
#ifdef EXACTSOLUTION_AVAILABLE
  ExactSolutionType u;
  DiscreteExactSolutionType discrete_exact_solution( "discrete exact solution ", u, gridPartFine );
#endif

  //! --------------------------------------------------------------------------------------




  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  RightHandSideAssembler< DiscreteFunctionType > rhsassembler;


  //----------------------------------------------------------------------------------------------//
  //----------------------- THE DISCRETE FEM OPERATOR -----------------------------------//
  //----------------------------------------------------------------------------------------------//

  //! define the discrete (elliptic) operator that describes our problem
  // ( effect of the discretized differential operator on a certain discrete function )
  EllipticOperatorType discrete_elliptic_op( finerDiscreteFunctionSpace, diffusion_op);

  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//

  RangeType size_of_domain = get_size_of_domain(discreteFunctionSpace);


#ifdef HOMOGENIZEDSOL_AVAILABLE

 #ifdef LINEAR_PROBLEM

  std :: string unit_cell_location = "../dune/multiscale/grids/cell_grids/unit_cube.dgf";
  FieldMatrix< RangeType, dimension, dimension > A_hom;

  //analytical homogenizer:

  #if 0
  AnalyticalHomogenizer< GridType, DiffusionType > ana_homogenizer( unit_cell_location );
  A_hom = ana_homogenizer.getHomTensor(diffusion_op);
  #endif


  // descretized homogenizer:

  #if 1
  Homogenizer< GridType, DiffusionType > disc_homogenizer( unit_cell_location );
  A_hom = disc_homogenizer.getHomTensor(diffusion_op);
  #endif

  // if known a-priori

  #if 0
  A_hom[0][0] = 1.73224;
  A_hom[1][1] = 2.0;
  A_hom[0][1] = 0.0;
  A_hom[1][0] = 0.0;
  #endif


  typedef FieldMatrix < RangeType, dimension, dimension > FieldMatrixType;
  Problem::HomDiffusion< FunctionSpaceType, FieldMatrixType > hom_diffusion_op( A_hom);

  typedef DiscreteEllipticOperator< DiscreteFunctionType, Problem::HomDiffusion<FunctionSpaceType, FieldMatrixType>, MassTermType > HomEllipticOperatorType;

  HomEllipticOperatorType hom_discrete_elliptic_op( finerDiscreteFunctionSpace, hom_diffusion_op);

  FEMMatrix hom_stiff_matrix( "homogenized stiffness matrix", finerDiscreteFunctionSpace, finerDiscreteFunctionSpace );

  DiscreteFunctionType hom_rhs( "homogenized rhs", finerDiscreteFunctionSpace );
  hom_rhs.clear();

  DiscreteFunctionType homogenized_solution( filename_ + " Homogenized Solution", finerDiscreteFunctionSpace );
  homogenized_solution.clear();
  hom_discrete_elliptic_op.assemble_matrix( hom_stiff_matrix );

  rhsassembler.assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( f , hom_rhs);

  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType hom_endit = finerDiscreteFunctionSpace.end();
  for( IteratorType fine_it = finerDiscreteFunctionSpace.begin(); fine_it != hom_endit; ++fine_it )
      boundaryTreatment( *fine_it , hom_rhs );

  InverseFEMMatrix hom_biCGStab( hom_stiff_matrix, 1e-8, 1e-8, 20000, VERBOSE );
  hom_biCGStab( hom_rhs, homogenized_solution );

 #else

 #endif

#endif


#ifdef FINE_SCALE_REFERENCE

 #ifdef FSR_COMPUTE
//! *******************************************************************

  // starting value for the Newton method
  DiscreteFunctionType zero_func( filename_ + " constant zero function ", finerDiscreteFunctionSpace );
  zero_func.clear();


  //! *************************** Assembling the reference problem ****************************
  // ( fine scale reference solution = fem_newton_solution )

  //! (stiffness) matrix
  FEMMatrix fem_newton_matrix( "FEM Newton stiffness matrix", finerDiscreteFunctionSpace, finerDiscreteFunctionSpace );

  //! right hand side vector
  // right hand side for the finite element method with Newton solver:
  // ( also right hand side for the finer discrete function space )
  DiscreteFunctionType fem_newton_rhs( "fem newton rhs", finerDiscreteFunctionSpace );
  fem_newton_rhs.clear();

  //! solution vector
  // solution of the finite element method, where we used the Newton method to solve the non-linear system of equations
  // in general this will be an accurate approximation of the exact solution, that is why we it also called reference solution
  DiscreteFunctionType fem_newton_solution( filename_ + " Reference (FEM Newton) Solution", finerDiscreteFunctionSpace );
  fem_newton_solution.clear();
  // By fem_newton_solution, we denote the "fine scale reference solution" (used for comparison)
  // ( if the elliptic problem is linear, the 'fem_newton_solution' is determined without the Newton method )

#ifdef LINEAR_PROBLEM

  std :: cout << "Solving linear problem." << std :: endl;
  if (data_file.is_open())
    {
      data_file << "Solving linear problem with standard FEM and resolution level " << problem_info.getRefinementLevelReferenceProblem() << "." << std :: endl;
      data_file << "------------------------------------------------------------------------------" << std :: endl;
    }

  // to assemble the computational time
  Dune::Timer assembleTimer;

  // assemble the stiffness matrix
  discrete_elliptic_op.assemble_matrix( fem_newton_matrix );

  std::cout << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
  if (data_file.is_open())
    {
      data_file << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
    }

  // assemble right hand side
  rhsassembler.assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( f , fem_newton_rhs);

  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType fine_endit = finerDiscreteFunctionSpace.end();
  for( IteratorType fine_it = finerDiscreteFunctionSpace.begin(); fine_it != fine_endit; ++fine_it )
         boundaryTreatment( *fine_it , fem_newton_rhs );


  InverseFEMMatrix fem_biCGStab( fem_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE );
  fem_biCGStab( fem_newton_rhs, fem_newton_solution );

  if (data_file.is_open())
    {
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
      data_file << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
    }

// if non-linear problem
#else

  std :: cout << "Solving non-linear problem." << std :: endl;
  if (data_file.is_open())
    {
      data_file << "Solving nonlinear problem with FEM + Newton-Method. Resolution level of grid = " << problem_info.getRefinementLevelReferenceProblem() << "." << std :: endl;
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
    }

  Dune::Timer assembleTimer;

  //! residual vector
  // current residual
  DiscreteFunctionType fem_newton_residual( filename_ + "FEM Newton Residual", finerDiscreteFunctionSpace );
  fem_newton_residual.clear();

  RangeType relative_newton_error_finescale = 10000.0;
  RangeType rhs_L2_norm = 10000.0;

  int iteration_step = 1;
  // the Newton step for the FEM reference problem (solved with Newton Method):
  // L2-Norm of residual < tolerance ?
  double tolerance = 1e-06;
  while( relative_newton_error_finescale > tolerance )
   {
     // (here: fem_newton_solution = solution from the last iteration step)

     std :: cout << "Newton iteration " << iteration_step << ":" << std :: endl;
     if (data_file.is_open())
      { data_file << "Newton iteration " << iteration_step << ":" << std :: endl; }

     Dune::Timer stepAssembleTimer;

     // assemble the stiffness matrix
     discrete_elliptic_op.assemble_jacobian_matrix( fem_newton_solution, fem_newton_matrix );

     std::cout << "Time to assemble FEM Newton stiffness matrix for current iteration: " << stepAssembleTimer.elapsed() << "s" << std::endl;

     // assemble right hand side
     rhsassembler.assemble_for_Newton_method< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( f , diffusion_op, fem_newton_solution, fem_newton_rhs);

     rhs_L2_norm = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( fem_newton_rhs, zero_func );

     if ( rhs_L2_norm < 1e-10 )
      {
        // residual solution almost identical to zero: break
        if (data_file.is_open())
           {
             data_file << "Residual solution almost identical to zero. Therefore: break loop." << std :: endl;
             data_file << "(L^2-Norm of current right hand side = " << rhs_L2_norm << " < 1e-10)" << std :: endl;
           }
        break;
      }

     // set Dirichlet Boundary to zero 
     typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
     IteratorType fine_endit = finerDiscreteFunctionSpace.end();
     for( IteratorType fine_it = finerDiscreteFunctionSpace.begin(); fine_it != fine_endit; ++fine_it )
         boundaryTreatment( *fine_it , fem_newton_rhs );


     InverseFEMMatrix fem_newton_biCGStab( fem_newton_matrix, 1e-8, 1e-8, 20000, true );
     fem_newton_biCGStab( fem_newton_rhs, fem_newton_residual );

     if( fem_newton_residual.dofsValid() )
      {

        fem_newton_solution += fem_newton_residual;

        relative_newton_error_finescale = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( fem_newton_residual, zero_func );
        relative_newton_error_finescale /= l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( fem_newton_solution, zero_func );

        std :: cout << "Relative L2-Newton Error = " << relative_newton_error_finescale << std :: endl;
        // residual solution almost identical to zero: break
        if (data_file.is_open())
           {
             data_file << "Relative L2-Newton Error = " << relative_newton_error_finescale << std :: endl;
             if ( relative_newton_error_finescale <= tolerance )
              {
                data_file << "Since tolerance = " << tolerance << ": break loop." << std :: endl;
              }
           }

        fem_newton_residual.clear();

      }
     else
      {
        std :: cout << "WARNING! Invalid dofs in 'fem_newton_residual'." << std :: endl;
        break;
      }

     iteration_step += 1;

   }

  std :: cout << "Problem with FEM + Newton-Method solved in " << assembleTimer.elapsed() << "s." << std :: endl << std :: endl;
  if (data_file.is_open())
    {
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
      data_file << "Problem with FEM + Newton-Method solved in " << assembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
    }

#endif // end '#ifdef LINEAR_PROBLEM <-> #else'

  //! ********************** End of assembling the reference problem ***************************
#endif
//end FSR_COMPUTE

#ifdef FSR_LOAD

  DiscreteFunctionType fem_newton_solution( filename_ + " Reference (FEM Newton) Solution", finerDiscreteFunctionSpace );
  fem_newton_solution.clear();

  char modeprob[50];
  sprintf( modeprob, "/Model_Problem_%d", problem_info.get_Number_of_Model_Problem() );
  std::string modeprob_s(modeprob);

  char reference_solution_directory[50];
  sprintf( reference_solution_directory, "/reference_solution_ref_%d", refinement_level_referenceprob_ );
  std::string reference_solution_directory_s(reference_solution_directory);

  char reference_solution_name[50];
  sprintf( reference_solution_name, "/finescale_solution_discFunc_refLevel_%d", refinement_level_referenceprob_ );
  std::string reference_solution_name_s(reference_solution_name);

  std :: string location_fine_scale_ref = "data/HMM/" + modeprob_s + reference_solution_directory_s + reference_solution_name_s;

  bool reader_is_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_ref( (location_fine_scale_ref).c_str() );
  discrete_function_reader_ref.open();

  discrete_function_reader_ref.read( 0, fem_newton_solution );
  std :: cout << "fine scale reference read." << std :: endl;

#endif
//end FSR_LOAD

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

  char reference_hmm_solution_directory[50];
  sprintf( reference_hmm_solution_directory, ".......", 10 );
  std::string reference_solution_directory_s(reference_solution_directory);

  char reference_solution_name[50];
  sprintf( reference_solution_name, "....", refinement_level_referenceprob_ );
  std::string reference_solution_name_s(reference_solution_name);

  std :: string location_hmm_ref = "data/HMM/" + modeprob_s + reference_solution_directory_s + reference_solution_name_s;
#endif

  std :: string location_hmm_ref = "data/HMM/Model_Problem_1/Macro_10_Micro_8/hmm_solution_discFunc_refLevel_10";

  bool hmm_ref_reader_is_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_hmm_ref( (location_hmm_ref).c_str() );
  discrete_function_reader_hmm_ref.open();

  discrete_function_reader_hmm_ref.read( 0, hmm_reference_solution );
  std :: cout << "HMM reference read." << std :: endl;

#endif


  //! ************************* Assembling and solving the HMM problem ****************************

#ifdef ADAPTIVE
// number of the loop cycle of the while-loop
int loop_cycle = 1;
double total_hmm_time = 0.0;
bool repeat = true;
while ( repeat == true )
{

  if (data_file.is_open())
    {
      data_file << "########################### LOOP CYCLE " << loop_cycle << " ###########################" << std :: endl << std :: endl << std :: endl;
    }
#endif

  std::cout << std :: endl << "Solving HMM-macro-problem for " << discreteFunctionSpace.size()
            << " unkowns and polynomial order "
            << DiscreteFunctionSpaceType :: polynomialOrder << "." 
            << std :: endl << std :: endl;

  //----------------------------------------------------------------------------------------------//
  //----------------------- THE DISCRETE HMM OPERATOR -----------------------------------//
  //----------------------------------------------------------------------------------------------//

  // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain cell problem
  CellProblemNumberingManagerType cp_num_manager(discreteFunctionSpace);


  //! define the elliptic hmm operator that describes our 'homogenized' macro problem
  // ( effect of the elliptic hmm operator on a certain discrete function )
  EllipticHMMOperatorType discrete_elliptic_hmm_op( discreteFunctionSpace, periodicDiscreteFunctionSpace, diffusion_op, cp_num_manager, filename_);

  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//

  //! matrix
  FEMMatrix hmm_newton_matrix( "HMM Newton stiffness matrix", discreteFunctionSpace, discreteFunctionSpace );

  //! right hand side vector
  // right hand side for the hm finite element method with Newton solver:
  DiscreteFunctionType hmm_newton_rhs( "hmm rhs", discreteFunctionSpace );
  hmm_newton_rhs.clear();

  //! solution vector
  // solution of the heterogeneous multiscale finite element method, where we used the Newton method to solve the non-linear system of equations
  DiscreteFunctionType hmm_solution( filename_ + " HMM (+Newton) Solution", discreteFunctionSpace );
  hmm_solution.clear();

#ifdef ADAPTIVE
  // one fot the discreteFunctionSpace
  RestrictProlongOperatorType rp( hmm_solution );
  AdaptationManagerType adaptationManager( grid, rp );
#endif

  // starting value for the Newton method
  DiscreteFunctionType zero_func_coarse( filename_ + " constant zero function coarse ", discreteFunctionSpace );
  zero_func_coarse.clear();

#ifdef LINEAR_PROBLEM


 // solve cell problems in a preprocess, if AD_HOC_COMPUTATION is not defined
 #ifndef AD_HOC_COMPUTATION

  //! -------------- solve and save the cell problems for the base function set --------------------------------------

  CellProblemSolver< PeriodicDiscreteFunctionType, DiffusionType > cell_problem_solver(periodicDiscreteFunctionSpace, diffusion_op, data_file /*optinal*/);

  int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);

  std :: cout << "Solving cell problems for " << number_of_grid_elements << " leaf entities." << std :: endl;

  // generate directory for cell problem data output
  if (mkdir(("data/HMM/" + filename_ + "/cell_problems/").c_str() DIRMODUS) == -1)
   {
    std::cout << "WARNING! Directory for the solutions of the cell problems already exists!";
   }
  else
   {
     mkdir(("data/HMM/" + filename_ + "/cell_problems/").c_str() DIRMODUS);
   }

  // -------------- solve cell problems for the macro basefunction set ------------------------------
  // save the solutions of the cell problems for the set of macroscopic base functions

  cell_problem_solver.saveTheSolutions_baseSet< DiscreteFunctionType >( discreteFunctionSpace, cp_num_manager, filename_ + "/cell_problems/");

  // ------------- end solving and saving cell problems for the macro basefunction set --------------

  //! --------------- end solving and saving cell problems -----------------------------------------
 #endif

  std :: cout << "Solving linear HMM problem." << std :: endl;
  if (data_file.is_open())
    {
      data_file << "Solving linear HMM problem." << std :: endl;
      data_file << "------------------------------------------------------------------------------" << std :: endl;
    }

  // to assemble the computational time
  Dune::Timer hmmAssembleTimer;

  // assemble the hmm stiffness matrix
  discrete_elliptic_hmm_op.assemble_matrix( hmm_newton_matrix );
  // to print the matrix, use:   hmm_newton_matrix.print();

  std::cout << "Time to assemble HMM macro stiffness matrix: " << hmmAssembleTimer.elapsed() << "s" << std::endl;
  if (data_file.is_open())
    {
      data_file << "Time to assemble HMM macro stiffness matrix: " << hmmAssembleTimer.elapsed() << "s" << std::endl;
    }

  // assemble right hand side
  rhsassembler.assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( f , hmm_newton_rhs);

  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
     boundaryTreatment( *it , hmm_newton_rhs );

  InverseFEMMatrix hmm_biCGStab( hmm_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE );
  hmm_biCGStab( hmm_newton_rhs, hmm_solution );

  std :: cout << "Linear HMM problem solved in " << hmmAssembleTimer.elapsed() << "s." << std :: endl << std :: endl;
  if (data_file.is_open())
    {
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
      data_file << "Linear HMM problem solved in " << hmmAssembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
    }

// if non-linear problem
#else

  // the nonlinear case

  // solve cell problems in a preprocess, if AD_HOC_COMPUTATION is not defined
  #ifndef AD_HOC_COMPUTATION


   //! -------------- solve and save the cell problems for the macroscopic base function set --------------------------------------

   CellProblemSolver< PeriodicDiscreteFunctionType, DiffusionType > cell_problem_solver(periodicDiscreteFunctionSpace, diffusion_op, data_file /*optinal*/);

   int number_of_grid_elements = periodicDiscreteFunctionSpace.grid().size(0);

   std :: cout << "Start solving cell problems for " << number_of_grid_elements << " leaf entities..." << std :: endl;

   // generate directory for cell problem data output
   if (mkdir(("data/HMM/" + filename_ + "/cell_problems/").c_str() DIRMODUS) == -1)
    {
     std::cout << "WARNING! Directory for the solutions of the cell problems already  exists!";
    }
   else
    {
      mkdir(("data/HMM/" + filename_ + "/cell_problems/").c_str() DIRMODUS);
    }

   // only for the case with test function reconstruction:
   #ifdef TFR
   // -------------- solve cell problems for the macro basefunction set ------------------------------
   // save the solutions of the cell problems for the set of macroscopic base functions

   cell_problem_solver.saveTheSolutions_baseSet< DiscreteFunctionType >(  discreteFunctionSpace, cp_num_manager, filename_ + "/cell_problems/");

   std :: cout << "Solving the cell problems for the base function set succeeded." << std :: endl;
   //  end solving and saving cell problems
   #endif

   //! --------------- end solving and saving cell problems -----------------------------------------
  #endif

  std :: cout << "Solving nonlinear HMM problem." << std :: endl;
  if (data_file.is_open())
    {
      data_file << "Solving nonlinear HMM problem with Newton method." << std :: endl;
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
    }

  Dune::Timer hmmAssembleTimer;

  // just to provide some information
  PeriodicDiscreteFunctionType dummy_periodic_func( "a periodic dummy", periodicDiscreteFunctionSpace );
  dummy_periodic_func.clear();

  //! residual vector
  // current residual
  DiscreteFunctionType hmm_newton_residual( filename_ + "HMM Newton Residual", discreteFunctionSpace );
  hmm_newton_residual.clear();

  RangeType relative_newton_error = 10000.0;
  RangeType hmm_rhs_L2_norm = 10000.0;

  // number of HMM Newton step (1 = first step)
  // HMM_NEWTON_ITERATION_STEP' Netwon steps have been already performed,
  // the next one is 'HMM_NEWTON_ITERATION_STEP+1' = hmm_iteration_step
  int hmm_iteration_step = HMM_NEWTON_ITERATION_STEP + 1;


  #ifdef RESUME_TO_BROKEN_COMPUTATION
  
  //std :: string location_hmm_newton_step_solution = "data/HMM/test/hmm_solution_discFunc_refLevel_5_NewtonStep_2";

  char fnewtonname[50];
  sprintf( fnewtonname, "/hmm_solution_discFunc_refLevel_%d_NewtonStep_%d", refinement_level_macrogrid_, HMM_NEWTON_ITERATION_STEP );
  std :: string fnewtonname_s( fnewtonname );
  std :: string location_hmm_newton_step_solution = "data/HMM/" + filename_ + fnewtonname_s;

  bool reader_open = false;

  // reader for the cell problem data file:
  DiscreteFunctionReader discrete_function_reader_hmm_newton_ref( (location_hmm_newton_step_solution).c_str() );
  discrete_function_reader_hmm_newton_ref.open();

  discrete_function_reader_hmm_newton_ref.read( 0, hmm_solution );

  #endif


  double old_error = 100.0;
  double error_decay = 0.0;


  // the Newton step for the nonlinear HMM problem:
  // L2-Norm of residual < tolerance ?
  #ifdef STOCHASTIC_PERTURBATION
  double hmm_tolerance = 1e-01 * VARIANCE;
  #else
  double hmm_tolerance = 1e-05;
  #endif
  while( relative_newton_error > hmm_tolerance )
   {
     // (here: hmm_solution = solution from the last iteration step)

     long double newton_step_time = clock();

     std :: cout << "HMM Newton iteration " << hmm_iteration_step << ":" << std :: endl;
     if (data_file.is_open())
      { data_file << "HMM Newton iteration " << hmm_iteration_step << ":" << std :: endl; }

     #ifndef AD_HOC_COMPUTATION
     // solve cell problems for the solution of the last iteration step
     cell_problem_solver.saveTheSolutions_discFunc< DiscreteFunctionType >( hmm_solution, filename_ + "/cell_problems/");
     cell_problem_solver.saveTheJacCorSolutions_baseSet_discFunc< DiscreteFunctionType >( hmm_solution, cp_num_manager, filename_ + "/cell_problems/");
     #endif

     // to assemble the computational time
     Dune::Timer stepHmmAssembleTimer;

     // assemble the stiffness matrix
     discrete_elliptic_hmm_op.assemble_jacobian_matrix( hmm_solution, hmm_newton_matrix );

     std::cout << "Time to assemble HMM stiffness matrix for current Newton iteration: " << stepHmmAssembleTimer.elapsed() << "s" << std::endl;
     if (data_file.is_open())
      {
        data_file << "Time to assemble HMM stiffness matrix for current Newton iteration: " << stepHmmAssembleTimer.elapsed() << "s" << std::endl;
      }

     std::cout << "Assemble right hand side..." << std::endl;
     // assemble right hand side
     rhsassembler.assemble_for_HMM_Newton_method< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( f , diffusion_op, hmm_solution, cp_num_manager, dummy_periodic_func, hmm_newton_rhs, filename_);
     std::cout << "Right hand side assembled!" << std::endl;

     if ( !(hmm_newton_rhs.dofsValid()) )
       {
         std::cout << "Right hand side invalid!" << std::endl;
         data_file << "Right hand side invalid!" << std::endl;
         std :: abort();
       }
     else
       {
         std::cout << "Right hand side valid ";
       }

     hmm_rhs_L2_norm = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( zero_func_coarse, hmm_newton_rhs );

     std::cout << "with L^2-Norm = " << hmm_rhs_L2_norm << "." << std::endl;
     data_file << "Assembled right hand side, with L^2-Norm of RHS = " << hmm_rhs_L2_norm << "." << std::endl;


     if (  hmm_rhs_L2_norm < 1e-10 )
      {
        // residual solution almost identical to zero: break
        if (data_file.is_open())
           {
             data_file << "HMM residual solution almost identical to zero. Therefore: break loop." << std :: endl;
             data_file << "(L^2-Norm of current right hand side = " << hmm_rhs_L2_norm << " < 1e-10)" << std :: endl;
           }
        break;
      }

     // set Dirichlet Boundary to zero 
     typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
     IteratorType endit = discreteFunctionSpace.end();
     for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
         boundaryTreatment( *it , hmm_newton_rhs );

     #ifndef AD_HOC_COMPUTATION
     double hmm_biCG_tolerance = 1e-8;
     bool hmm_solution_convenient = false;
     while( hmm_solution_convenient == false )
       {

         hmm_newton_residual.clear();
         InverseFEMMatrix hmm_newton_biCGStab( hmm_newton_matrix,
                                               1e-8, hmm_biCG_tolerance, 20000, VERBOSE );

         hmm_newton_biCGStab( hmm_newton_rhs, hmm_newton_residual );

         if ( hmm_newton_residual.dofsValid() )
             { hmm_solution_convenient = true; }

         if ( hmm_biCG_tolerance > 1e-4 )
          {
            std :: cout << "WARNING! Iteration step " << hmm_iteration_step << ". Invalid dofs in 'hmm_newton_residual'." << std :: endl;
            std :: abort();
          }
         hmm_biCG_tolerance *= 10.0;

       }
     #else
     InverseFEMMatrix hmm_newton_biCGStab( hmm_newton_matrix, 1e-8, 1e-8, 20000, VERBOSE );
     hmm_newton_biCGStab( hmm_newton_rhs, hmm_newton_residual );
     #endif

     if( hmm_newton_residual.dofsValid() )
      {

        hmm_solution += hmm_newton_residual;
		

        // write the solution after the current HMM Newton step to a file
	#ifdef WRITE_HMM_SOL_TO_FILE

          // for adaptive computations, the saved solution is not suitable for a later usage
          #ifndef ADAPTIVE

	  bool writer_open = false;

          char fname[50];
          sprintf( fname, "/hmm_solution_discFunc_refLevel_%d_NewtonStep_%d", refinement_level_macrogrid_, hmm_iteration_step );
          std :: string fname_s( fname );

          std :: string location = "data/HMM/" + filename_ + fname_s;
          DiscreteFunctionWriter dfw( (location).c_str() );
          writer_open = dfw.open();
          if ( writer_open )
          dfw.append( hmm_solution );

          // if you want an utput for all newton steps, even for an adaptive computation, use:
 	  // #endif

          // writing paraview data output

          // general output parameters
          myDataOutputParameters outputparam;
          outputparam.set_path( "data/HMM/" + filename_ );

	  // sequence stamp
          std::stringstream outstring;

          // create and initialize output class
	  IOTupleType hmm_solution_newton_step_series( &hmm_solution );
          #ifdef ADAPTIVE
          char hmm_prefix[50];
          sprintf( hmm_prefix, "hmm_solution_%d_NewtonStep_%d", loop_cycle, hmm_iteration_step );
          #else
          char hmm_prefix[50];
          sprintf( hmm_prefix, "hmm_solution_NewtonStep_%d", hmm_iteration_step );		  
          #endif
	  outputparam.set_prefix( hmm_prefix );
	  DataOutputType hmmsol_dataoutput( gridPart.grid(), hmm_solution_newton_step_series, outputparam );

          // write data
          outstring << "hmm-solution-NewtonStep";
          hmmsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
          // clear the std::stringstream:
          outstring.str(std::string());

	  #endif

        #endif

        // || u^(n+1) - u^(n) ||_L2
        relative_newton_error = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( hmm_newton_residual, zero_func_coarse );
        // || u^(n+1) - u^(n) ||_L2 / || u^(n+1) ||_L2
        relative_newton_error = relative_newton_error / l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( hmm_solution, zero_func_coarse );

        std :: cout << "Relative L2 HMM Newton iteration error = " << relative_newton_error << std :: endl;

        // residual solution almost identical to zero: break
        if (data_file.is_open())
           {
             data_file << "Relative L2 HMM Newton iteration error = " << relative_newton_error << std :: endl;
             if ( relative_newton_error <= hmm_tolerance )
              {
                newton_step_time = clock() - newton_step_time;
                newton_step_time = newton_step_time / CLOCKS_PER_SEC;
                if (data_file.is_open())
                 {
                   data_file << std :: endl << "Total time for current HMM Newton step = " << newton_step_time << "s." << std :: endl << std :: endl;
                 }
                data_file << "Since HMM-tolerance = " << hmm_tolerance << ": break loop." << std :: endl;
                data_file << "....................................................." << std :: endl << std :: endl;
              }
           }

        hmm_newton_residual.clear();

      }
     else
      {
        std :: cout << "WARNING! Invalid dofs in 'hmm_newton_residual'." << std :: endl;
        break;
      }

     hmm_iteration_step += 1;

     if ( relative_newton_error > hmm_tolerance )
      {
        newton_step_time = clock() - newton_step_time;
        newton_step_time = newton_step_time / CLOCKS_PER_SEC;
        if (data_file.is_open())
         {
           data_file << std :: endl << "Total time for current HMM Newton step = " << newton_step_time << "s." << std :: endl << std :: endl;

           error_decay = relative_newton_error / old_error;
           old_error = relative_newton_error;
           // maximum number of Newton iterations
           if ( (hmm_iteration_step >= 20) || (error_decay >= 0.95) )
            {
              data_file << std :: endl << "Reached constant or inceasing error decay or maximum number of Newton iterations:  break loop." << std :: endl;
              data_file << "....................................................." << std :: endl << std :: endl;
              break;
            }
           data_file << "....................................................." << std :: endl << std :: endl;
         }
      }

   }

  std :: cout << "HMM problem with Newton method solved in " << hmmAssembleTimer.elapsed() << "s." << std :: endl << std :: endl;
  if (data_file.is_open())
    {
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
      data_file << "HMM problem with Newton method solved in " << hmmAssembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
    }
  #ifdef ADAPTIVE
  total_hmm_time += hmmAssembleTimer.elapsed();
  #endif
#endif
//end if not defined LINEAR_PROBLEM


//! ---------- Error Estimation ----------
#ifdef ERRORESTIMATION

  // Notation: u_H = hmm_solution

  // to load the solutions of the cell problems:

  // location of the solutions of the cell problems for the base function set:
  std :: string cell_solution_location_baseSet;
  // location of the solutions of the cell problems for the discrete function u_H:
  std :: string cell_solution_location_discFunc;

  cell_solution_location_baseSet = "data/HMM/"+filename_+"/cell_problems/_cellSolutions_baseSet";
  cell_solution_location_discFunc = "data/HMM/"+filename_+"/cell_problems/_cellSolutions_discFunc";

  // reader for the cell problem data file (for tha macro base set):
  DiscreteFunctionReader discrete_function_reader_baseSet( (cell_solution_location_baseSet).c_str() );
  discrete_function_reader_baseSet.open();

  // reader for the cell problem data file (for u_H):
  DiscreteFunctionReader discrete_function_reader_discFunc( (cell_solution_location_discFunc).c_str() );
  discrete_function_reader_discFunc.open();


  ErrorEstimatorType error_estimator( periodicDiscreteFunctionSpace,
                                      discreteFunctionSpace,
                                      auxiliaryDiscreteFunctionSpace,
                                      diffusion_op );

  RangeType estimated_source_error = 0.0;
  RangeType estimated_approximation_error = 0.0;
  RangeType estimated_residual_error = 0.0;

  // estimated_residual_error = sqrt( estimated_residual_error_micro_jumps^2 + estimated_residual_error_macro_jumps^2 )
  RangeType estimated_residual_error_micro_jumps = 0.0;
  RangeType estimated_residual_error_macro_jumps = 0.0;
#ifdef TFR
  RangeType estimated_tfr_error = 0.0;
#endif

  RangeType local_error_indicator[discreteFunctionSpace.grid().size(0)];

  RangeType minimal_loc_indicator = 10000.0;
  RangeType maximal_loc_indicator = 0.0;
  RangeType average_loc_indicator = 0.0;

  int element_number = 0;
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
     {

       const EntityType& entity = *it;

       // corrector of u_H^(n-1) \approx u_H on the macro element T
       PeriodicDiscreteFunctionType corrector_u_H_on_entity( "Corrector of u_H", periodicDiscreteFunctionSpace );
       corrector_u_H_on_entity.clear();


       // in the linear case, we still need to compute the corrector of u_H:
       #ifdef LINEAR_PROBLEM
       PeriodicDiscreteFunctionType corrector_of_base_func( "Corrector of macro base function", periodicDiscreteFunctionSpace );
       corrector_of_base_func.clear();

       DiscreteFunctionType::LocalFunctionType local_hmm_solution = hmm_solution.localFunction( entity );

       const BaseFunctionSetType &baseSet = discreteFunctionSpace.baseFunctionSet(entity);
       const unsigned int numMacroBaseFunctions = baseSet.numBaseFunctions();
       int cell_problem_id [ numMacroBaseFunctions ];
       for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {
          cell_problem_id[ i ] = cp_num_manager.get_number_of_cell_problem( it, i );

          discrete_function_reader_baseSet.read( cell_problem_id[ i ], corrector_of_base_func );

          corrector_of_base_func *= local_hmm_solution[ i ];

          corrector_u_H_on_entity += corrector_of_base_func;

          corrector_of_base_func.clear();

        }
       #else
       // in the nonlinear case this corrector is already available
       discrete_function_reader_discFunc.read( element_number, corrector_u_H_on_entity );
       #endif


       // contribution of the local source error

       RangeType local_source_indicator = error_estimator.indicator_f( f , entity );

       estimated_source_error += local_source_indicator;


       // contribution of the local approximation error

       RangeType local_approximation_indicator = error_estimator.indicator_app_1( entity, hmm_solution, corrector_u_H_on_entity );

       local_approximation_indicator += error_estimator.indicator_app_2( entity, hmm_solution, corrector_u_H_on_entity );

       estimated_approximation_error += local_approximation_indicator;


       // contribution of the local residual error

       RangeType local_residual_indicator = error_estimator.indicator_res_T( entity, hmm_solution, corrector_u_H_on_entity );

       estimated_residual_error_micro_jumps += local_residual_indicator;

       typedef GridPartType :: IntersectionIteratorType IntersectionIteratorType;
       IntersectionIteratorType endnit = gridPart.iend(entity);
       for(IntersectionIteratorType nit = gridPart.ibegin(entity); nit != endnit ; ++nit)
         {

           if ( nit->neighbor() ) //if there is a neighbor entity
            {

              // corrector of u_H^(n-1) \approx u_H on the neighbor element
              PeriodicDiscreteFunctionType corrector_u_H_on_neighbor_entity( "Corrector of u_H", periodicDiscreteFunctionSpace );
              corrector_u_H_on_neighbor_entity.clear();

              EntityPointerType it_outside = nit->outside();
              const EntityType& entity_outside = *it_outside;

              // in the linear case, we still need to compute the corrector of u_H:
              #ifdef LINEAR_PROBLEM
              PeriodicDiscreteFunctionType corrector_of_base_func_neighbor( "Corrector of macro base function", periodicDiscreteFunctionSpace );
              corrector_of_base_func_neighbor.clear();

              DiscreteFunctionType::LocalFunctionType local_hmm_solution_neighbor = hmm_solution.localFunction( entity_outside );

              const BaseFunctionSetType &baseSet_neighbor = discreteFunctionSpace.baseFunctionSet(entity_outside);
              const unsigned int numMacroBaseFunctions_neighbor = baseSet_neighbor.numBaseFunctions();
              int cell_problem_id_neighbor [ numMacroBaseFunctions_neighbor ];
              for( unsigned int i = 0; i < numMacroBaseFunctions_neighbor; ++i )
                {
                  cell_problem_id_neighbor[ i ] = cp_num_manager.get_number_of_cell_problem( it_outside, i );

                  discrete_function_reader_baseSet.read( cell_problem_id_neighbor[ i ], corrector_of_base_func_neighbor );

                  corrector_of_base_func_neighbor *= local_hmm_solution_neighbor[ i ];

                  corrector_u_H_on_neighbor_entity += corrector_of_base_func_neighbor;

                  corrector_of_base_func_neighbor.clear();

                }
              #else
              int neighbor_element_number = cp_num_manager.get_number_of_cell_problem( it_outside );
              // in the nonlinear case this corrector is already available
              discrete_function_reader_discFunc.read( neighbor_element_number, corrector_u_H_on_neighbor_entity );
              #endif

              RangeType val = error_estimator.indicator_res_E( *nit, hmm_solution, corrector_u_H_on_entity, corrector_u_H_on_neighbor_entity );

              local_residual_indicator += val;

              estimated_residual_error_macro_jumps += val;

            }

         }

       estimated_residual_error += local_residual_indicator;

#ifdef TFR
       // use 'indicator_effective_tfr' or 'indicator_tfr_1'
       // contribution of the local tfr error:

       RangeType local_tfr_indicator = error_estimator.indicator_tfr_1( entity, hmm_solution, corrector_u_H_on_entity );

       estimated_tfr_error += local_tfr_indicator;
#endif

       local_error_indicator[element_number] = local_source_indicator +
                                               local_approximation_indicator +
                                               local_residual_indicator;

#ifdef TFR
       local_error_indicator[element_number] += local_tfr_indicator;
#endif

       if ( local_error_indicator[element_number] < minimal_loc_indicator )
        { minimal_loc_indicator = local_error_indicator[element_number]; }

       if ( local_error_indicator[element_number] > maximal_loc_indicator )
        { maximal_loc_indicator = local_error_indicator[element_number]; }

       average_loc_indicator += local_error_indicator[element_number];

       element_number += 1;

     }

  average_loc_indicator /= discreteFunctionSpace.grid().size(0);

  estimated_source_error = sqrt(estimated_source_error);
  estimated_approximation_error = sqrt(estimated_approximation_error);
  estimated_residual_error = sqrt(estimated_residual_error);

  estimated_residual_error_micro_jumps = sqrt( estimated_residual_error_micro_jumps );
  estimated_residual_error_macro_jumps = sqrt( estimated_residual_error_macro_jumps );

  RangeType estimated_error = estimated_source_error + estimated_approximation_error + estimated_residual_error;

#ifdef TFR
  estimated_tfr_error = sqrt(estimated_tfr_error);
  estimated_error += estimated_tfr_error;
#endif

  #ifdef ADAPTIVE

  // maximum variation (up) from average 
  double max_variation = average_loc_indicator / maximal_loc_indicator;
  double min_variation = average_loc_indicator / minimal_loc_indicator;

  #endif

//! -------- End Error Estimation --------
#endif





#ifdef WRITE_HMM_SOL_TO_FILE

 // for adaptive computations, the saved solution is not suitable for a later usage
 #ifndef ADAPTIVE

  bool writer_is_open = false;

  char fname[40];
  sprintf( fname, "/hmm_solution_discFunc_refLevel_%d", refinement_level_macrogrid_ );
  std :: string fname_s( fname );

  std :: string location = "data/HMM/" + filename_ + fname_s;
  DiscreteFunctionWriter dfw( (location).c_str() );
  writer_is_open = dfw.open();
  if ( writer_is_open )
    dfw.append( hmm_solution );

 #endif

#endif


#ifdef FINE_SCALE_REFERENCE

 #ifdef WRITE_FINESCALE_SOL_TO_FILE

  bool fine_writer_is_open = false;

  char fine_fname[50];
  sprintf( fine_fname, "/finescale_solution_discFunc_refLevel_%d", refinement_level_referenceprob_ );
  std :: string fine_fname_s( fine_fname );

  std :: string fine_location = "data/HMM/" + filename_ + fine_fname_s;
  DiscreteFunctionWriter fine_dfw( (fine_location).c_str() );
  fine_writer_is_open = fine_dfw.open();
  if ( fine_writer_is_open )
    fine_dfw.append( fem_newton_solution );

 #endif

#endif

  //! ******************** End of assembling and solving the HMM problem ***************************

  std :: cout << std :: endl << "The L2 errors:" << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "The L2 errors:" << std :: endl << std :: endl; }

  //! ----------------- compute L2-errors -------------------

#ifdef FINE_SCALE_REFERENCE
  long double timeadapt = clock();

  RangeType hmm_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(hmm_solution,fem_newton_solution);

  std :: cout << "|| u_hmm - u_fine_scale ||_L2 =  " << hmm_error << std :: endl << std :: endl;
  if (data_file.is_open())
         { data_file << "|| u_hmm - u_fine_scale ||_L2 =  " << hmm_error << std :: endl; }

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


#ifdef HMM_REFERENCE
  long double timeadapthmmref = clock();

  RangeType hmm_vs_hmm_ref_error = impL2error.norm_adaptive_grids_2< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >(hmm_solution,hmm_reference_solution);

  std :: cout << "|| u_hmm - u_hmm_ref ||_L2 =  " << hmm_vs_hmm_ref_error << std :: endl << std :: endl;
  if (data_file.is_open())
         { data_file << "|| u_hmm - u_hmm_ref ||_L2 =  " << hmm_vs_hmm_ref_error << std :: endl; }

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

  #ifdef FINE_SCALE_REFERENCE
  // not yet modified according to a generalized L2-error, here, homogenized_solution and fem_newton_solution still need to be defined on the same grid!
  RangeType hom_error = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( homogenized_solution, fem_newton_solution );

  std :: cout << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_hom - u_fine_scale ||_L2 =  " << hom_error << std :: endl; }
  #endif

  RangeType hom_hmm_error = l2error.norm2<2 * DiscreteFunctionSpaceType :: polynomialOrder + 2>( hmm_solution, homogenized_solution );

  std :: cout << "|| u_hom - u_hmm ||_L2 =  " << hom_hmm_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_hom - u_hmm ||_L2 =  " << hom_hmm_error << std :: endl; }
#endif


#ifdef EXACTSOLUTION_AVAILABLE
  RangeType exact_hmm_error = l2error.norm< ExactSolutionType >( u, hmm_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_hmm - u_exact ||_L2 =  " << exact_hmm_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_hmm - u_exact ||_L2 =  " << exact_hmm_error << std :: endl; }

 #ifdef FINE_SCALE_REFERENCE
  RangeType fem_newton_error = l2error.norm< ExactSolutionType >( u, fem_newton_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_fem_newton - u_exact ||_L2 =  " << fem_newton_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_fem_newton - u_exact ||_L2 =  " << fem_newton_error << std :: endl; }
 #endif

#endif


#ifdef ERRORESTIMATION
  std :: cout << "Estimated error = " << estimated_error << "." << std :: endl;
  std :: cout << "In detail:" << std :: endl;
  std :: cout << "   Estimated source error = " << estimated_source_error << "." << std :: endl;
  std :: cout << "   Estimated approximation error = " << estimated_approximation_error << "." << std :: endl;
  std :: cout << "   Estimated residual error = " << estimated_residual_error << ", where:" << std :: endl;
  std :: cout << "        contribution of macro jumps = " << estimated_residual_error_macro_jumps << " and " << std :: endl;
  std :: cout << "        contribution of micro jumps = " << estimated_residual_error_micro_jumps << " and " << std :: endl;
  #ifdef TFR
   std :: cout << "   Estimated tfr error = " << estimated_tfr_error << "." << std :: endl;
  #endif
  if (data_file.is_open())
    {
      data_file << "Estimated error = " << estimated_error << "." << std :: endl;
      data_file << "In detail:" << std :: endl;
      data_file << "   Estimated source error = " << estimated_source_error << "." << std :: endl;
      data_file << "   Estimated approximation error = " << estimated_approximation_error << "." << std :: endl;
      data_file << "   Estimated residual error = " << estimated_residual_error << ", where:" << std :: endl;
      data_file << "        contribution of macro jumps = " << estimated_residual_error_macro_jumps << " and " << std :: endl;
      data_file << "        contribution of micro jumps = " << estimated_residual_error_micro_jumps << " and " << std :: endl;
      #ifdef TFR
       data_file << "   Estimated tfr error = " << estimated_tfr_error << "." << std :: endl;
      #endif
    }
#endif

//! -------------------------------------------------------


//! --------------- writing data output ---------------------


  // general output parameters
  myDataOutputParameters outputparam;
  outputparam.set_path( "data/HMM/" + filename_ );

  // sequence stamp
  std::stringstream outstring;


  // --------- data output hmm solution --------------

  // create and initialize output class
  IOTupleType hmm_solution_series( &hmm_solution );
  #ifdef ADAPTIVE
  char hmm_prefix[30];
  sprintf( hmm_prefix, "hmm_solution_%d", loop_cycle );
  outputparam.set_prefix( hmm_prefix );
  #else
  outputparam.set_prefix("hmm_solution");
  #endif
  DataOutputType hmmsol_dataoutput( gridPart.grid(), hmm_solution_series, outputparam );

  // write data
  outstring << "hmm-solution";
  hmmsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------


#ifdef EXACTSOLUTION_AVAILABLE
  // --------- data output discrete exact solution --------------

  // create and initialize output class
  ExSolIOTupleType exact_solution_series( &discrete_exact_solution );
  outputparam.set_prefix("exact_solution");
  ExSolDataOutputType exactsol_dataoutput( gridPartFine.grid(), exact_solution_series, outputparam );

  // write data
  outstring << "exact-solution";
  exactsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());
  // -------------------------------------------------------
#endif


#ifdef WRITE_FINESCALE_SOL_TO_FILE
  // --------- data output reference solution (fine fem newton computation) --------------

  // create and initialize output class
  IOTupleType fem_newton_solution_series( &fem_newton_solution );
  outputparam.set_prefix("reference_solution");
  DataOutputType fem_newton_dataoutput( gridPartFine.grid(), fem_newton_solution_series, outputparam );

  // write data
  outstring << "reference_solution";
  fem_newton_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------
#endif


#ifdef HOMOGENIZEDSOL_AVAILABLE
  // --------- data output homogenized solution --------------

  // create and initialize output class
  IOTupleType homogenized_solution_series( &homogenized_solution );
  outputparam.set_prefix("homogenized_solution");
  DataOutputType homogenized_solution_dataoutput( gridPartFine.grid(), homogenized_solution_series, outputparam );

  // write data
  outstring << "homogenized_solution";
  homogenized_solution_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------
#endif


//!-------------------------------------------------------------




#ifdef ADAPTIVE


  int default_refinement = 0;

  //double error_tolerance_ = 0.2;

  // int(...) rundet ab zum naechsten Integer
  // assuming we had a quadratic order of convergence of the error estimator and 
  // that we have a certain estimated error for the current uniform grid, than we can compute how many additional uniform refinements are required to get under the error (estimator) tolerance:
  // number_of_uniform_refinments = int( sqrt( (\eta_have) / (\eta_want) ) )
  // (obtained from the EOC formula)
  // uniform contribution only for the first loop cycle
  if ( loop_cycle == 1 )
   {

     // "divided by 2.0" we go half the way with a uniform computation 
     int number_of_uniform_refinements = 2*int(int( sqrt( estimated_error / error_tolerance_ ) ) / 2.0);

      if ( data_file.is_open() )
       {
         data_file << std :: endl << "Uniform default refinement:" << std :: endl << std :: endl;
         data_file << "sqrt( estimated_error / error_tolerance_ ) = " << sqrt( estimated_error / error_tolerance_ ) << std :: endl;
         data_file << "number_of_uniform_refinements = " << number_of_uniform_refinements << std :: endl;
         data_file << "***************" << std :: endl << std :: endl;
       }

     default_refinement = number_of_uniform_refinements;
   }


  int number_of_areas;
  if ( loop_cycle == 1 )
   {
     number_of_areas = 1;
   }
  else
   {
     if ( loop_cycle == 3 )
      { number_of_areas = 2; }
     else
      { number_of_areas = 2; }
   }

#if 0
  if ( loop_cycle == 2 )
     number_of_areas = 2;
#endif

  double border[number_of_areas-1];
  border[0] = 0.5;
  for( int bo = 1; bo < (number_of_areas-1) ; ++bo )
    { border[bo] = border[bo-1] + ((1.0 - border[bo-1])/2.0); }

  // 3 areas: 1: |0-30%| 2: |30-80%| 3: |80-100%|
  //border[0] = 0.3;
  //border[1] = 0.8;
  //border[2] = 0.95;

  int refinements_in_area[number_of_areas];
  for( int bo = 0; bo < number_of_areas ; ++bo )
    { refinements_in_area[bo] = default_refinement + bo + 1; }

#if 0
if ( loop_cycle > 1 ) { refinements_in_area[bo] = 6; }
  if ( loop_cycle == 2 )
   {
     //go the full way 
     int number_of_uniform_refinements = 2*int(int( sqrt( estimated_error / error_tolerance_ ) ) );
     refinements_in_area[0] = number_of_uniform_refinements + 2.0;
     refinements_in_area[1] = number_of_uniform_refinements + 2.0;
   }
#endif

  if ( data_file.is_open() )
    {
      data_file << "Adaption strategy:" << std :: endl << std :: endl;
      data_file << "Define 'variation = (indicator_on_element - average_indicator) / (maximum_indicator - average_indicator)'" << std :: endl;
      data_file << "Subdivide the region [average_indicator,maximum_indicator] into " << number_of_areas << " areas." << std :: endl;
      if ( number_of_areas == 1 )
       {
         data_file << "1.: [average_indicator,maximum_indicator]. Mark elements for " << refinements_in_area[0] << " refinements." << std :: endl;
       }
      else
       {
         data_file << "1.: [average_indicator," << border[0] << "*maximum_indicator]. If 'variance' in area: mark elements for " << refinements_in_area[0] << " refinements." << std :: endl;
         for( int bo = 1; bo < (number_of_areas-1) ; ++bo )
           data_file << bo+1 << ".: [" << border[bo-1] << "*average_indicator," << border[bo] << "*maximum_indicator]. If 'variance' in area: mark elements for " << refinements_in_area[bo] << " refinements." << std :: endl;
         data_file << number_of_areas << ".: [" << border[number_of_areas-2] << "*average_indicator,maximum_indicator]. If 'variance' in area: mark elements for " << refinements_in_area[number_of_areas-1] << " refinements." << std :: endl;
       }
      data_file << "Default refinement for elements with 'variance <= 0 ': " << default_refinement << std :: endl;
    }


  if ( estimated_error < error_tolerance_ )
   {
     repeat = false;
     std :: cout << "Total HMM time = " << total_hmm_time << "s." << std :: endl;
     data_file << std :: endl << std :: endl << "Total HMM time = " << total_hmm_time << "s." << std :: endl << std :: endl;

   }
  else
   {

    element_number = 0;
    typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
    IteratorType endit_test = discreteFunctionSpace.end();
    for( IteratorType it = discreteFunctionSpace.begin(); it != endit_test; ++it )
      {

        int additional_refinement;

        if ( local_error_indicator[element_number] <= average_loc_indicator )
          { additional_refinement = default_refinement; }
        else
          {

            double variation = (local_error_indicator[element_number] - average_loc_indicator ) / 
                        ( maximal_loc_indicator - average_loc_indicator );

            if ( number_of_areas == 1 )
             {
               additional_refinement = refinements_in_area[0];
             }
            else
             {
               if ( variation <= border[0] )
                 { additional_refinement = refinements_in_area[0]; }

               for( int bo = 1; bo <= (number_of_areas-2) ; ++bo )
                 if ( ( variation > border[bo-1] ) && ( variation <= border[bo] ) )
                  { additional_refinement = refinements_in_area[bo]; }

               if ( variation > border[number_of_areas-2] )
                 { additional_refinement = refinements_in_area[number_of_areas-1]; }
             }

          }

        grid.mark( additional_refinement , *it );
        element_number += 1;

      }

    adaptationManager.adapt();

  }

  if (data_file.is_open())
    {
      data_file << std :: endl << "#########################################################################" << std :: endl << std :: endl << std :: endl;
    }

  loop_cycle += 1;

}// end of repeat loop (for the adaptive cicles)
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

  // generate directories for data output
  if (mkdir(("data/HMM/" + filename_).c_str() DIRMODUS) == -1)
   {
    std::cout << "Directory already exists! Overwrite? y/n: ";
    char answer;
    std :: cin >> answer;
    if (!(answer=='y'))
     {std :: abort();}
   }
  else
   {
     mkdir(("data/HMM/" + filename_).c_str() DIRMODUS);
   }



  #ifdef RESUME_TO_BROKEN_COMPUTATION
  // man koennte hier noch den genauen Iterationsschritt in den Namen mit einfliessen lassen:
  // (vorlauefig sollte diese Variante aber reichen) 
  std :: string save_filename = "data/HMM/" + filename_ + "/problem-info-resumed-computation.txt";  
  #else
  std :: string save_filename = "data/HMM/" + filename_ + "/problem-info.txt";
  #endif
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
  int refinement_level_cellgrid;
  std :: cout << "Enter refinement level for the cell problems: ";
  std :: cin >> refinement_level_cellgrid;
  std :: cout << std :: endl;

#ifdef ADAPTIVE
  // error tolerance (comparison with upper bound determined by the a-posteriori error estimate)
  std :: cout << "Enter global error tolerance for program abort: ";
  std :: cin >> error_tolerance_;
  std :: cout << std :: endl;
#endif


  // (starting) grid refinement level for solving the reference problem
  refinement_level_referenceprob_ = info.getRefinementLevelReferenceProblem();
  // in general: for the homogenized case = 11 and for the high resolution case = 14
  // Note that this depends on the model problem!
#ifndef FINE_SCALE_REFERENCE
  refinement_level_referenceprob_ = 8;
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
  std :: string UnitCubeName( "../dune/multiscale/grids/cell_grids/unit_cube_0_centered.dgf" ); // --> the 0-centered unit cube, i.e. [-1/2,1/2]^2
  // note that the centering is fundamentaly important for the implementation. Do NOT change it to e.g. [0,1]^2!!!
  // to solve the cell problems, we always need a periodic gridPart.
  // Here it is always the unit cube that needs to be used (after transformation, cell problems are always formulated on such a grid )
  GridPtr< GridType > periodic_grid_pointer( UnitCubeName );
  periodic_grid_pointer->globalRefine( refinement_level_cellgrid );


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
               data_file << "Refinement Level for Periodic Micro Grid = " << refinement_level_cellgrid << std :: endl << std :: endl;
#ifdef TFR
               data_file << "We use TFR-HMM (HMM with test function reconstruction)." << std :: endl;
#else
               data_file << "We use HMM without test function reconstruction (NO TFR)." << std :: endl;
#endif
#ifdef AD_HOC_COMPUTATION
               data_file << "Cell problems are solved ad hoc (where required)." << std :: endl << std :: endl;
#else
               data_file << "Cell problems are solved and saved (in a pre-process)." << std :: endl << std :: endl;
               #ifdef ERRORESTIMATION
               data_file << "Error estimation activated!" << std :: endl << std :: endl;
               #endif
#endif
               data_file << "Epsilon = " << epsilon_ << std :: endl;
               data_file << "Estimated Epsilon = " << epsilon_est_ << std :: endl;
               data_file << "Delta (edge length of cell-cube) = " << delta_ << std :: endl;
#ifdef STOCHASTIC_PERTURBATION
               data_file << std :: endl << "Stochastic perturbation added. Variance = " << VARIANCE << std :: endl;
#endif
#ifdef ADAPTIVE
               data_file << std :: endl << "Adaptive computation. Global error tolerance for program abort = " << error_tolerance_ << std :: endl;
#endif
               data_file << std :: endl << std :: endl;
            }

     algorithm( UnitCubeName,
                macro_grid_pointer,
                fine_macro_grid_pointer,
                periodic_grid_pointer,
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
