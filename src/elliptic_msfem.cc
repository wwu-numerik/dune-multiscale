
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
// MsFEM currently only works for linear problems
#define LINEAR_PROBLEM

//! MsFEM in Petrov-Galerkin formulation (PGF) or standard MsFEM formulation?
#define PGF


//! is the homogenized solution available?
// this information should be provided by the 'problem specification file'
// there we define or don't define the macro HOMOGENIZEDSOL_AVAILABLE
// (if HOMOGENIZEDSOL_AVAILABLE == true, it means that it can be computed. It still needs to be determined by using a homogenizer )
#define HOMOGENIZEDSOL_AVAILABLE

// compute the L2 errors? (might be expensive)
#define ERROR_COMPUTATION

//! Do we have/want a fine-scale reference solution?
//#define FINE_SCALE_REFERENCE
#ifdef FINE_SCALE_REFERENCE

  // load the precomputed fine scale reference from a file
  // #define FSR_LOAD

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

//! Do we write the discrete MsFEM solution to a file? (for later usage)
#define WRITE_MSFEM_SOL_TO_FILE

//! Do we want to add a stochastic perturbation on the data?
//#define STOCHASTIC_PERTURBATION
#ifdef STOCHASTIC_PERTURBATION
 // size of variance:
 #define VARIANCE 0.01
 //! Do we want to force the algorithm to come to an end?
 // (was auch immer der Grund war, dass das Programm zuvor endlos lange weiter gelaufen ist. z.B. Tolerenzen nicht erreicht etc.)
 #define FORCE
#endif



//!NOTE: All the multiscale code requires an access to the 'ModelProblemData' class (typically defined in problem_specification.hh), which provides us with information about the problem itself.
// MsFEM Assembler, Error Estimator, ... they all hark back to 'ModelProblemData'. Probably there is a better solution, but for me, it works perfectly.





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

#if 1
//!---------
#include <dune/fem/solver/oemsolver/oemsolver.hh>


#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
//#include <dune/fem/io/visual/grape/datadisp/errordisplay.hh>
//!-----------
#endif

#include <dune/fem/solver/inverseoperators.hh>




//! local (dune-multiscale) includes
#include <dune/multiscale/problems/elliptic_problems/model_problem_6/problem_specification.hh>
//#include <dune/multiscale/problems/elliptic_problems/model_problem_easy/problem_specification.hh>

#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>

#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>

#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/deprecated/localproblemsolver.hh>

#include <dune/multiscale/tools/meanvalue.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/tools/assembler/matrix_assembler/deprecated/elliptic_msfem_matrix_assembler.hh>
#include <dune/multiscale/tools/misc/h1error.hh>

//! (very restrictive) homogenizer
#ifdef LINEAR_PROBLEM
#include <dune/multiscale/tools/homogenizer/elliptic_analytical_homogenizer.hh>
#include <dune/multiscale/tools/homogenizer/elliptic_homogenizer.hh>
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



//! ------------ cell problem solver and numbering manager -----------------------------------------

typedef LocalProblemNumberingManager< DiscreteFunctionSpaceType > LocProbNumberingManagerType;

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


//! --------------- the discrete operators (standard FEM and MsFEM) ------------------------

// discrete elliptic operator (corresponds with FEM Matrix)
typedef DiscreteEllipticOperator< DiscreteFunctionType, DiffusionType, MassTermType > EllipticOperatorType;

// discrete elliptic HMM operator (corresponds with HMM (or HMFEM) Matrix)
typedef DiscreteEllipticMsFEMOperator< DiscreteFunctionType, DiffusionType, LocProbNumberingManagerType > EllipticMsFEMOperatorType;

//! --------------------------------------------------------------------------------------







//! -------------------------------- important variables ---------------------------------

enum { dimension = GridType :: dimension};

//name of the error file in which the data will be saved
std :: string filename_;

double epsilon_; // 'epsilon' in for instance A^{epsilon}(t,x) = A(t,x/epsilon)
double epsilon_est_; // estimated epsilon in case epsilon is unknown
double delta_; //edge length of the cells in the cell proplems,

int refinement_level_macrogrid_;
int refinement_level_referenceprob_;

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
  H1Error< DiscreteFunctionType > h1error;

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

  #if 1
  AnalyticalHomogenizer< GridType, DiffusionType > ana_homogenizer( unit_cell_location );
  A_hom = ana_homogenizer.getHomTensor(diffusion_op);
  #endif

  // descretized homogenizer:

  #if 0
  Homogenizer< GridType, DiffusionType > disc_homogenizer( unit_cell_location );
  A_hom = disc_homogenizer.getHomTensor(diffusion_op);
  #endif

  // if known a-priori

  #if 0
  A_hom[0][0] = 1.0;
  A_hom[1][1] = 1.0;
  A_hom[0][1] = 0.0;
  A_hom[1][0] = 0.0;
  #endif
  //A_hom[0][0] = 0.667207;


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

 #endif

#endif


#ifdef FINE_SCALE_REFERENCE

 #ifdef FSR_COMPUTE
//! *******************************************************************

  // starting value for the Newton method
  DiscreteFunctionType zero_func( filename_ + " constant zero function ", finerDiscreteFunctionSpace );
  zero_func.clear();


  //! *************************** Assembling the reference problem ****************************
  // ( fine scale reference solution = fem_solution )

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
  DiscreteFunctionType fem_solution( filename_ + " Reference (FEM Newton) Solution", finerDiscreteFunctionSpace );
  fem_solution.clear();
  // By fem_solution, we denote the "fine scale reference solution" (used for comparison)
  // ( if the elliptic problem is linear, the 'fem_solution' is determined without the Newton method )

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
  fem_biCGStab( fem_newton_rhs, fem_solution );

  if (data_file.is_open())
    {
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
      data_file << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
    }

#endif

  //! ********************** End of assembling the reference problem ***************************
#endif
//end FSR_COMPUTE

#ifdef FSR_LOAD

  DiscreteFunctionType fem_solution( filename_ + " Reference (FEM Newton) Solution", finerDiscreteFunctionSpace );
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


  //! ************************* Assembling and solving the MsFEM problem ****************************

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

  std::cout << std :: endl << "Solving MsFEM-macro-problem for " << discreteFunctionSpace.size()
            << " unkowns and polynomial order "
            << DiscreteFunctionSpaceType :: polynomialOrder << "." 
            << std :: endl << std :: endl;

  //----------------------------------------------------------------------------------------------//
  //----------------------- THE DISCRETE HMM OPERATOR -----------------------------------//
  //----------------------------------------------------------------------------------------------//

  // to identify (macro) entities and basefunctions with a fixed global number, which stands for a certain local problem
  LocProbNumberingManagerType lp_num_manager( discreteFunctionSpace, refSimplexDiscreteFunctionSpace, filename_);


  //! define the elliptic hmm operator that describes our 'homogenized' macro problem
  // ( effect of the elliptic hmm operator on a certain discrete function )
  EllipticMsFEMOperatorType discrete_elliptic_msfem_op( discreteFunctionSpace, refSimplexDiscreteFunctionSpace, diffusion_op, lp_num_manager, filename_);

  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//

  //! matrix
  FEMMatrix msfem_matrix( "MsFEM stiffness matrix", discreteFunctionSpace, discreteFunctionSpace );

  //! right hand side vector
  // right hand side for the hm finite element method with Newton solver:
  DiscreteFunctionType msfem_rhs( "MsFEM right hand side", discreteFunctionSpace );
  msfem_rhs.clear();

  //! solution vector
  // solution of the heterogeneous multiscale finite element method, where we used the Newton method to solve the non-linear system of equations
  DiscreteFunctionType msfem_solution( filename_ + " MsFEM Solution", discreteFunctionSpace );
  msfem_solution.clear();

  // starting value for the Newton method
  DiscreteFunctionType zero_func_coarse( filename_ + " constant zero function coarse ", discreteFunctionSpace );
  zero_func_coarse.clear();

#ifdef LINEAR_PROBLEM

  // solve cell problems in a preprocess
  //! -------------- solve and save the cell problems for the base function set --------------------------------------

  MsFEMLocalProblemSolver< DiscreteFunctionType, DiffusionType > local_problem_solver(refSimplexDiscreteFunctionSpace, diffusion_op, data_file /*optinal*/);

  int number_of_grid_elements = refSimplexDiscreteFunctionSpace.grid().size(0);

  std :: cout << "Solving cell problems for " << number_of_grid_elements << " leaf entities." << std :: endl;

  // generate directory for cell problem data output
  if (mkdir(("data/MsFEM/" + filename_ + "/local_problems/").c_str() DIRMODUS) == -1)
   {
    std::cout << "WARNING! Directory for the solutions of the cell problems already exists!" << std :: endl;
   }
  else
   {
     mkdir(("data/MsFEM/" + filename_ + "/local_problems/").c_str() DIRMODUS);
   }

  // -------------- solve cell problems for the macro basefunction set ------------------------------
  // save the solutions of the cell problems for the set of macroscopic base functions

  local_problem_solver.saveTheSolutions_baseSet< DiscreteFunctionType >( discreteFunctionSpace, lp_num_manager, filename_ + "/local_problems/");

  // ------------- end solving and saving cell problems for the macro basefunction set --------------

  //! --------------- end solving and saving cell problems -----------------------------------------


  std :: cout << "Solving linear MsFEM problem." << std :: endl;
  if (data_file.is_open())
    {
      data_file << "Solving linear MsFEM problem." << std :: endl;
      data_file << "------------------------------------------------------------------------------" << std :: endl;
    }

  // to assemble the computational time
  Dune::Timer msfemAssembleTimer;

  // assemble the msfem stiffness matrix
  discrete_elliptic_msfem_op.assemble_matrix( msfem_matrix );
  // to print the matrix, use:   msfem_matrix.print();

  std::cout << "Time to assemble MsFEM macro stiffness matrix: " << msfemAssembleTimer.elapsed() << "s" << std::endl;
  if (data_file.is_open())
    {
      data_file << "Time to assemble MsFEM macro stiffness matrix: " << msfemAssembleTimer.elapsed() << "s" << std::endl;
    }

  // assemble right hand side
  rhsassembler.assemble_msfem< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( lp_num_manager, f , msfem_rhs);

  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
     boundaryTreatment( *it , msfem_rhs );

  InverseFEMMatrix msfem_biCGStab( msfem_matrix, 1e-8, 1e-8, 20000, VERBOSE );
  msfem_biCGStab( msfem_rhs, msfem_solution );

  std :: cout << "Linear MsFEM problem solved in " << msfemAssembleTimer.elapsed() << "s." << std :: endl << std :: endl;
  if (data_file.is_open())
    {
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
      data_file << "Linear MsFEM problem solved in " << msfemAssembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
    }

#endif




#ifdef WRITE_MSFEM_SOL_TO_FILE

  bool writer_is_open = false;

  char fname[40];
  sprintf( fname, "/msfem_solution_discFunc_refLevel_%d", refinement_level_macrogrid_ );
  std :: string fname_s( fname );

  std :: string location = "data/MsFEM/" + filename_ + fname_s;
  DiscreteFunctionWriter dfw( (location).c_str() );
  writer_is_open = dfw.open();
  if ( writer_is_open )
    dfw.append( msfem_solution );

#endif



#ifdef HOMOGENIZEDSOL_AVAILABLE

  bool hom_writer_is_open = false;

  char fhom_name[40];
  sprintf( fhom_name, "/homogenized_solution_discFunc_refLevel_%d", refinement_level_referenceprob_ );
  std :: string fhom_name_s( fhom_name );

  std :: string hom_location = "data/MsFEM/" + filename_ + fhom_name_s;
  DiscreteFunctionWriter hom_dfw( (hom_location).c_str() );
  hom_writer_is_open = hom_dfw.open();
  if ( hom_writer_is_open )
    hom_dfw.append( homogenized_solution );

#endif



#ifdef FINE_SCALE_REFERENCE

 #ifdef WRITE_FINESCALE_SOL_TO_FILE

  bool fine_writer_is_open = false;

  char fine_fname[50];
  sprintf( fine_fname, "/finescale_solution_discFunc_refLevel_%d", refinement_level_referenceprob_ );
  std :: string fine_fname_s( fine_fname );

  std :: string fine_location = "data/MsFEM/" + filename_ + fine_fname_s;
  DiscreteFunctionWriter fine_dfw( (fine_location).c_str() );
  fine_writer_is_open = fine_dfw.open();
  if ( fine_writer_is_open )
    fine_dfw.append( fem_solution );

 #endif

#endif

  //! ******************** End of assembling and solving the MsFEM problem ***************************

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


  RangeType h1_msfem_error(0.0);
  h1_msfem_error = h1error.semi_norm< ExactSolutionType >( u, msfem_solution );
  h1_msfem_error += exact_msfem_error;

  std :: cout << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std :: endl; }


 #ifdef FINE_SCALE_REFERENCE
  RangeType fem_error = l2error.norm< ExactSolutionType >( u, fem_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl; }
 #endif

 #ifdef HOMOGENIZEDSOL_AVAILABLE
  RangeType hom_exact_error = l2error.norm< ExactSolutionType >( u, homogenized_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_hom - u_exact ||_L2 =  " << hom_exact_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_hom - u_exact ||_L2 =  " << hom_exact_error << std :: endl; }
 #endif

#endif

// endif for macro ERROR_COMPUTATION
#endif



//! -------------------------------------------------------


//! --------------- writing data output ---------------------


  // general output parameters
  myDataOutputParameters outputparam;
  outputparam.set_path( "data/MsFEM/" + filename_ );

  // sequence stamp
  std::stringstream outstring;


  // --------- data output hmm solution --------------

  // create and initialize output class
  IOTupleType msfem_solution_series( &msfem_solution );
  #ifdef ADAPTIVE
  char msfem_prefix[30];
  sprintf( msfem_prefix, "msfem_solution_%d", loop_cycle );
  outputparam.set_prefix( msfem_prefix );
  #else
  outputparam.set_prefix("msfem_solution");
  #endif
  DataOutputType msfemsol_dataoutput( gridPart.grid(), msfem_solution_series, outputparam );

  // write data
  outstring << "msfem-solution";
  msfemsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
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
  IOTupleType fem_solution_series( &fem_solution );
  outputparam.set_prefix("reference_solution");
  DataOutputType fem_dataoutput( gridPartFine.grid(), fem_solution_series, outputparam );

  // write data
  outstring << "reference_solution";
  fem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
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


  std :: string save_filename = "data/MsFEM/" + filename_ + "/problem-info.txt";
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
  refinement_level_referenceprob_ = 12;
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
