
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

#include <dune/subgrid/subgrid.hh>


//! local (dune-multiscale) includes
#include <dune/multiscale/problems/elliptic_problems/model_problem_easy/problem_specification.hh>

#include <dune/multiscale/operators/righthandside_assembler.hh>

#include <dune/multiscale/operators/disc_func_writer/discretefunctionwriter.hh>

#include <dune/multiscale/operators/meanvalue.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/multiscale/operators/matrix_assembler/elliptic_fem_matrix_assembler.hh>



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


//! --------- typedefs for the local grid and the corresponding local ('sub') )discrete space -------------

typedef SubGrid< WORLDDIM , GridType > SubGridType; 

typedef LeafGridPart< SubGridType > SubGridPartType; 

typedef LagrangeDiscreteFunctionSpace < FunctionSpaceType, SubGridPartType, 1 > //1=POLORDER
   SubDiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction < SubDiscreteFunctionSpaceType > SubDiscreteFunctionType;

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





//! ---------  typedefs for the local function space (subgrid space) -------------

typedef SubGridType :: Codim<0> :: Entity SubgridEntityType; 
typedef SubGridType :: Codim<0> :: EntityPointer SubgridEntityPointerType; 
typedef SubGridType :: Codim<0> :: Geometry SubgridEntityGeometryType; 
typedef SubGridType :: Codim<1> :: Geometry SubgridFaceGeometryType;

typedef SubDiscreteFunctionSpaceType :: IteratorType SubgridIteratorType;

#if 0
typedef DiscreteFunctionSpaceType     :: BaseFunctionSetType      BaseFunctionSetType;

typedef CachingQuadrature < GridPartType , 0 > EntityQuadratureType;
typedef CachingQuadrature < GridPartType , 1 > FaceQuadratureType;

typedef DiscreteFunctionSpaceType :: DomainFieldType TimeType;
typedef DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;


typedef DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
typedef DiscreteFunctionType :: DofIteratorType DofIteratorType;
#endif

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




//! --------------- the matrix traits for the subgrid matrix -----------------------------

struct SubgridMatrixTraits
{
  typedef SubDiscreteFunctionSpaceType RowSubSpaceType;
  typedef SubDiscreteFunctionSpaceType ColumnSubSpaceType;
  typedef LagrangeMatrixSetup< false > StencilType;
  typedef ParallelScalarProduct< SubDiscreteFunctionSpaceType > SubgridParallelScalarProductType;

  template< class M >
  struct Adapter
  {
    typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
  };
};

//! --------------------------------------------------------------------------------------




//! --------------------- type of fem stiffness matrix -----------------------------------

typedef SparseRowMatrixOperator< DiscreteFunctionType, DiscreteFunctionType, MatrixTraits > FEMMatrix;

typedef SparseRowMatrixOperator< SubDiscreteFunctionType,
                                 SubDiscreteFunctionType,
			          SubgridMatrixTraits > SubgridFEMMatrix;

//! --------------------------------------------------------------------------------------


//! --------------- solver for the linear system of equations ----------------------------

// use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric ones. GMRES seems to be more stable, but is extremely slow!
typedef OEMBICGSQOp/*OEMBICGSTABOp*/< DiscreteFunctionType, FEMMatrix > InverseFEMMatrix;

typedef OEMBICGSQOp/*OEMBICGSTABOp*/< SubDiscreteFunctionType, SubgridFEMMatrix > InverseSubgridFEMMatrix;

//! --------------------------------------------------------------------------------------


//! --------------- the discrete operators (standard FEM) ------------------------

// discrete elliptic operator (corresponds with FEM Matrix)
typedef DiscreteEllipticOperator< DiscreteFunctionType, DiffusionType, MassTermType > EllipticOperatorType;

typedef DiscreteEllipticOperator< SubDiscreteFunctionType, DiffusionType, MassTermType > SubgridEllipticOperatorType;

//! --------------------------------------------------------------------------------------







//! -------------------------------- important variables ---------------------------------

enum { dimension = GridType :: dimension};

//name of the error file in which the data will be saved
std :: string filename_;

int refinement_level_macrogrid_;


//! -----------------------------------------------------------------------------



//! --------- typedefs and classes for data output -----------------------------------------

typedef Tuple<DiscreteFunctionType*> IOTupleType;
typedef DataOutput< GridType, IOTupleType> DataOutputType;


//! loeschen:
typedef Tuple< SubDiscreteFunctionType* > SubIOTupleType;
typedef DataOutput<SubGridType, SubIOTupleType> SubDataOutputType;





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







//! set the dirichlet points to zero
template< class EntityType, class HostDiscreteFunctionType, class DiscreteFunctionType >
void boundaryTreatment( const EntityType &entity, 
                        const typename DiscreteFunctionType :: DiscreteFunctionSpaceType :: GridType &sub_grid,
                        DiscreteFunctionType &rhs )
{
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
    LagrangePointSetType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: GridType GridType;


  typedef typename HostDiscreteFunctionType :: DiscreteFunctionSpaceType
    HostDiscreteFunctionSpaceType;

  typedef typename HostDiscreteFunctionSpaceType :: GridType HostGridType;

  typedef typename HostDiscreteFunctionSpaceType :: GridPartType HostGridPartType;

  typedef typename HostDiscreteFunctionSpaceType :: IteratorType :: Entity HostEntityType;

  typedef typename HostEntityType :: EntityPointer HostEntityPointerType; 

  HostEntityPointerType host_entity = sub_grid.template getHostEntity<0>( entity );


  enum { faceCodim = 1 };

  typedef typename GridPartType :: IntersectionIteratorType
    IntersectionIteratorType;

  typedef typename LagrangePointSetType :: template Codim< faceCodim >
                                        :: SubEntityIteratorType
    FaceDofIteratorType;

  const DiscreteFunctionSpaceType &discreteFunctionSpace = rhs.space();


  //! man muss das subgrid nicht uebergeben!!!
  const GridType &saaasdsasadsad = discreteFunctionSpace.grid();

  const GridPartType &gridPart = discreteFunctionSpace.gridPart();



#if 0

    LocalFunctionType rhsLocal = rhs.localFunction( entity );
    const LagrangePointSetType &lagrangePointSet
      = discreteFunctionSpace.lagrangePointSet( entity );

    const int face = 0;
    FaceDofIteratorType faceIterator
      = lagrangePointSet.template beginSubEntity< faceCodim >( face );
    const FaceDofIteratorType faceEndIterator
      = lagrangePointSet.template endSubEntity< faceCodim >( face );
    for( ; faceIterator != faceEndIterator; ++faceIterator )
      rhsLocal[ *faceIterator ] = 0;

#endif
#if 0
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
#endif


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



void algorithm ( GridPointerType &macro_grid_pointer, // grid pointer that belongs to the macro grid
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
  DiscreteExactSolutionType discrete_exact_solution( "discrete exact solution ", u, gridPart );
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
  EllipticOperatorType discrete_elliptic_op( discreteFunctionSpace, diffusion_op);

  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//


//! *******************************************************************

  // starting value for the Newton method
  DiscreteFunctionType zero_func( filename_ + " constant zero function ", discreteFunctionSpace );
  zero_func.clear();


  //! *************************** Assembling the reference problem ****************************
  // ( fine scale reference solution = fem_solution )

  //! (stiffness) matrix
  FEMMatrix fem_matrix( "FEM stiffness matrix", discreteFunctionSpace, discreteFunctionSpace );

  //! right hand side vector
  // right hand side for the finite element method:
  DiscreteFunctionType fem_rhs( "fem newton rhs", discreteFunctionSpace );
  fem_rhs.clear();

  //! solution vector
  // solution of the finite element method, where we used the Newton method to solve the non-linear system of equations
  // in general this will be an accurate approximation of the exact solution, that is why we it also called reference solution
  DiscreteFunctionType fem_solution( filename_ + " FEM Solution", discreteFunctionSpace );
  fem_solution.clear();
  // By fem_solution, we denote the "fine scale reference solution" (used for comparison)
  // ( if the elliptic problem is linear, the 'fem_solution' is determined without the Newton method )

  std :: cout << "Solving linear problem." << std :: endl;
  if (data_file.is_open())
    {
      data_file << "Solving linear problem with standard FEM and resolution level " << problem_info.getRefinementLevelReferenceProblem() << "." << std :: endl;
      data_file << "------------------------------------------------------------------------------" << std :: endl;
    }

  // to assemble the computational time
  Dune::Timer assembleTimer;

  // assemble the stiffness matrix
  discrete_elliptic_op.assemble_matrix( fem_matrix );

  std::cout << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
  if (data_file.is_open())
    {
      data_file << "Time to assemble standard FEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;
    }

  // assemble right hand side
  rhsassembler.assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 >( f , fem_rhs);

  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
         boundaryTreatment( *it , fem_rhs );


  InverseFEMMatrix fem_biCGStab( fem_matrix, 1e-8, 1e-8, 20000, VERBOSE );
  fem_biCGStab( fem_rhs, fem_solution );

  if (data_file.is_open())
    {
      data_file << "---------------------------------------------------------------------------------" << std :: endl;
      data_file << "Standard FEM problem solved in " << assembleTimer.elapsed() << "s." << std :: endl << std :: endl << std :: endl;
    }




  //! ********************** End of assembling the standard fem problem ***************************






















//! BATTLE FIELD
#if 1


   typedef GridType :: Codim< 0 > :: Partition< All_Partition > :: LevelIterator LevelEntityIteratorType;
   typedef GridPartType::IndexSetType IndexSetType;


   std :: vector< bool > cell_mark;
   

   const int codim = 0;
   const int gridlevel = 0;
   const int maxlevel = grid.maxLevel();

   LevelEntityIteratorType level_0_iterator_end = grid.lend< codim >( gridlevel );
   LevelEntityIteratorType level_0_iterator_begin = grid.lbegin< codim >( gridlevel );

   LevelEntityIteratorType a_level_0_entity = level_0_iterator_begin;
   ++a_level_0_entity;

//   for( ; level_0_iterator_begin != level_0_iterator_end; ++level_0_iterator_begin )


#if 1
   
   // create subgrid:
   SubGrid< dimension , GridType > subGrid(grid);
   subGrid.createBegin();

   subGrid.insert( *a_level_0_entity );
   ++a_level_0_entity;
   subGrid.insert( *a_level_0_entity );
   ++a_level_0_entity;
   subGrid.insert( *a_level_0_entity );
   ++a_level_0_entity;
   subGrid.insert( *a_level_0_entity );
   subGrid.createEnd();

   const EntityType& host_entity = *level_0_iterator_begin;//*a_level_0_entity;
   std :: cout << "subGrid.contains( host_entity ) = " << subGrid.contains<0>( host_entity ) << std :: endl;



   // subGrid.report();
    
   SubGridPartType subGridPart( subGrid );

   SubDiscreteFunctionSpaceType subDiscreteFunctionSpace( subGridPart );

   SubDiscreteFunctionType local_solution( filename_ + " Sub FEM Solution", subDiscreteFunctionSpace );
   local_solution.clear();

   SubgridIteratorType sub_endit = subDiscreteFunctionSpace.end();
   for( SubgridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
      {
        std :: cout << "subGridPart.indexSet().index( *sub_it ) = " << subGridPart.indexSet().index( *sub_it ) << std :: endl;

        EntityPointerType host_entity_2 = subGrid.getHostEntity<0>( *sub_it );

      }


#if 0
   typedef SubGridType :: Codim< 0 > :: Partition< All_Partition > :: LevelIterator SubgridLevelEntityIteratorType;
   SubgridLevelEntityIteratorType sub_max_level_iterator_end = subGrid.lend< 1 >( 0 );
   SubgridLevelEntityIteratorType sub_max_level_iterator_begin = subGrid.lbegin< 1 >( 0 );
   for( SubgridLevelEntityIteratorType max_level_it = sub_max_level_iterator_begin ; max_level_it != sub_max_level_iterator_end; ++max_level_it )
        std :: cout << " " << std :: endl;
#endif

#endif

   int id = 0;

   LevelEntityIteratorType max_level_iterator_end = grid.lend< codim >( maxlevel );
   LevelEntityIteratorType max_level_iterator_begin = grid.lbegin< codim >( maxlevel );
   for( LevelEntityIteratorType max_level_it = max_level_iterator_begin ; max_level_it != max_level_iterator_end; ++max_level_it )
     {

       //std :: cout << "indexSet.index( *max_level_it ) = " << indexSet.index( *max_level_it ) << std:: endl;

        std :: cout << "gridPart.indexSet().index( *max_level_it ) = " << gridPart.indexSet().index( *max_level_it ) << std :: endl;
        std :: cout << "id = " << id << std :: endl;

        EntityPointerType fine_father_entity = max_level_it;
        for (int lev = 0; lev < maxlevel; ++lev)
            fine_father_entity = fine_father_entity->father();

        if ( fine_father_entity == a_level_0_entity )
         {
           cell_mark.push_back(true);
           std :: cout << "true" << std :: endl;
         }
        else
         {
           cell_mark.push_back(false);
           std :: cout << "false" << std :: endl;
         }

       id += 1;

     }




   std :: cout << "Grid MaxLevel = " << grid.maxLevel() << std :: endl;


   int size = 0;
   //grid.size();
   std :: cout << "Grid Size = " << grid.size( gridlevel, codim ) << std :: endl;




  RightHandSideAssembler< SubDiscreteFunctionType > subgridrhsassembler;
  SubgridEllipticOperatorType sub_discrete_elliptic_op( subDiscreteFunctionSpace, diffusion_op);
  SubDiscreteFunctionType sub_zero_func( filename_ + " constant zero function ", subDiscreteFunctionSpace );
  sub_zero_func.clear();
  SubgridFEMMatrix subgrid_fem_matrix( "Dubgrid FEM stiffness matrix", subDiscreteFunctionSpace, subDiscreteFunctionSpace );

  SubDiscreteFunctionType sub_fem_rhs( "subgrid fem rhs", subDiscreteFunctionSpace );
  sub_fem_rhs.clear();


  sub_discrete_elliptic_op.assemble_matrix( subgrid_fem_matrix , false );


  subgridrhsassembler.assemble< 2 * SubDiscreteFunctionSpaceType :: polynomialOrder + 2 >( f , sub_fem_rhs);



#if 1
   for( SubgridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
       {
          boundaryTreatment< SubgridEntityType, DiscreteFunctionType, SubDiscreteFunctionType >( *sub_it , subGrid, sub_fem_rhs );
       }
#endif

  InverseSubgridFEMMatrix sub_fem_biCGStab( subgrid_fem_matrix, 1e-8, 1e-8, 20000, VERBOSE );
  //sub_fem_biCGStab( sub_fem_rhs, local_solution );












#endif































  //! -- write discrete solution to file ---

  bool writer_is_open = false;

  char fem_fname[50];
  sprintf( fem_fname, "/fem_solution_discFunc_refLevel_%d", refinement_level_macrogrid_ );
  std :: string fem_fname_s( fem_fname );

  std :: string fine_location = "data/MsFEM/" + filename_ + fem_fname_s;
  DiscreteFunctionWriter fem_dfw( (fine_location).c_str() );
  writer_is_open = fem_dfw.open();
  if ( writer_is_open )
    fem_dfw.append( fem_solution );

  //! -----------------------------


  std :: cout << std :: endl << "The L2 errors:" << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "The L2 errors:" << std :: endl << std :: endl; }



  //! ----------------- compute L2-errors -------------------
#if 0
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




// endif for macro ERROR_COMPUTATION
#endif

#ifdef EXACTSOLUTION_AVAILABLE

  RangeType fem_error = l2error.norm< ExactSolutionType >( u, fem_solution, 2 * DiscreteFunctionSpaceType :: polynomialOrder + 2 );

  std :: cout << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl << std :: endl;
  if (data_file.is_open())
    { data_file << "|| u_fem - u_exact ||_L2 =  " << fem_error << std :: endl; }

#endif

//! -------------------------------------------------------


//! --------------- writing data output ---------------------


  // general output parameters
  myDataOutputParameters outputparam;
  outputparam.set_path( "data/MsFEM/" + filename_ );

  // sequence stamp
  std::stringstream outstring;

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
  
  // --------- data output standard solution --------------

  // create and initialize output class
  IOTupleType fem_solution_series( &fem_solution );
  outputparam.set_prefix("fem_solution");
  DataOutputType fem_dataoutput( gridPart.grid(), fem_solution_series, outputparam );

  // write data
  outstring << "fem_solution";
  fem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------
  

#if 0
  // --------- data output local solution --------------

  // create and initialize output class
  SubIOTupleType sub_fem_solution_series( &local_solution );
  outputparam.set_prefix("local_solution");
 
  SubDataOutputType sub_fem_dataoutput( subGrid, sub_fem_solution_series, outputparam );

  // write data
  outstring << "local_solution";
  sub_fem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
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

  // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
  refinement_level_macrogrid_ = atoi( argv[ 1 ] );

  //name of the grid file that describes the macro-grid:
  std :: string macroGridName;
  info.getMacroGridFile( macroGridName );
  std :: cout << "loading dgf: " << macroGridName << std :: endl;

  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values for the parameters:

  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer( macroGridName );
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine( refinement_level_macrogrid_ );


  // to save all information in a file
  std :: ofstream data_file( (save_filename).c_str() );
  if (data_file.is_open())
            {
               data_file << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << "." << std :: endl << std :: endl;
               data_file << "Computations were made for:" << std :: endl << std :: endl;
               data_file << "Refinement Level for (uniform) Macro Grid = " << refinement_level_macrogrid_ << std :: endl;
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
