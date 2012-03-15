
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


//#define USE_GRAPE HAVE_GRAPE
#ifndef USE_GRAPE
#define USE_GRAPE HAVE_GRAPE
#endif

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>


// polynomial order of discrete space
#define POLORDER 1

// grid type
#define GRIDTYPE ALBERTAGRID
//#define GRIDTYPE YASPGRID

// computational domain is a subset of \R^{GRIDDIM}
#define GRIDDIM 2
#define WORLDDIM GRIDDIM

#ifndef POLORDER
  #define POLORDER 1
#endif

//- system includes
#include <iostream>
#include <sstream>

#include <stdio.h>
//#include <cstdio>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>



#define USE_TWISTFREE_MAPPER
#define VERBOSE false





//-----------------------------



#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions


// for alberta grid:
#include <dune/grid/albertagrid/dgfparser.hh>
// for yasp grid
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#if USE_GRAPE
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

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>


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

#include <dune/fem/misc/mpimanager.hh>

//- local includes
#include <dune/multiscale/space/ownbasisspace/examples/sinespace.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

using namespace Dune;

//! check for GridType:
#if GRIDTYPE==ALBERTAGRID
 typedef AlbertaGrid< GRIDDIM, WORLDDIM > GridType;
#endif

#if GRIDTYPE==YASPGRID
  // the Grid Type ( Yasp-Grid )
 typedef YaspGrid< GRIDDIM > GridType;
#endif

const int polOrder = POLORDER;






//! --------- typedefs for ... -------------

typedef AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;

typedef GridPtr< GridType > GridPointerType;

typedef FunctionSpace < double , double , WORLDDIM , 1 > FunctionSpaceType;

typedef LagrangeDiscreteFunctionSpace < FunctionSpaceType, GridPartType, 1 > //1=POLORDER
   DiscreteFEMSpaceType;

typedef AdaptiveDiscreteFunction< DiscreteFEMSpaceType >        FEMFunctionType;

//!-----------------------------------------------------------------------------------------

typedef SineReducedBasisSpace< DiscreteFEMSpaceType, 4 >    AdvancedDiscreteFunctionSpaceType;
typedef AdaptiveDiscreteFunction< AdvancedDiscreteFunctionSpaceType >        AdvancedDiscreteFunctionType;


//! --------- typedefs and classes for data output -----------------------------------------

typedef Tuple<AdvancedDiscreteFunctionType*> IOTupleType;
typedef DataOutput<GridType, IOTupleType> DataOutputType;

typedef Tuple<FEMFunctionType*> FEMIOTupleType;
typedef DataOutput<GridType, FEMIOTupleType> FEMDataOutputType;


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
        return "data_output_l2projection";
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




template< class DiscreteFunctionImp >
class L2Projection
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;

private:
  typedef L2Projection< DiscreteFunctionType > ThisType;

public:
  typedef typename DiscreteFunctionType
            :: DiscreteFunctionSpaceType                             DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType :: DomainType           DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeType            RangeType;

  typedef typename DiscreteFunctionSpaceType :: DomainFieldType      DomainFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType       RangeFieldType;

public:
  template< class FunctionType >
  static inline void project ( const FunctionType &function,
                               DiscreteFunctionType &discreteFunction )
  {
    typedef typename DiscreteFunctionType :: LocalFunctionType       LocalFunctionType;

    typedef typename LocalFunctionType :: BaseFunctionSetType        BaseFunctionSetType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType       IteratorType;

    typedef typename GridPartType :: GridType :: template Codim< 0 > :: Entity :: Geometry
      GeometryType;

    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef typename QuadratureType :: CoordinateType QuadraturePointType;

    const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();

    discreteFunction.clear();

    const IteratorType end = dfSpace.end();
    for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      LocalFunctionType localFunction = discreteFunction.localFunction( *it );

      const BaseFunctionSetType &baseFunctionSet = localFunction.baseFunctionSet();
      const unsigned int numBaseFunctions = baseFunctionSet.numBaseFunctions();

      QuadratureType quadrature( *it, 2*dfSpace.order() + 2);
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {

        const QuadraturePointType &point = quadrature.point( pt );

        const GeometryType &geometry = it->geometry();

        RangeFieldType weight = quadrature.weight( pt ) * geometry.integrationElement( point );

        RangeType y;
        function.evaluate( geometry.global( point ), y );

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {

          RangeType phi = 0.0;
          baseFunctionSet.evaluate( i, quadrature[ pt ], phi );
          localFunction[ i ] += weight * (y * phi);
        }
      }
    }
  }
};



template< class DiscreteFunctionImp >
class L2ErrorInternal
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;


private:
  typedef L2Error< DiscreteFunctionType > ThisType;

public:
  typedef typename DiscreteFunctionType
            :: DiscreteFunctionSpaceType                             DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: DomainType           DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeType            RangeType;
  typedef typename DiscreteFunctionSpaceType :: DomainFieldType      DomainFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType       RangeFieldType;

  static const unsigned int DimDomain = DiscreteFunctionSpaceType :: DimDomain;
  static const unsigned int DimRange = DiscreteFunctionSpaceType :: DimRange;

public:
  template< class FunctionType >
  static inline void norm ( const FunctionType &function,
                            const DiscreteFunctionType &discreteFunction,
                            RangeType &error )
  {
    typedef typename DiscreteFunctionType :: LocalFunctionType       LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType       IteratorType;

    typedef typename GridPartType :: GridType
              :: template Codim< 0 > :: Entity :: Geometry           GeometryType;

    typedef CachingQuadrature< GridPartType, 0 >                     QuadratureType;
    typedef typename QuadratureType :: CoordinateType                QuadraturePointType;

    const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();

    error = 0;

    const IteratorType end = dfSpace.end();
    for( IteratorType it = dfSpace.begin(); it != end ; ++it )
    {
      LocalFunctionType localFunction = discreteFunction.localFunction( *it );

      QuadratureType quadrature( *it, 2*dfSpace.order() + 2 );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const QuadraturePointType &point = quadrature.point( pt );

        const GeometryType &geometry = it->geometry();

        RangeFieldType weight = quadrature.weight( pt ) * geometry.integrationElement( point );

        RangeType y;
        function.evaluate( geometry.global( point ), y );

        RangeType phi = 0.0;
        localFunction.evaluate( quadrature[ pt ], phi );
        for( unsigned int i = 0; i < DimRange; ++i )
          error[ i ] += weight * SQR( y[ i ] - phi[ i ] );
      }
    }

    for( unsigned int i = 0; i < DimRange; ++i)
      error[ i ] = sqrt( error[ i ] );
  }
};


void algorithm ( GridPartType &gridPart )
{

  DiscreteFEMSpaceType femSpace( gridPart );

//std :: cout << "..." << std :: endl;

  AdvancedDiscreteFunctionSpaceType advancedDiscFunctionSpace( femSpace );

  int numBaseFuncs = advancedDiscFunctionSpace.numBaseFunctions ();
  std :: cout << "Number of basis functions = " << numBaseFuncs << std :: endl;

  typedef AdvancedDiscreteFunctionSpaceType::BaseFunctionType BaseFunctionType;

#if 1

  FEMFunctionType some_fem_function( "some fem function", femSpace );
  some_fem_function.clear();

  FEMFunctionType base_function_i( "base function i", femSpace );
  for( unsigned int i = 0; i < numBaseFuncs; ++i)
    {

      advancedDiscFunctionSpace.getBaseFunction( i, base_function_i );

      double coeff_i = (1.0 / (i+1) );
      //std :: cout << "coeff_i = " << coeff_i << std :: endl;

      base_function_i *= coeff_i;
      some_fem_function += base_function_i;

    }

  typedef AdvancedDiscreteFunctionSpaceType::RangeType RangeType;
  typedef AdvancedDiscreteFunctionSpaceType::DomainType DomainType;

#if 0
  const int BlockSize = 10;
  const int N = BlockSize;

  typedef FieldVector< double, BlockSize > VectorType;
  typedef FieldMatrix< double, BlockSize, BlockSize > MatrixType;


  // build little blocks
  // if 'BlockSize = 1':
  //   diagonal_block = 4
  MatrixType A;
  for (int i = 0; i < BlockSize; i++ )
    for (int j = 0; j < BlockSize; j++ )
      {
        if ( i == j )
          { A[i][j] = 4; }
        else
          { A[i][j] = 0; }
      }


#endif

#if 1
  // define Types
  const int BlockSize = 1;

  typedef FieldVector< double, BlockSize > VectorBlockType;
  typedef FieldMatrix< double, BlockSize, BlockSize > MatrixBlockType;
  typedef BlockVector< VectorBlockType > VectorType;
  typedef BCRSMatrix< MatrixBlockType > MatrixType;

  // build little blocks
  // if 'BlockSize = 1':
  //   diagonal_block = 4
  MatrixBlockType diagonal_block = 0;
  for (int i = 0; i < BlockSize; i++ )
    for (int j = 0; j < BlockSize; j++ )
      {
        if ( i == j )
          { diagonal_block[i][j] = 4 + ( BlockSize - 1 ); }
        else
          { diagonal_block[i][j] = -1; }
      }

  // if 'BlockSize = 1':
  //   subdiagonal_block = -1
  MatrixBlockType subdiagonal_block = 0;
  for ( int i = 0; i < BlockSize; i++ )
    subdiagonal_block[i][i] = -1;


  // make a block compressed row matrix with five point stencil
  const int BW2 = 31;
  const int N = BW2 * BW2;
  // matrix (rows, columns, number of 'non-zero' entries, build mode )
  //   (use 'build mode = row_wise' to get a sequential order )
//!  MatrixType A ( N, N, 5*N, Dune::BCRSMatrix< MatrixBlockType >::row_wise );
  MatrixType A ( BW2, BW2, Dune::BCRSMatrix< MatrixBlockType >::row_wise );

#if 1
  // iterator ueber die rows:
  for (MatrixType::CreateIterator i = A.createbegin(); i != A.createend(); ++i )
	{

#if 0
          // number of the row
          int row = i.index()/BW2;

          // number of the column
          int col = i.index()%BW2;
#endif

          // number of the row
          int row = i.index();

          // number of the column
          int col = i.index();
          i.insert( col );

          std :: cout << "i.index() vorher = " << i.index() << std :: endl;
          std :: cout << "col = " << col << std :: endl;
          std :: cout << "row = " << row << std :: endl << std :: endl;

	  i.insert(i.index());
#if 1
          // insert = put column index in row 

	  if ( col - 1 >= 0) i.insert( i.index() - 1 );

	  if ( col + 1 < BW2) i.insert( i.index()+1 );
	  if ( row-1>= 0 ) i.insert( i.index()-BW2 );
	  if ( row+1< BW2 ) i.insert( i.index()+BW2 );
#endif
          //std :: cout << "i.index() nachher = " << i.index() << std :: endl;

	}
#endif
abort();

  for ( MatrixType::RowIterator i = A.begin(); i != A.end(); ++i )
    for ( MatrixType::ColIterator j = (*i).begin(); j != (*i).end(); ++j )
      if ( i.index() == j.index() )
        (*j) = diagonal_block;
      else
        (*j) = subdiagonal_block;

#endif


  //  printmatrix(std::cout,A,"system matrix","row",10,2);

  // set up system
  VectorType x(N), b(N);

  // prescribe known solution:
  x = 0;
  x[0] = 1;
  x[N-1] = 2;

  // set right hand side accordingly:
  b = 0;
  A.umv(x,b);

  // initial guess
  x = 1;
  for (int i = 0; i < N; i++ )
    x[i] = i * 0.1;

  // set up the high-level solver objects
  Dune::MatrixAdapter< MatrixType, VectorType, VectorType> op(A);        // make linear operator from A

  Dune::SeqJac  < MatrixType, VectorType, VectorType > jac(A,1,1);         // Jacobi preconditioner
  Dune::SeqGS   < MatrixType, VectorType, VectorType > gs(A,1,1);          // GS preconditioner
  Dune::SeqSOR  < MatrixType, VectorType, VectorType > sor(A,1,1.9520932); // SSOR preconditioner
  Dune::SeqSSOR < MatrixType, VectorType, VectorType > ssor(A,1,1.0);      // SSOR preconditioner
  Dune::SeqILU0 < MatrixType, VectorType, VectorType > ilu0(A,1.0);        // preconditioner object
  Dune::SeqILUn < MatrixType, VectorType, VectorType > ilu1(A,1,0.92);     // preconditioner object

  Dune::LoopSolver< VectorType > loop( op, jac, 1E-4, 18000, 2 );      // an inverse operator 
  Dune::CGSolver< VectorType > cg( op, ilu0, 1E-4, 8000, 2 );          // an inverse operator 
  Dune::BiCGSTABSolver< VectorType > bcgs( op, ilu1, 1E-8, 8000, 2 );  // an inverse operator 
  Dune::GradientSolver< VectorType > gras( op, jac, 1E-4, 18000, 2 );  // an inverse operator 

  // call the solver
  Dune::InverseOperatorResult r;
  loop.apply(x,b,r);

#if 0
  typedef Matrix< AdvancedDiscreteFunctionSpaceType::RangeType > MatrixType;

  MatrixType stiffness_matrix( numBaseFuncs, numBaseFuncs );
  for( unsigned int i = 0; i < numBaseFuncs; ++i)
    {
       for( unsigned int j = 0; j < numBaseFuncs; ++j)
         {
           stiffness_matrix[ i ][ j ] = 1.0;
         }
    }

  for( unsigned int i = 0; i < numBaseFuncs; ++i)
    {
       for( unsigned int j = 0; j < numBaseFuncs; ++j)
         {
           std :: cout << stiffness_matrix[ i ][ j ] << " ";
         }
       std :: cout << std :: endl;
    }

   typedef AssembledLinearOperator< MatrixType, DomainType, RangeType > LinearOperatorType;
#endif




// Naechster Schritt: LGS mit istl-Matrizen aufstellen, loesen und schauen ob das Gleiche rauskommt, wie mit dem standard-Loeser.


#endif

  ExactFunction< FunctionSpaceType > someContinuousFunction;
  // continuous function in new basis space
  AdvancedDiscreteFunctionType discContinuousFunction( "discrete function", advancedDiscFunctionSpace );


  L2Projection< AdvancedDiscreteFunctionType > :: project( someContinuousFunction, discContinuousFunction );

  FunctionSpaceType :: RangeType error = 0.0;
  L2ErrorInternal< AdvancedDiscreteFunctionType > :: norm( someContinuousFunction, discContinuousFunction, error );
  std :: cout << "L2 Error: " << error << std :: endl;






#if 1

//! --------------- writing data output ---------------------

  // general output parameters
  myDataOutputParameters outputparam;
  outputparam.set_path( "data/l2projection/test" );

  // sequence stamp
  std::stringstream outstring;


  // --------- data output for an advancved discrete function --------------

  // create and initialize output class
  IOTupleType l2projection_solution_series( &discContinuousFunction );
  outputparam.set_prefix("advanced_discrete_function");
  DataOutputType l2projection_dataoutput( gridPart.grid(), l2projection_solution_series, outputparam );

  // write data
  outstring << "advanced_discrete_function";
  l2projection_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------


  // --------- data output for a fem function --------------

  // create and initialize output class
  FEMIOTupleType l2projection_fem_solution_series( &some_fem_function );
  outputparam.set_prefix("fem_discrete_function");
  FEMDataOutputType l2projection_fem_dataoutput( gridPart.grid(), l2projection_fem_solution_series, outputparam );

  // write data
  outstring << "fem_discrete_function";
  l2projection_fem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str(std::string());

  // -------------------------------------------------------


//!-------------------------------------------------------------

#endif


  #if USE_GRAPE
    GrapeDataDisplay< GridSelector :: GridType > grape( gridPart );
    grape.dataDisplay( solution );
  #endif

}



int main ( int argc, char **argv )
{

  MPIManager :: initialize( argc, argv );
  if( argc != 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
    return 1;
  }


  try
  {

    unsigned int level = atoi( argv[ 1 ] );
    double error = 0.0;

    const unsigned int step = DGFGridInfo< GridType > :: refineStepsForHalf();

    std :: ostringstream macroGridNameStream;
    macroGridNameStream << "../dune/multiscale/space/ownbasisspace/examples/2dgrid.dgf";
    std :: string macroGridName = macroGridNameStream.str();

    GridPtr< GridType > gridptr( macroGridName );
    GridType &grid = *gridptr;
    GridPartType gridPart( grid );

    grid.globalRefine( level );
    algorithm( gridPart);

    system("rm -r -f ReducedSpace_BaseFunctions/*");
    system("rm -r ReducedSpace_BaseFunctions");

    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }


}
