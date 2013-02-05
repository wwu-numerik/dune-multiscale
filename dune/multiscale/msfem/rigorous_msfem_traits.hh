#ifndef RIGOROUS_MSFEM_TRAITS_HH
#define RIGOROUS_MSFEM_TRAITS_HH

#include <dune/fem/space/reducedbasisspace/reducedbasisspace.hh>

//#include <dune/multiscale/tools/errorestimation/MsFEM/msfem_elliptic_error_estimator.hh>
//#include <dune/multiscale/tools/solver/MsFEM/msfem_solver.hh>
//#include <dune/subgrid/subgrid.hh>
//#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/fem/io/file/dataoutput.hh>

namespace Dune {

template < class T >
class MacroMicroGridSpecifier;


#if 1
  template< class FunctionSpaceImp >
  class ExactFunction
  : public Fem::Function< FunctionSpaceImp, ExactFunction< FunctionSpaceImp > >
  {
  private:
    typedef ExactFunction< FunctionSpaceImp > ThisType;
    typedef Function< FunctionSpaceImp, ThisType > BaseType;

  public:
    typedef FunctionSpaceImp                                        FunctionSpaceType;

  public:
    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    typedef typename FunctionSpaceType :: RangeType                  RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType            DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType             RangeFieldType;

    static const unsigned int dimDomain = FunctionSpaceType :: dimDomain;
    static const unsigned int dimRange = FunctionSpaceType :: dimRange;

  public:

    inline void evaluate ( const DomainType &x, RangeType &y ) const
    {
      y = 1;
      for( unsigned int i = 0; i < dimDomain; ++i )
       {
         const DomainFieldType &xi = x[ i ];
         y *= xi - xi * xi;
       }
    }

    inline void evaluate ( const DomainType &x, const RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }

  };

  template< class FunctionSpaceImp >
  class SineBaseFunction
  : public Fem::Function< FunctionSpaceImp, SineBaseFunction< FunctionSpaceImp > >
  {
  private:
    typedef SineBaseFunction< FunctionSpaceImp >                     ThisType;
    typedef Function< FunctionSpaceImp, ThisType > BaseType;

  public:
    typedef FunctionSpaceImp                                         FunctionSpaceType;

  public:
    typedef typename FunctionSpaceType :: DomainType                 DomainType;
    typedef typename FunctionSpaceType :: RangeType                  RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType            DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType             RangeFieldType;

    static const unsigned int dimDomain = FunctionSpaceType :: dimDomain;
    static const unsigned int dimRange = FunctionSpaceType :: dimRange;

    typedef FieldVector< int, dimDomain > CoefficientType;
  public:
    explicit SineBaseFunction ( const CoefficientType coefficient )
    : coefficient_( coefficient )
    {
      base_func_counter_ += 1;
      base_func_id_ = base_func_counter_;
      // std :: cout << "base_func_id_ = " << base_func_id_ << std :: endl;
    }

    inline void evaluate ( const DomainType &x, RangeType &y ) const
    {

      y = 1;

      for( unsigned int i = 0; i < dimDomain; ++i )
      {
        y *= sqrt( 2 ) * sin( M_PI * coefficient_[ i ] * x[ i ] );
      }

      // the fifth basis function - counting starts with zero
      if ( (coefficient_[0]==2) && (coefficient_[1]==2) )
        { y = 1; }

    }

    inline void evaluate ( const DomainType &x, const RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }

  protected:
    const CoefficientType coefficient_;
    int base_func_id_;

    static int base_func_counter_;
  };

  template< class FunctionSpaceImp >
  int SineBaseFunction< FunctionSpaceImp > :: base_func_counter_ = -1;


  // number of basis functions = maxCoefficient^worlddim
  template< class FEMSpaceImp, unsigned int maxCoefficient >
  class SineReducedBasisSpace
  : public Dune::ReducedBasisSpace< AdaptiveDiscreteFunction< FEMSpaceImp > >
  {
  private:
    typedef SineReducedBasisSpace< FEMSpaceImp,
                                   maxCoefficient >                  ThisType;
    typedef Dune :: ReducedBasisSpace
                    < AdaptiveDiscreteFunction
                     < FEMSpaceImp > >                      BaseType;
  public:
    typedef FEMSpaceImp                                     FEMSpaceType;

    typedef AdaptiveDiscreteFunction< FEMSpaceType >        FEMFunctionType;

    typedef SineBaseFunction< typename FEMSpaceType
                             :: FunctionSpaceType >                  ContinuousBaseFunctionType;

  private:
    typedef typename ContinuousBaseFunctionType :: CoefficientType   CoefficientType;

  public:
    inline explicit SineReducedBasisSpace ( FEMSpaceType &FEMSpace )
    : BaseType( FEMSpace )
    {
      // CoefficientType coefficient( -maxCoefficient );
      CoefficientType coefficient( 1.0 );

      for( unsigned int i = 0; i < maxCoefficient; ++i )
       {
         for( unsigned int j = 0; j < maxCoefficient; ++j )
           {
             coefficient[ 0 ] = i+1;
             coefficient[ 1 ] = j+1;
             addBaseFunction( coefficient );
             // std :: cout << "coefficient = " << coefficient << std :: endl;
           }
       }
    }

  private:
    static inline int abs( const CoefficientType &coefficient )
    {
      int value = 0;
      for( unsigned int i = 0; i < CoefficientType :: dimension; ++i )
        value += (coefficient[ i ] < 0 ? -coefficient[ i ] : coefficient[ i ]);
      return value;
    }

    inline void addBaseFunction( const CoefficientType &coefficient )
    {
      FEMFunctionType discreteBaseFunction( "base function", baseFunctionSpace_ );
      ContinuousBaseFunctionType continuousBaseFunction( coefficient );
      LagrangeInterpolation< FEMFunctionType >
        :: interpolateFunction( continuousBaseFunction, discreteBaseFunction );
      BaseType :: addBaseFunction( discreteBaseFunction );
    }


  public:

#if 1
    inline void getBaseFunction( unsigned int i, FEMFunctionType &base_func ) const
    {

      BaseType :: baseFunction( i, base_func );

    }
#endif

    using BaseType :: baseFunctionSpace_;
  };


#endif


//! --------- typedefs for the coefficient and data functions ------------------
struct RigorousMsfemTraits {
    
  //! ----- typedefs for the macro grid and the corresponding discrete space -----
  typedef Dune::GridSelector::GridType
  GridType;
  // Dune::InteriorBorder_Partition or Dune::All_Partition >?
  // see:
  // http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
  typedef Dune::AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;

  typedef Dune::GridPtr< GridType > GridPointerType;

  typedef Dune::FunctionSpace< double, double, WORLDDIM, 1 > FunctionSpaceType;
  // type of first source term (right hand side of differential equation or type of 'f')
  typedef Problem::FirstSource< FunctionSpaceType > FirstSourceType;

  // type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
  typedef Problem::Diffusion< FunctionSpaceType > DiffusionType;

  // default type for any missing coefficient function (e.g. advection,...)
  typedef Problem::DefaultDummyFunction< FunctionSpaceType > DefaultDummyFunctionType;

  // type of exact solution (in general unknown)
  typedef Problem::ExactSolution< FunctionSpaceType > ExactSolutionType;
  typedef Dune::GridFunctionAdapter< ExactSolutionType, GridPartType >
    DiscreteExactSolutionType;     // for data output with paraview or grape


  //! ----  typedefs for the standard discrete function space (macroscopic) -----
  typedef FunctionSpaceType::DomainType DomainType;
  //! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
  typedef FunctionSpaceType::RangeType RangeType;
  typedef std::vector< RangeType > RangeVector;
  typedef std::vector< RangeVector > RangeVectorVector;
  //! defines the function space to which the numerical solution belongs to
  //! see dune/fem/lagrangebase.hh
  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER >  // 1 = POLORDER
    DiscreteFunctionSpaceType;

  typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  //!-----------------------------------------------------------------------------

  //!------------------------- for adaptive grid refinement ---------------------------------
  //! For adaption:
  //! type of restrict-prolong operator
  typedef Dune::RestrictProlongDefault< DiscreteFunctionType >
    RestrictProlongOperatorType;
  //! type of the adaption manager
  typedef Dune::AdaptationManager< GridType, RestrictProlongOperatorType >
    AdaptationManagerType;
  //!---------------------------------------------------------------------------------------

#if 0
  typedef Dune::MacroMicroGridSpecifier< DiscreteFunctionSpaceType >                          MacroMicroGridSpecifierType;
  typedef Dune::SubGrid< GridType::dimension, GridType >                                      SubGridType;
  typedef Dune::SubGridList< DiscreteFunctionType, SubGridType, MacroMicroGridSpecifierType > SubGridListType;

  //! -------------------------- MsFEM error estimator ----------------------------
  typedef Dune::MsFEMErrorEstimator< DiscreteFunctionType,
                               DiffusionType,
                               FirstSourceType,
                               MacroMicroGridSpecifierType,
                               SubGridListType >
  MsFEMErrorEstimatorType;
  //! -----------------------------------------------------------------------------
#endif
  //! ------------------ typedefs and classes for data output ---------------------
  typedef Dune::tuple< const DiscreteFunctionType* >      IOTupleType;
  typedef Dune::DataOutput< GridType, IOTupleType > DataOutputType;
  // just for the discretized exact solution (in case it is available)
  typedef Dune::tuple< const DiscreteExactSolutionType* > ExSolIOTupleType;
  // just for the discretized exact solution (in case it is available)
  typedef Dune::DataOutput< GridType, ExSolIOTupleType > ExSolDataOutputType;
  
  
/// NEW RB TEST AREA
#if 1

  typedef SineReducedBasisSpace< DiscreteFunctionSpaceType, 4 > RBSpace;
  typedef AdaptiveDiscreteFunction< RBSpace > RBFunction;
  
  struct RBMatrixTraits
   {
    typedef RBSpace                                RowSpaceType;
    typedef RBSpace                                ColumnSpaceType;
    typedef Dune::LagrangeMatrixSetup< false >     StencilType;
    typedef Dune::ParallelScalarProduct< RBSpace > ParallelScalarProductType;

    template< class M >
    struct Adapter
     {
       typedef Dune::LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
     };
   };
  
   
  //typedef Dune::SparseRowMatrixObject< RBSpace, RBSpace, RBMatrixTraits > RBMatrix;
  //typedef Dune::SparseRowMatrix< double > RBMatrix;
  //Dune::Matrix< double > RBMatrix;
  
#if 0
  typedef //Dune::OEMBICGSQOp
          Dune::OEMBICGSTABOp
          //  OEMGMRESOp
          < RBFunction, RBMatrix > InverseRBMatrix;
#endif
#endif

};
} //namespace DUNE

#endif // MSFEM_TRAITS_HH
