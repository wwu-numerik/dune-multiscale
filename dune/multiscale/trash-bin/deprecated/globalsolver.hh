#ifndef DUNE_GLOBALSOLVER_HH
#define DUNE_GLOBALSOLVER_HH

//- Dune includes
#include <dune/fem/quadrature/quadrature.hh>

//- local includes
#include "globalfeop.hh"

namespace Dune 
{

// Imp = Implementation, this ending is to emphasize that we are dealing with template parameters
template< class DiscreteFunctionImp,
          class PeriodicDiscreteFunctionImp,
          class TensorImp,
          class AdvectionImp,
          class MassImp,
          class MatrixAssemblerImp >
class EllipticFEOp
  : public FEOp< DiscreteFunctionImp,
                 SparseRowMatrix< typename DiscreteFunctionImp :: RangeFieldType >,
                 EllipticFEOp< DiscreteFunctionImp,
                               PeriodicDiscreteFunctionImp,
                               TensorImp,
                               AdvectionImp,
                               MassImp,
                               MatrixAssemblerImp > >
  {
  public:

    //! Necessary typedefs for the PeriodicDiscreteFunctionImp:
    /*  The following typedefs are only important if the MatrixAssembler is a HM-FEM Assembler. In case of standard  *  FEM Assemblers, the following typedefs are simply dummies! */

    //! type of discrete functions
    typedef PeriodicDiscreteFunctionImp PeriodicDiscreteFunctionType;

    //! type of discrete function space
    typedef typename PeriodicDiscreteFunctionImp :: DiscreteFunctionSpaceType
      PeriodicDiscreteFunctionSpaceType;

    //! type of grid partition
    typedef typename PeriodicDiscreteFunctionSpaceType :: GridPartType PeriodicGridPartType;

    //! type of grid
    typedef typename PeriodicDiscreteFunctionSpaceType :: GridType PeriodicGridType;


    //! Necessary typedefs for the DiscreteFunctionImp:

    //! type of discrete functions
    typedef DiscreteFunctionImp DiscreteFunctionType;

    typedef MatrixAssemblerImp MatrixAssemblerType;

    //! type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of an element of the jacobian matrix
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    //! field type of range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;

    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef DomainFieldType TimeType;

    //! type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    //! set of base functions
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;


    //! polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

    //! The grid's dimension
    enum { dimension = GridType :: dimension };

    //! type of tensor
    typedef TensorImp TensorType;

    //! type of advection-term
    typedef AdvectionImp AdvectionType;

    //! type of mass-term 
    typedef MassImp MassTermType;

    //! type of system matrix
    typedef SparseRowMatrix< RangeFieldType > MatrixType;

    //! type of quadrature to be used
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;


  private:
    typedef EllipticFEOp< DiscreteFunctionType, PeriodicDiscreteFunctionImp, TensorType, AdvectionType, MassTermType, MatrixAssemblerImp > ThisType;
    typedef FEOp< DiscreteFunctionType, MatrixType, ThisType > BaseType;

  public:
    //! Operation mode for Finite Element Operator
    typedef typename BaseType :: OpMode OpMode;

  private:

    TensorType *stiffTensor_;
    AdvectionType *advection_;
    MassTermType *massTerm_;
    const TimeType *t_;

    const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace_;

  public:
    //! constructor - no Tensor (A(x)), no MassTerm (a(x)) and no \lambda-term
    EllipticFEOp( const DiscreteFunctionSpaceType &discreteFunctionSpace,
                  const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                  OpMode opMode)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( NULL ),
      advection_( NULL ),
      massTerm_( NULL ),
      t_( NULL ),
      periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    {
    }

    //! constructor - Tensor (A(x) but no MassTerm (a(x))
    EllipticFEOp( TensorType &stiff,
                  const DiscreteFunctionSpaceType &discreteFunctionSpace,
                  const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                  OpMode opMode)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( &stiff ),
      advection_( NULL ), 
      massTerm_( NULL ),
      t_( NULL ),
      periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    { 
    }

    //! constructor - no Tensor (A(x) but MassTerm (a(x))
    EllipticFEOp( MassTermType &mass,
                 const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                 OpMode opMode)
    : BaseType( discreteFunctionSpace, opMode),
      stiffTensor_( NULL ),
      advection_( NULL ),
      massTerm_( &mass ),
      t_( NULL ),
      periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    { 
    }

    //! constructor - both: Tensor (A(x) and MassTerm (a(x))
    EllipticFEOp( TensorType &stiff,
                  MassTermType &mass,
                  const DiscreteFunctionSpaceType &discreteFunctionSpace,
                  const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                  OpMode opMode)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( &stiff ),
      advection_( NULL ),
      massTerm_( &mass ),
      t_( NULL ),
      periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    { 
    }

    //! constructor - all: Tensor A(x), MassTerm a(x) and advection term b(x)
    EllipticFEOp( TensorType &stiff,
                  AdvectionType &advection,
                  MassTermType &mass,
                  const DiscreteFunctionSpaceType &discreteFunctionSpace,
                  const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                  OpMode opMode)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( &stiff ),
      advection_( &advection ),
      massTerm_( &mass ),
      t_( NULL ),
      periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    { 
    }

    //! constructor - all: Tensor A(x), MassTerm a(x), advection term b(x) and time t
    EllipticFEOp( TensorType &stiff,
                  AdvectionType &advection,
                  MassTermType &mass,
                  const TimeType &t,
                  const DiscreteFunctionSpaceType &discreteFunctionSpace,
                  const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                  OpMode opMode)
    : BaseType( discreteFunctionSpace, opMode ),
      stiffTensor_( &stiff ),
      advection_( &advection ),
      massTerm_( &mass ),
      t_( &t ),
      periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    { 
    }

    //! Returns the actual matrix if it is assembled
    const MatrixType* getMatrix () const
    {
      assert( this->matrix_ );
      return this->matrix_;
    }

    //! Creates a new empty matrix
    MatrixType* newEmptyMatrix () const 
    {
      return new MatrixType( this->functionSpace_.size(),
                             this->functionSpace_.size(), 
                             15 * dimension );
    }

    //! Prepares the local operator before calling apply()
    void prepareGlobal ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest )
    {
      this->arg_  = &arg;
      this->dest_ = &dest;
      this->dest_.clear();
    }


    template< class  EntityType, class LocalMatrixType > //Corresponding matrix to only one element of the grid
    void getLocalStiffMatrix( const EntityType &entity,
                              const int matrixSize,
                              LocalMatrixType &localStiffMatrix ) const
    {

     const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;

//Using a standard FEM, use:
  // MatrixAssemblerType lhsMatrix(discreteFunctionSpace, *stiffTensor_ , *massTerm_ ); 

     if ( t_ )
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_, *t_ ); 

       lhsMatrix.template getTheLocalStiffMatrix< EntityType, LocalMatrixType >
           (entity, matrixSize, localStiffMatrix );
      }
     else
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_ ); 

       lhsMatrix.template getTheLocalStiffMatrix< EntityType, LocalMatrixType >
           (entity, matrixSize, localStiffMatrix );
      }


    }


    template< class  EntityType, class LocalMatrixType >
    void getLocalMixedMatrix( const EntityType &entity,
                              const int matrixSize,
                              LocalMatrixType &localMixedMatrix ) const
    {

     const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;

     if ( t_ )
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_, *t_ ); 

       lhsMatrix.template getTheLocalMixedMatrix< EntityType, LocalMatrixType >
           (entity, matrixSize, localMixedMatrix );
      }
     else
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_ ); 

       lhsMatrix.template getTheLocalMixedMatrix< EntityType, LocalMatrixType >
           (entity, matrixSize, localMixedMatrix );
      }

    }


    template< class  EntityType, class LocalMatrixType >
    void getLocalMassMatrix( const EntityType &entity,
                             const int matrixSize,
                             LocalMatrixType &localMassMatrix) const
    {

      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;

//Using a standard FEM, use:
   // MatrixAssemblerType lhsMatrix(discreteFunctionSpace, *stiffTensor_ , *massTerm_ ); 

     if ( t_ )
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_, *t_ ); 

       lhsMatrix.template getTheLocalMassMatrix< EntityType, LocalMatrixType >
           (entity, matrixSize, localMassMatrix );
      }
     else
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_ ); 

       lhsMatrix.template getTheLocalMassMatrix< EntityType, LocalMatrixType >
           (entity, matrixSize, localMassMatrix );
      }

    }


    template< class  EntityType, class LocalMatrixType >
    void getLocalMatrix( const EntityType &entity,
                         const int matrixSize,
                         LocalMatrixType &localMatrix) const
    {
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;

      LocalMatrixType localStiffMatrix, localMassMatrix;

//Using a standard FEM, use:
   // MatrixAssemblerType lhsMatrix(discreteFunctionSpace, *stiffTensor_ , *massTerm_ ); 

     if ( t_ )
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_, *t_ );

       lhsMatrix.template getLocalLHS< EntityType, LocalMatrixType >
          (entity, matrixSize, localMatrix );
      }
     else
      {
       MatrixAssemblerType lhsMatrix
           (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_ ); 

       lhsMatrix.template getLocalLHS< EntityType, LocalMatrixType >
          (entity, matrixSize, localMatrix );
      }

    }


  }; // end class


} // end namespace 

#endif
