#ifndef DUNE_GLOBALSOLVER_HH
#define DUNE_GLOBALSOLVER_HH

//- Dune includes
#include <dune/fem/quadrature/quadrature.hh>

//- local includes
#include "globalfeop.hh"
#include "cellproblemsolver.hh"

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
    typedef typename PeriodicDiscreteFunctionImp :: FunctionSpaceType
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
    typedef typename DiscreteFunctionType :: FunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of an element of the jacobian matrix
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    //! field type of range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;

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

     MatrixAssemblerType lhsMatrix
         (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_ ); 

     lhsMatrix.template getTheLocalStiffMatrix< EntityType, LocalMatrixType >
         (entity, matrixSize, localStiffMatrix );

    }


    template< class  EntityType, class LocalMatrixType >
    void getLocalMixedMatrix( const EntityType &entity,
                              const int matrixSize,
                              LocalMatrixType &localMixedMatrix ) const
    {

     const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;

     MatrixAssemblerType lhsMatrix
         (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_ ); 

     lhsMatrix.template getTheLocalMixedMatrix< EntityType, LocalMatrixType >
         (entity, matrixSize, localMixedMatrix );

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

      MatrixAssemblerType lhsMatrix
          (discreteFunctionSpace, periodicDiscreteFunctionSpace_, *stiffTensor_ , *advection_, *massTerm_ ); 

      lhsMatrix.template getTheLocalMassMatrix< EntityType, LocalMatrixType >
          (entity, matrixSize, localMassMatrix );

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

      MatrixAssemblerType 
         lhsMatrix(discreteFunctionSpace, periodicDiscreteFunctionSpace_,*stiffTensor_ , *advection_ , *massTerm_ ); 

      lhsMatrix.template getLocalLHS< EntityType, LocalMatrixType >
         (entity, matrixSize, localMatrix );

    }


  }; // end class

  // Assembler for right rand side
  template< class DiscreteFunctionImp , class TensorImp, class SecondSourceImp >
  class RightHandSideAssembler
  {  
  public:

    typedef DiscreteFunctionImp DiscreteFunctionType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionType :: LocalFunctionType
      LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

    typedef typename GridPartType :: GridType GridType;

    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    enum { dimension = GridType :: dimension };

    typedef TensorImp TensorType;

    typedef SecondSourceImp SecondSourceType;
 
  private:

    const TensorType *tensor_;
    const SecondSourceType *G_;

  public:

    //! constructor - Tensor (A(x)), MassTerm (a(x)) and (-\div G)-Term
    template <class IrgendwasImp >
    RightHandSideAssembler(IrgendwasImp &h)
    : tensor_( NULL ),
      G_( NULL )
    {
    }

    //! constructor - Tensor (A(x)), MassTerm (a(x)) and (-\div G)-Term
    template <class TensorType, class SecondSourceType >
    RightHandSideAssembler( const TensorType &tensor,
                            const SecondSourceType &G)
    : tensor_( &tensor ),
      G_( &G )
    {
    }

    //! constructor - Tensor (A(x)), MassTerm (a(x)) and (-\div G)-Term
    RightHandSideAssembler()
    : tensor_( NULL ),
      G_( NULL )
    {
    }

  public:

    // discreteFunction is an output parameter (kind of return value)
    template< int polOrd, class FunctionType >
    void assemble( const FunctionType &function,
                         DiscreteFunctionType &discreteFunction)
    //discreteFunction ist der Rueckgabewert der funktion 'assemble'. Hierin wird sozusagen die rechte Seite gespeichert

    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace    
        = discreteFunction.space();

      // set discreteFunction to zero:
      discreteFunction.clear();

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {
        //it* Pointer auf ein Element der Entity
        const GeometryType &geometry = (*it).geometry(); //Referenz auf Geometrie

        LocalFunctionType elementOfRHS = discreteFunction.localFunction( *it ); //*it zeigt auf ein bestimmtes Element der entity
	//hier wird sozusagen ein Pointer von localFunction auf discreteFunction erzeugt. Befinden wir uns auf einer bestimmten entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in discreteFunction(aktuelleEntity)	

        const BaseFunctionSetType baseSet //BaseFunctions leben immer auf Refernzelement!!!
          = discreteFunctionSpace.baseFunctionSet( *it ); //*it Referenz auf eine bestimmtes Element der entity. In der ersten Klasse war das Element fest, deshalb konnte man sich dort Pointer sparen. //loeschen: discreteFunctionSpace statt functionSpace

        CachingQuadrature< GridPartType, 0 > quadrature( *it, polOrd ); //0 --> codim 0

        const int numDofs = elementOfRHS.numDofs(); //Dofs = Freiheitsgrade (also die Unbekannten)
	for( int i = 0; i < numDofs; ++i ) //Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
        { 
          //std :: cout << i << "discreteFunction.localFunction( *it ) " << i << ": " << discreteFunction.localFunction( *it )[ i ] << std :: endl;

          // the return values:
          RangeType y, z;

          RangeType a[dimension][dimension];

          JacobianRangeType gradientPhi;

          // to save: A \nabla PHI_H
          RangeType w[dimension];

          // to save: A \nabla PHI_H * \nabla phi_h;
          RangeType t = 0; 

          const int numQuadraturePoints = quadrature.nop();
          for( int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
          {
            const double det
              = geometry.integrationElement( quadrature.point( quadraturePoint ) );

	    // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'y':
	    function.evaluate( geometry.global( quadrature.point( quadraturePoint ) ) , y );

	    // evaluate the current base function at the current quadrature point and save its value in 'z':
	    baseSet.evaluate( i, quadrature, quadraturePoint, z ); //i = i'te Basisfunktion;

            if ( tensor_ )
            {
	     // evaluate the gradient of the current base function at the current quadrature point and save its value in 'returnGradient':
	     baseSet.jacobian( i, quadrature, quadraturePoint, gradientPhi );
             // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:	    

             const FieldMatrix< double, dimension, dimension > &inv
                 = geometry.jacobianInverseTransposed( quadrature.point( quadraturePoint ) );

             // multiply with transpose of jacobian inverse 
             gradientPhi[ 0 ] = FMatrixHelp :: mult( inv, gradientPhi[ 0 ] );

             // set all entries of w to zero to delete old data of a former loop cycle
             for( int k = 0; k < dimension; ++k ) 
               {  
                 w[ k ] = 0;
               }

             // the same for t:
             t = 0;

             // evaluate the tensor at the current quadrature point and save its value in 'w':
              for( int k = 0; k < dimension; ++k ) 
                {
                 for( int l = 0; l < dimension; ++l )
                    {
                      tensor_->evaluate( k , l , geometry.global( quadrature.point( quadraturePoint ) ) , a[k][l] );
                    }
                 }

            // to illustrate what we are doing, see the following (dimension-restricted) alternatives:

            #if 0
             // 1-D:
             tensor.evaluate( 0 , 0 , geometry.global( quadrature.point( quadraturePoint ) ) , a[0][0] );
 
             // 2-D:
             tensor.evaluate( 0 , 0 , geometry.global( quadrature.point( quadraturePoint ) ) , a[0][0] );
             tensor.evaluate( 0 , 1 , geometry.global( quadrature.point( quadraturePoint ) ) , a[0][1] );
             tensor.evaluate( 1 , 0 , geometry.global( quadrature.point( quadraturePoint ) ) , a[1][0] );
             tensor.evaluate( 1 , 1 , geometry.global( quadrature.point( quadraturePoint ) ) , a[1][1] );
 
             // independent of dimension:
             for( int k = 0; k < dimension; ++k ) 
               {  
                 for( int l = 0; l < dimension; ++l )
                    {
                      tensor.evaluate( k , l , geometry.global( quadrature.point( quadraturePoint ) ) , a[k][l] );
                    }
               }
            #endif


            // evaluate the gradient of Phi_H at the current quadrature point and save its value in 'v':
            for( int k = 0; k < dimension; ++k ) 
              {  
                G_->evaluate( k , geometry.global( quadrature.point( quadraturePoint ) ) , w[k] );
              }

            // illustrating alternatives:
            #if 0
              // 1-D:
              gradPhiH.evaluate( 0 , geometry.global( quadrature.point( quadraturePoint ) ) , v[0] );
 
              // 2-D:
              gradPhiH.evaluate( 0 , geometry.global( quadrature.point( quadraturePoint ) ) , v[0] );
              gradPhiH.evaluate( 1 , geometry.global( quadrature.point( quadraturePoint ) ) , v[1] );
 
              // independent of dimension:
              for( int k = 0; k < dimension; ++k ) 
               {  
                gradPhiH.evaluate( k , geometry.global( quadrature.point( quadraturePoint ) ) , v[k] );
               }
            #endif

            for( int k = 0; k < dimension; ++k ) 
                  {  
                    t += w[ k ] * gradientPhi[ 0 ][ k ];
                  }

            // illustrating alternatives:
            #if 0
             // 1-D:
             // t = w[ 0 ] * gradientPhi[ 0 ][ 0 ];
 
             // 2-D:
             // t = w[ 0 ] * gradientPhi[ 0 ][ 0 ] + w[ 1 ] * gradientPhi[ 0 ][ 1 ];
             // t = w[ 0 ] + w[ 1 ];
 
             // independent of dimension:
             for( int k = 0; k < dimension; ++k ) 
                  {  
                    t += w[ k ] * gradientPhi[ 0 ][ k ];
                  }
             #endif
            } //end of if-loop

           elementOfRHS[ i ] += det * quadrature.weight( quadraturePoint ) * (y * z);

           if (tensor_)
           {
      	     elementOfRHS[ i ] += det * quadrature.weight( quadraturePoint ) * (t);
           } // end of if-loop  

          }

        }
      }
    }  // end method


    static void printRHS(const DiscreteFunctionType &rhs)    
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace    
        = rhs.space();

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {

        LocalFunctionType elementOfRHS = rhs.localFunction( *it ); 

        const int numDofs = elementOfRHS.numDofs(); 
	for( int i = 0; i < numDofs; ++i )
        { 
         std :: cout << "Number of Dof: " << i << " ; " << rhs.name() << " : " << elementOfRHS[ i ] << std :: endl;
        }

      }

    }  // end method


//! build the right hand side for a discrete parabolic problem that is solved with backward Euler:
// Note that in comparison to the elliptic case:
//   1. f needs to be muliplied with the time step size and
//   2. the solution u_H^{(k)} of the preceeding time step needs to be added to the right hand side

   // Wir haben die rechte Seite = f, wollen darauf aber noch eine discrete Function aufaddieren, die wir discFuncToAdd nennen.
   // am Ende soll also nicht nur die rechte Seite durch f gebildet werden, sondern durch f+discFuncToAdd
   // der Wert wird in discreteFunction gespeichert
   // add discreteFunction is an output parameter (kind of return value)
    template< int polOrd, class FunctionType >
    void assembleParabolic( const FunctionType &function,
                            const RangeType &time_step_size,
                            const DiscreteFunctionType &u_H_k,
                            DiscreteFunctionType &rhs)
    // rhs is the vector that occurs on the right hand side of the linear system of equations that is to solve

    //rhs ist der Rueckgabewert der funktion 'assemble'. Hierin wird sozusagen die rechte Seite gespeichert

    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace    
        = rhs.space();

      // set rhs to zero:
      rhs.clear();

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {
        //it* Pointer auf ein Element der Entity
        const GeometryType &geometry = (*it).geometry(); //Referenz auf Geometrie

        LocalFunctionType elementOfRHS = rhs.localFunction( *it ); //*it zeigt auf ein bestimmtes Element der entity
	//hier wird sozusagen ein Pointer von localFunction auf rhs erzeugt. Befinden wir uns auf einer bestimmten entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in rhs(aktuelleEntity)	
        LocalFunctionType elementOf_u_H_k = u_H_k.localFunction( *it );

        const BaseFunctionSetType baseSet //BaseFunctions leben immer auf Refernzelement!!!
          = discreteFunctionSpace.baseFunctionSet( *it ); //*it Referenz auf eine bestimmtes Element der entity. In der ersten Klasse war das Element fest, deshalb konnte man sich dort Pointer sparen. //loeschen: discreteFunctionSpace statt functionSpace

        CachingQuadrature< GridPartType, 0 > quadrature( *it, polOrd ); //0 --> codim 0

        const int numDofs = elementOfRHS.numDofs(); //Dofs = Freiheitsgrade (also die Unbekannten)
	for( int i = 0; i < numDofs; ++i ) //Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
        { 

          // the return values:
          RangeType y, z;

          RangeType a[dimension][dimension];

          JacobianRangeType gradientPhi;

          // to save: A \nabla PHI_H
          RangeType w[dimension];

          // to save: A \nabla PHI_H * \nabla phi_h;
          RangeType t = 0; 

          const int numQuadraturePoints = quadrature.nop();
          for( int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
          {
            const double det
              = geometry.integrationElement( quadrature.point( quadraturePoint ) );

	    // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'y':
	    function.evaluate( geometry.global( quadrature.point( quadraturePoint ) ) , y );

	    // evaluate the current base function at the current quadrature point and save its value in 'z':
	    baseSet.evaluate( i, quadrature, quadraturePoint, z ); //i = i'te Basisfunktion;

            if ( tensor_ )
            {
	     // evaluate the gradient of the current base function at the current quadrature point and save its value in 'returnGradient':
	     baseSet.jacobian( i, quadrature, quadraturePoint, gradientPhi );
             // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:	    

             const FieldMatrix< double, dimension, dimension > &inv
                 = geometry.jacobianInverseTransposed( quadrature.point( quadraturePoint ) );

             // multiply with transpose of jacobian inverse 
             gradientPhi[ 0 ] = FMatrixHelp :: mult( inv, gradientPhi[ 0 ] );

             // set all entries of w to zero to delete old data of a former loop cycle
             for( int k = 0; k < dimension; ++k ) 
               {  
                 w[ k ] = 0;
               }

             // the same for t:
             t = 0;

             // evaluate the tensor at the current quadrature point and save its value in 'w':
              for( int k = 0; k < dimension; ++k ) 
                {
                 for( int l = 0; l < dimension; ++l )
                    {
                      tensor_->evaluate( k , l , geometry.global( quadrature.point( quadraturePoint ) ) , a[k][l] );
                    }
                 }

            // to illustrate what we are doing, see the following (dimension-restricted) alternatives:

            #if 0
             // 1-D:
             tensor.evaluate( 0 , 0 , geometry.global( quadrature.point( quadraturePoint ) ) , a[0][0] );
 
             // 2-D:
             tensor.evaluate( 0 , 0 , geometry.global( quadrature.point( quadraturePoint ) ) , a[0][0] );
             tensor.evaluate( 0 , 1 , geometry.global( quadrature.point( quadraturePoint ) ) , a[0][1] );
             tensor.evaluate( 1 , 0 , geometry.global( quadrature.point( quadraturePoint ) ) , a[1][0] );
             tensor.evaluate( 1 , 1 , geometry.global( quadrature.point( quadraturePoint ) ) , a[1][1] );
 
             // independent of dimension:
             for( int k = 0; k < dimension; ++k ) 
               {  
                 for( int l = 0; l < dimension; ++l )
                    {
                      tensor.evaluate( k , l , geometry.global( quadrature.point( quadraturePoint ) ) , a[k][l] );
                    }
               }
            #endif


            // evaluate the gradient of Phi_H at the current quadrature point and save its value in 'v':
            for( int k = 0; k < dimension; ++k ) 
              {  
                G_->evaluate( k , geometry.global( quadrature.point( quadraturePoint ) ) , w[k] );
              }

            // illustrating alternatives:
            #if 0
              // 1-D:
              gradPhiH.evaluate( 0 , geometry.global( quadrature.point( quadraturePoint ) ) , v[0] );
 
              // 2-D:
              gradPhiH.evaluate( 0 , geometry.global( quadrature.point( quadraturePoint ) ) , v[0] );
              gradPhiH.evaluate( 1 , geometry.global( quadrature.point( quadraturePoint ) ) , v[1] );
 
              // independent of dimension:
              for( int k = 0; k < dimension; ++k ) 
               {  
                gradPhiH.evaluate( k , geometry.global( quadrature.point( quadraturePoint ) ) , v[k] );
               }
            #endif

            for( int k = 0; k < dimension; ++k ) 
                  {  
                    t += w[ k ] * gradientPhi[ 0 ][ k ];
                  }

            // illustrating alternatives:
            #if 0
             // 1-D:
             // t = w[ 0 ] * gradientPhi[ 0 ][ 0 ];
 
             // 2-D:
             // t = w[ 0 ] * gradientPhi[ 0 ][ 0 ] + w[ 1 ] * gradientPhi[ 0 ][ 1 ];
             // t = w[ 0 ] + w[ 1 ];
 
             // independent of dimension:
             for( int k = 0; k < dimension; ++k ) 
                  {  
                    t += w[ k ] * gradientPhi[ 0 ][ k ];
                  }
             #endif
            } //end of if-loop

           // value of u_H_k:
           RangeType val_to_add;
           elementOf_u_H_k.evaluate( quadrature , quadraturePoint , val_to_add );

           elementOfRHS[ i ] += time_step_size * det * quadrature.weight( quadraturePoint ) * (y * z);


           //add the local value of discFuncToAdd to the right hand side
           elementOfRHS[ i ] += det * quadrature.weight( quadraturePoint ) * val_to_add * z;

           if (tensor_)
           {
      	     elementOfRHS[ i ] += det * quadrature.weight( quadraturePoint ) * (t);
           } // end of if-loop  

          }

        }
      }
    }  // end method


  };

} // end namespace 

#endif
