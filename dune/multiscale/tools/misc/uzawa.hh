#ifndef DUNE_FEM_UZAWAINVERSEOPERATORS_HH
#define DUNE_FEM_UZAWAINVERSEOPERATORS_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/solver/inverseoperators.hh>

namespace Dune
{

  /** \class ConjugateGradientSolver
   *  \ingroup OEMSolver
   *  \brief   linear solver using the CG algorithm
   *
   *  \param  Operator  type of the operator to invert
   */
  template< class AType, class BType >
  struct UzawaSolver
  {
    //! type of the operators to invert
    typedef AType AOperatorType;
    typedef BType BOperatorType;

    //! field type of the operator's domain vectors
 // typedef typename AOperatorType :: DomainFieldType ADomainFieldType;
    //! field type of the operator's range vectors
    typedef double ARangeFieldType;
//      typename AOperatorType :: RangeFieldType ARangeFieldType;
    
    /*
    //! type of the operator's domain vectors
    typedef typename AOperatorType :: DomainType ADomainType;
    //! type of the operator's range vectors
    typedef typename AOperatorType :: RangeType ARangeType;


    //! field type of the operator's domain vectors
    typedef typename BOperatorType :: DomainFieldType BDomainFieldType;
    //! field type of the operator's range vectors
    typedef typename BOperatorType :: RangeFieldType BRangeFieldType;

    //! type of the operator's domain vectors
    typedef typename BOperatorType :: DomainType BDomainType;
    //! type of the operator's range vectors
    typedef typename BOperatorType :: RangeType BRangeType;
    */

  private:
    //dune_static_assert( (Conversion< ADomainType, ARangeType >::sameType), "DomainType must equal RangeType for Operator A" );
    //dune_static_assert( (Conversion< ADomainType, BDomainType >::sameType), "B has to map from A DomainType to B RangeType" );


  protected:
    const ARangeFieldType epsilon_;
    const unsigned int maxIterations_;
    const bool verbose_;
    mutable double averageCommTime_;
    mutable int realCount_;
    
  public:
    /** \brief constructor
     *
     *  \param[in]  epsilon        tolerance
     *  \param[in]  maxIterations  maximum number of CG iterations
     *  \param[in]  verbose        verbose output
     */
    UzawaSolver ( ARangeFieldType epsilon,
                  unsigned int maxIterations,
                  bool verbose )
    : epsilon_( epsilon ),
      maxIterations_( maxIterations ),
      verbose_( verbose ),
      averageCommTime_( 0.0 ),
      realCount_(0)
    {}

    /** \brief constructor
     *
     *  \param[in]  epsilon        tolerance
     *  \param[in]  maxIterations  maximum number of CG iterations
     */
    UzawaSolver ( ARangeFieldType epsilon,
                  unsigned int maxIterations )
    : epsilon_( epsilon ),
      maxIterations_( maxIterations ),
      verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
      averageCommTime_( 0.0 ),
      realCount_(0)
    {}

  private:
    // prohibit copying
    UzawaSolver ( const UzawaSolver & );

  public:
    /** \brief solve \f$op( x ) = b\f$
     *
     *  \note The CG algorithm also works for positive semidefinite operators.
     *        In this case, \f$x \cdot v = b \cdot v\f$ for all \f$v\f$ in the
     *        operator's kernel.
     *
     *  \param[in]   op  linear operator to invert (must be symmetic and
     *                   positive definite)
     *  \param[in]   b   right hand side
     *  \param       x   solution (must be initialized to a start value)
     */
    template <class ARangeType, class BRangeType, class ADomainType>
    void solve (const AOperatorType &inverseLaplace, const BOperatorType &B,  
                const ARangeType &b, const BRangeType &g, 
                ADomainType &x, BRangeType &lambda ) const
    {
      const bool verbose = (verbose_ && (b.space().grid().comm().rank() == 0));
      
      const ARangeFieldType tolerance = SQR( epsilon_ ) * b.scalarProductDofs( b ); 

      averageCommTime_ = 0.0;

      ARangeType p(b);
      B.applyTransposed(lambda , p);
      p -= b;
      p *= -1.;

      ARangeType h(b);
      h.clear();

      // apply inverse laplace
      inverseLaplace( p, x );

      BRangeType d(g);
      B(x,d);
      d -= g;


      BRangeType r(d);

      ARangeFieldType prevResiduum = 0;
      ARangeFieldType residuum = r.scalarProductDofs( r );
   
      realCount_ = 0;
      for( unsigned int count = 0; (residuum > tolerance) && (count < maxIterations_); ++count )
      {
        if( count > 0 )
        {
          d *= (residuum / prevResiduum);
          d += r;
        }

        p.clear();
        B.applyTransposed( d, p);
        // apply inverse laplace 
        inverseLaplace( p, h );

        const ARangeFieldType alpha = residuum / p.scalarProductDofs( h );
        lambda.axpy( alpha, d );
        x.axpy( (-1.*alpha), h );

        B(x,r);
        r -= g;

        prevResiduum = residuum;
        residuum = r.scalarProductDofs( r );
        
        double exchangeTime = h.space().communicator().exchangeTime();
        if( verbose )
        {
          std :: cerr << "Uzawa-Iteration: " << count << ", Residuum: " << residuum
                      << std :: endl;
          // only for parallel apps 
          if( b.space().grid().comm().size() > 1 )
            std :: cerr << "Communication needed: " << exchangeTime << " sec " << std::endl;
        }
        
        averageCommTime_ += exchangeTime;
        ++realCount_;
      }
    }

    //! number of iterations needed for last solve 
    int iterations () const 
    {
      return realCount_;
    }
    
    //! return average communication time during last solve 
    double averageCommTime() const 
    {
      return averageCommTime_;
    }
  };

  /** \class   UzawaInverseOp
   *  \ingroup OEMSolver
   *  \brief   Inversion operator using Uzawa algorithm, to solve systems in the form
   *            A B^T  u  =  x
   *            B 0    v  =  y
   *            Where A operates on u and B mappes u to v
   */
  template< class DFV, class DFP, class InvLaplace, class Div >
  struct UzawaInverseOp
  : public Operator< typename DFV::RangeFieldType, typename DFV::RangeFieldType, DFV, DFV >
  {
    typedef DFV DiscreteVelocityFunctionType;
    typedef DFP DiscretePreasureFunctionType;
    typedef InvLaplace AType;
    typedef Div BType;

    /** \brief constructor of UzawaInverseOperator
     *
     *  \param[in] Laplace A matrix
     *  \param[in] Div B matrix 
     *  \param[in] redEps reduction epsilon
     *  \param[in] absLimit absolut limit of residual
     *  \param[in] maxIter maximal iteration steps
     *  \param[in] verbose verbosity
     */
    UzawaInverseOp( const AType &A,
                    const BType &B,
                    double redEps,
                    double absLimit,
                    int maxIter,
                    bool verbose )
    : A_( A ),
      B_( B ),
      solver_( absLimit, maxIter, verbose )
    {} 

    /** \brief constructor of UzawaInverseOperator
     *
     *  \param[in] op Mapping describing operator to invert
     *  \param[in] redEps reduction epsilon
     *  \param[in] absLimit absolut limit of residual
     *  \param[in] maxIter maximal iteration steps
     */
    UzawaInverseOp( const AType &A,
                    const BType &B,
                    double redEps,
                    double absLimit,
                    int maxIter = std::numeric_limits< int >::max() )
    : A_( A ),
      B_( B ),
      solver_( absLimit, maxIter )
    {} 

    /** \brief solve the system 
        \param[in] arg right hand side 
        \param[out] dest solution 
    */
    virtual void operator() ( const DiscreteVelocityFunctionType &arg1,
                              const DiscretePreasureFunctionType &arg2,
                              DiscreteVelocityFunctionType &dest1,
                              DiscretePreasureFunctionType &dest2 ) const
    {
      solver_.solve( A_, B_, arg1, arg2, dest1, dest2 );
    }

    virtual void operator() ( const DiscreteVelocityFunctionType &arg1,
                              DiscreteVelocityFunctionType &dest1) const
    {

    }

    
    //! number of iterations needed for last solve 
    int iterations () const 
    {
      return solver_.iterations();
    }

    //! return average communication time during last solve 
    double averageCommTime() const 
    {
      return solver_.averageCommTime();
    }

  protected:
    const AType &A_;
    const BType &B_;
    const UzawaSolver< AType, BType > solver_;
  };

} // end namespace Dune

#endif
