#ifndef DUNE_RB_REDUCEDBASISSPACE_REDUCEDBASISSPACE_INLINE_HH
#define DUNE_RB_REDUCEDBASISSPACE_REDUCEDBASISSPACE_INLINE_HH

#include "reducedbasisspace.hh"

namespace Dune
{

namespace Multiscale
{

  // ReducedBasisSpace
  // -----------------

  template< class BaseFunction >
  inline ReducedBasisSpace< BaseFunction >
  :: ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace,
                         const InterfaceType commInterface,
                         const CommunicationDirection commDirection)
    : BaseType( baseFunctionSpace.gridPart(), commInterface, commDirection  ),
      baseFunctionSpace_( baseFunctionSpace ),
      baseFunctionList_( baseFunctionSpace_, "ReducedSpace_BaseFunctions" ),
      //    baseFunctionList_(),
      mapper_( baseFunctionList_ ),
      pipeline_(0)
  {}

  // ReducedBasisSpace
  // -----------------

  template< class BaseFunction >
  inline ReducedBasisSpace< BaseFunction >
  :: ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace,
                         const std::string &listName,
                         const InterfaceType commInterface,
                         const CommunicationDirection commDirection)
    : BaseType( baseFunctionSpace.gridPart(), commInterface, commDirection  ),
      baseFunctionSpace_( baseFunctionSpace ),
      baseFunctionList_( baseFunctionSpace_, listName, listName + ".lst" ),
      //    baseFunctionList_(),
      mapper_( baseFunctionList_ ),
      pipeline_(0)
  {}

  template< class BaseFunction >
  template< class StreamTraits >
  inline ReducedBasisSpace< BaseFunction >
  :: ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace,
                         InStreamInterface< StreamTraits > &in,
                         const InterfaceType commInterface,
                         const CommunicationDirection commDirection)
    : BaseType( baseFunctionSpace.gridPart(), commInterface, commDirection  ),
      baseFunctionSpace_( baseFunctionSpace ),
      //    baseFunctionList_(),
      baseFunctionList_( baseFunctionSpace_, "ReducedSpace_BaseFunctions" ),
      mapper_( baseFunctionList_ ),
      pipeline_(0)
      //   {
      //     read( in );
      //   }
  {}
  template< class BaseFunction >
  inline ReducedBasisSpace< BaseFunction > :: ~ReducedBasisSpace ()
  {
    clear();
  }


  template< class BaseFunction >
  inline void ReducedBasisSpace< BaseFunction >
  :: addBaseFunction ( const BaseFunctionType &baseFunction )
  {
    baseFunctionList_.push_back( baseFunction );
  }


  template< class BaseFunction >
  inline void
  ReducedBasisSpace< BaseFunction > :: baseFunction ( unsigned int i, BaseFunctionType &ret ) const
  {
    baseFunctionList_.getFunc(i, ret);
  }

  template< class BaseFunction >
  inline void
  ReducedBasisSpace< BaseFunction > :: setBaseFunction (unsigned int i, BaseFunctionType &arg)
  {
    baseFunctionList_.setFunc(i, arg);
  }


  template< class BaseFunction >
  inline const typename ReducedBasisSpace< BaseFunction > :: BaseFunctionSpaceType &
  ReducedBasisSpace< BaseFunction > :: baseFunctionSpace () const
  {
    return baseFunctionSpace_;
  }


  template< class BaseFunction >
  inline void ReducedBasisSpace< BaseFunction > :: clear ()
  {
    //! \todo (s_kaul01#7#): implement me!!!
    //    crop( 0 );
    //std :: cout << "Implement me!!" << std :: endl;
    runPipeline();
  }

  //crop would mean we know where which function is stored, but we don't in general
  //so I commented it out

  //   template< class BaseFunction >
  //   inline void ReducedBasisSpace< BaseFunction > :: crop ( unsigned int n )
  //   {
  //     const unsigned int size = baseFunctionList_.size();
  //     assert( n <= size );
  //     for( unsigned int i = n; i < size; ++i )
  //       delete baseFunctionList_[ i ];
  //     baseFunctionList_.resize( n );
  //   }


  template< class BaseFunction >
  inline unsigned int ReducedBasisSpace< BaseFunction > :: numBaseFunctions () const
  {
    return baseFunctionList_.size();
  }


  template< class BaseFunction >
  template< class DiscreteFunctionType >
  inline void ReducedBasisSpace< BaseFunction >
  :: project ( const DiscreteFunctionType &sourceFunction,
               BaseFunctionType &destFunction ) const
  {
    typedef typename DiscreteFunctionType :: RangeFieldType DofType;

    const unsigned int size = baseFunctionList_.size();

    destFunction.clear();
    BaseFunctionType baseFunction("baseFunction", baseFunctionSpace_);
    for ( unsigned int i = 0; i < size; ++i )
      {
        baseFunction.clear();
        baseFunctionList_.getFunc(i, baseFunction);
        const DofType &dof = sourceFunction.dof( i );
        destFunction.addScaled( baseFunction, dof );
      }
  }


//   template< class BaseFunction >
//   template< class DiscreteFunctionType >
//   inline void ReducedBasisSpace< BaseFunction >
//   :: restrictFunction ( const BaseFunctionType &sourceFunction,
//                      DiscreteFunctionType &destFunction ) const
//   {
//     const unsigned int size = baseFunctionList_.size();
//     assert( size == destFunction.size() );

//     BaseFunctionType baseFunction("baseFunction", baseFunctionSpace_);
//     for ( unsigned int i = 0; i < size; ++i )
//       {
//      baseFunction.clear();
//      baseFunctionList_.getFunc(i, baseFunction);
//      destFunction.dof( i )
//        = baseFunction.scalarProductDofs( sourceFunction );
//       }
//   }

  template< class BaseFunction >
  template< class SourceFunctionType, class DiscreteFunctionType >
  inline void ReducedBasisSpace< BaseFunction >
  :: restrictFunction ( const SourceFunctionType &sourceFunction,
                        DiscreteFunctionType &destFunction ) const
  {
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

    typedef typename BaseFunctionType::LocalFunctionType BaseLocalFunctionType;
    typedef typename SourceFunctionType::LocalFunctionType SourceLocalFunctionType;
    typedef typename Traits :: DiscreteFunctionSpaceType                              DiscreteFunctionSpaceType;
    typedef typename Traits :: IteratorType                        IteratorType;
    typedef typename Traits :: GridPartType                        GridPartType;
    typedef typename GridType :: template Codim< 0 > :: Entity     ElementType;
    typedef typename Traits :: RangeType                           RangeType;

#ifndef NDEBUG
    const unsigned int size = baseFunctionList_.size();
    assert( size == (unsigned)destFunction.size() );
#endif

    GridPartType &gridPart = baseFunctionSpace_.gridPart();

    BaseFunctionType baseFunction("baseFunction", baseFunctionSpace_);

    //integrate as mentioned in the long description

    //loop over grid
    IteratorType gridEnd = gridPart.template end<0>();
    for (IteratorType gridIt = gridPart.template begin<0>();
         gridIt != gridEnd; ++gridIt) {
      int i=0;
      DofIteratorType funcEnd = destFunction.dend();
      for ( DofIteratorType funcIt = destFunction.dbegin();
            funcIt != funcEnd; ++funcIt, ++i )
        {
          baseFunction.clear();
          baseFunctionList_.getFunc(i, baseFunction);
//        baseFunction.print(std::cout);
          double entityIntegral=0.0;
          //get a quadrature for this entity
          int order_ = 0; //order for the quadrature
          CachingQuadrature<GridPartType, 0> quadrature(*gridIt, order_);
          //loop over all quadrature points

          for (unsigned int quadPoint=0; quadPoint != quadrature.nop(); ++quadPoint) {
            double integrationElement = gridIt->geometry().integrationElement(quadrature.point(quadPoint));
            double weight = quadrature.weight(quadPoint);

            // evaluate the current base function
            BaseLocalFunctionType baseLocalFunction = baseFunction.localFunction(*gridIt);
            RangeType baseFunctionValue;
            baseLocalFunction.evaluate(quadrature.point(quadPoint), baseFunctionValue);

            // evaluate the source function
            SourceLocalFunctionType sourceLocalFunction = sourceFunction.localFunction(*gridIt);
            RangeType sourceFunctionValue;
            sourceLocalFunction.evaluate( quadrature.point(quadPoint), sourceFunctionValue );
            entityIntegral += integrationElement * weight * baseFunctionValue * sourceFunctionValue;

          } //loop over all quadrature points
          //accumulate integral for this base function
          *funcIt += entityIntegral;
        } // loop over all base functions
    } //loop over grid
  }

  template< class BaseFunction >
  template< class MatrixOfflineType, class MatrixOnlineType >
  inline void ReducedBasisSpace< BaseFunction >
  :: restrictMatrix ( const MatrixOfflineType &matrixOffline,
                      MatrixOnlineType &matrixOnline ) const
  {
    const unsigned int size = baseFunctionList_.size();
    assert( (size > 0) && (matrixOnline.cols() == size)
            && (matrixOnline.rows() == size) );
    assert( (matrixOffline.cols() == baseFunctionSpace_.size())
            && (matrixOffline.rows() == baseFunctionSpace_.size()) );

    assert(false);

    BaseFunctionType Sphi("Sphi", baseFunctionSpace_);
    baseFunctionList_.getFunc(0, Sphi); //(!)
    BaseFunctionType phi("phi", baseFunctionSpace_);
    BaseFunctionType psi("psi", baseFunctionSpace_);

    for ( unsigned int i = 0; i < size; ++i )
      {
        phi.clear();
        baseFunctionList_.getFunc(i, phi);
        matrixOffline( phi, Sphi );

        for ( unsigned int j = 0; j < size; ++j)
          {
            psi.clear();
            baseFunctionList_.getFunc(j, psi);
            matrixOnline.set( i, j, Sphi.scalarProductDofs( psi ) );
          }
      }
  }

  template< class BaseFunction >
  template< class OperatorComponentType, class MatrixOnlineType >
  inline void ReducedBasisSpace< BaseFunction >
  :: restrictOperator ( const OperatorComponentType &opOffline,
                        MatrixOnlineType &matrixOnline ) const
  {
    typedef typename BaseFunction :: RangeFieldType                  FieldType;
    const unsigned int size = baseFunctionList_.size();
    assert( (size > 0) && (matrixOnline.cols() == size)
            && (matrixOnline.rows() == size) );
    //         assert( (matrixOffline.cols() == baseFunctionSpace_.size())
    //                 && (matrixOffline.rows() == baseFunctionSpace_.size()) );

    BaseFunctionType op_of_phi("op_of_phi", baseFunctionSpace_);
    baseFunctionList_.getFunc(0, op_of_phi); //(!)
    BaseFunctionType phi("phi", baseFunctionSpace_);
    BaseFunctionType psi("psi", baseFunctionSpace_);

    for ( unsigned int i = 0; i < size; ++i )
      {
        phi.clear();
        baseFunctionList_.getFunc(i, phi);
        //! \todo we can evaluate the operator opOffline inside the l2ScalarProduct locally only (saves one grid iteration)
        opOffline( phi, op_of_phi );

        for ( unsigned int j = 0; j < size; ++j)
          {
            psi.clear();
            baseFunctionList_.getFunc(j, psi);
            FieldType scalar_product;
            l2ScalarProduct(op_of_phi, psi, scalar_product);
            matrixOnline.set( i, j, scalar_product );
          }
      }
  }

  template< class BaseFunction >
  template< class OperatorComponentType,
            class FuncComponentType,
            class DiscreteFunctionType >
  inline void ReducedBasisSpace< BaseFunction >
  :: restrictOperatorFunc ( const OperatorComponentType &opOffline,
                            const FuncComponentType &func,
                        DiscreteFunctionType &vectorOnline ) const
  {
    typedef typename FuncComponentType :: LocalFunctionType FuncLocalFunctionType;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

    typedef typename BaseFunctionType::LocalFunctionType BaseLocalFunctionType;

    typedef typename Traits :: DiscreteFunctionSpaceType                              DiscreteFunctionSpaceType;
    typedef typename Traits :: IteratorType                        IteratorType;
    typedef typename Traits :: GridPartType                        GridPartType;
    typedef typename GridType :: template Codim< 0 > :: Entity     ElementType;
    typedef typename Traits :: RangeType                           RangeType;

    const unsigned int size = baseFunctionList_.size();
    assert( (size > 0) && (vectorOnline.size() == size));

    BaseFunctionType Sphi("Sphi", baseFunctionSpace_);
    baseFunctionList_.getFunc(0, Sphi); //(!)

    BaseFunctionType baseFunction("baseFunction", baseFunctionSpace_);
    GridPartType &gridPart = baseFunctionSpace_.gridPart();

        // Integrate Sphi * func over Omega
        //loop over grid
        IteratorType gridEnd = gridPart.template end<0>();
        for (IteratorType gridIt = gridPart.template begin<0>();
             gridIt != gridEnd; ++gridIt) {
          int i=0;
          DofIteratorType funcEnd = vectorOnline.dend();
          for ( DofIteratorType funcIt = vectorOnline.dbegin();
                funcIt != funcEnd; ++funcIt, ++i )
            {
              //              ElementType &entity = *gridIt;
              baseFunction.clear();
              baseFunctionList_.getFunc(i, baseFunction);
              opOffline( baseFunction, Sphi );
              //          baseFunction.print(std::cout);
              double entityIntegral=0.0;
              //get a quadrature for this entity
              int order_ = 2; //order for the quadrature
              CachingQuadrature<GridPartType, 0> quadrature(*gridIt, order_);
              //loop over all quadrature points

              for (int quadPoint=0; quadPoint != quadrature.nop(); ++quadPoint) {
                double integrationElement = gridIt->geometry().integrationElement(quadrature.point(quadPoint));
                double weight = quadrature.weight(quadPoint);

                // evaluate the current base function
                BaseLocalFunctionType SphiLocalFunction = Sphi.localFunction(*gridIt);
                RangeType SphiFunctionValue;
                SphiLocalFunction.evaluate(quadrature.point(quadPoint), SphiFunctionValue);

                // evaluate the source function
                FuncLocalFunctionType funcLocal = func.localFunction(*gridIt);
                RangeType funcFunctionValue;
                funcLocal.evaluate( quadrature.point(quadPoint), funcFunctionValue );
                entityIntegral += integrationElement * weight * SphiFunctionValue * funcFunctionValue;
              } //loop over all quadrature points
              //accumulate integral for this base function
              *funcIt += entityIntegral;
            } // loop over all base functions
        } //loop over grid

  }

 template< class BaseFunction >
 template< class FuncComponentType1, class FuncComponentType2 >
  inline void ReducedBasisSpace< BaseFunction >
  :: l2ScalarProduct ( const FuncComponentType1 &func1,
                       const FuncComponentType2 &func2,
                       double &value ) const
  {
    typedef typename FuncComponentType1 :: LocalFunctionType FuncLocalFunctionType1;
    typedef typename FuncComponentType2 :: LocalFunctionType FuncLocalFunctionType2;

    typedef typename BaseFunctionType::LocalFunctionType BaseLocalFunctionType;

    typedef typename Traits :: DiscreteFunctionSpaceType                              DiscreteFunctionSpaceType;
    typedef typename Traits :: IteratorType                        IteratorType;
    typedef typename Traits :: GridPartType                        GridPartType;
    typedef typename GridType :: template Codim< 0 > :: Entity     ElementType;
    typedef typename Traits :: RangeType                           RangeType;

    GridPartType &gridPart = baseFunctionSpace_.gridPart();

    value = 0.0;

        // Integrate func * func over Omega
        //loop over grid
        IteratorType gridEnd = gridPart.template end<0>();
        for (IteratorType gridIt = gridPart.template begin<0>();
             gridIt != gridEnd; ++gridIt) {
          //          ElementType &entity = *gridIt;
              double entityIntegral=0.0;
              //get a quadrature for this entity
              int order_ = 2; //order for the quadrature
              CachingQuadrature<GridPartType, 0> quadrature(*gridIt, order_);
              //loop over all quadrature points

              for (int quadPoint=0; quadPoint != quadrature.nop(); ++quadPoint) {
                double integrationElement = gridIt->geometry().integrationElement(quadrature.point(quadPoint));
                double weight = quadrature.weight(quadPoint);

                // evaluate the source function
                FuncLocalFunctionType1 funcLocal1 = func1.localFunction(*gridIt);
                FuncLocalFunctionType2 funcLocal2 = func2.localFunction(*gridIt);
                RangeType funcFunctionValue1;
                RangeType funcFunctionValue2;
                funcLocal1.evaluate( quadrature.point(quadPoint), funcFunctionValue1 );
                funcLocal2.evaluate( quadrature.point(quadPoint), funcFunctionValue2 );
                double SfuncValue = funcFunctionValue1 * funcFunctionValue2;
                entityIntegral += integrationElement * weight * SfuncValue;
              } //loop over all quadrature points
              value += entityIntegral;
        } //loop over grid
  }


  template< class BaseFunction >
  template< class OperatorComponentType1,
            class OperatorComponentType2,
            class MatrixOnlineType >
  inline void ReducedBasisSpace< BaseFunction >
  :: restrictTwoOperators ( const OperatorComponentType1 &opOffline1,
                            const OperatorComponentType2 &opOffline2,
                        MatrixOnlineType &matrixOnline ) const
  {
    typedef typename BaseFunction :: RangeFieldType                  FieldType;
    const unsigned int size = baseFunctionList_.size();
    assert( (size > 0) && (matrixOnline.cols() == size)
            && (matrixOnline.rows() == size) );
    //         assert( (matrixOffline.cols() == baseFunctionSpace_.size())
    //                 && (matrixOffline.rows() == baseFunctionSpace_.size()) );

    BaseFunctionType op1_of_phi("op1_of_phi", baseFunctionSpace_);
    BaseFunctionType op2_of_psi("op2_of_psi", baseFunctionSpace_);
/*    baseFunctionList_.getFunc(0, op1_of_phi); //(!)
 *    baseFunctionList_.getFunc(0, op2_of_psi); //(!)*/
    BaseFunctionType phi("phi", baseFunctionSpace_);
    BaseFunctionType psi("psi", baseFunctionSpace_);

    for ( unsigned int i = 0; i < size; ++i )
      {
        phi.clear();
        baseFunctionList_.getFunc(i, phi);
        //! \todo we can evaluate the operator opOffline1 inside the l2ScalarProduct locally only (saves one grid iteration)
        opOffline1( phi, op1_of_phi );

        for ( unsigned int j = 0; j < size; ++j)
          {
            psi.clear();
            baseFunctionList_.getFunc(j, psi);
            //! \todo we can evaluate the operator opOffline2 inside the l2ScalarProduct locally only (saves one grid iteration)
            opOffline2( psi, op2_of_psi );

            FieldType scalar_product;
            l2ScalarProduct( op1_of_phi, op2_of_psi, scalar_product );

            matrixOnline.set( i, j, scalar_product );
          }
      }
  }


  //read doesn't quite match into  our current setting, do we need it?

  //   template< class BaseFunction >
  //   template< class StreamTraits >
  //   inline void ReducedBasisSpace< BaseFunction >
  //     :: read ( InStreamInterface< StreamTraits > &in )
  //   {
  //     clear();

  //     unsigned int size;
  //     in >> size;
  //     baseFunctionList_.resize( size );

  //     for( unsigned int i = 0; i < size; ++i )
  //     {
  //       BaseFunctionType *baseFunction
  //         = new BaseFunctionType( "BaseFunction", baseFunctionSpace() );
  //       in >> *baseFunction;
  //       baseFunctionList_.getFuncById( i ) = baseFunction;
  //     }
  //   }


  //same as for read

  //   template< class BaseFunction >
  //   template< class StreamTraits >
  //   inline void ReducedBasisSpace< BaseFunction >
  //     :: write ( OutStreamInterface< StreamTraits > &out ) const
  //   {
  //     const unsigned int size = numBaseFunctions();

  //     out << size;
  //     for( unsigned int i = 0; i < size; ++i )
  //       out << baseFunction( i );
  //   }



  // Stream Operators for ReducedBasisSpace
  // --------------------------------------

  /** \brief write a ReducedBasisSpace into an output stream
   *  \relates ReducedBasisSpace
   *
   *  \param[in]  out    stream to write to
   *  \param[in]  space  ReducedBasisSpace to write
   *
   *  \returns the output stream (for concatenation)
   */
  template< class StreamTraits, class BaseFunctionType >
  inline OutStreamInterface< StreamTraits > &
  operator<< ( OutStreamInterface< StreamTraits > &out,
               const ReducedBasisSpace< BaseFunctionType > &space )
  {
    space.write( out );
    return out;
  }


  /** \brief read a ReducedBasisSpace from an input stream
   *  \relates ReducedBasisSpace
   *
   *  \param[in]   in     stream to read from
   *  \param[out]  space  ReducedBasisSpace to read
   *
   *  \returns the input stream (for concatenation)
   */
  template< class StreamTraits, class BaseFunctionType >
  inline InStreamInterface< StreamTraits > &
  operator>> ( InStreamInterface< StreamTraits > &in,
               ReducedBasisSpace< BaseFunctionType > &space )
  {
    space.read( in );
    return in;
  }

} //end of namespace Dune::RB

} //end of namespace Dune

#endif
