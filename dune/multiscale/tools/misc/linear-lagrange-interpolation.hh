/**************************************************************************
**       Title: L2Error class
**    $RCSfile$
**   $Revision: 1723 $$Name$
**       $Date: 2007-06-20 15:20:54 +0000 (Wed, 20 Jun 2007) $
**   Copyright: GPL $Author: dune community $
** Description: L2 error class, which computes the error between a function
**              and a discrete function. Extracted class from
**              Roberts poisson-example 
**
**************************************************************************/

#ifndef DUNE_LINEARLAGRANGEINTERPOLATION_HH
#define DUNE_LINEARLAGRANGEINTERPOLATION_HH



namespace Dune 
{

  // 2D:
  // easy way to obtain a linear function by either steting 3 point and associated values
  // or giving a picewiese linear discrete function on an entity.
  // this class provides a general evaluate method

  // NOTE if you want to use the method with the 'entity-Version' FunctionSpaceImp needs to be a DiscreteFunctionSpaceImp
  template< class FunctionSpaceImp >
  class LinearLagrangeFunction2D
  : public Function< FunctionSpaceImp, LinearLagrangeFunction2D< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef LinearLagrangeFunction2D< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef FunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;
    typedef typename EntityType :: EntityPointer EntityPointerType; 

    typedef DomainFieldType TimeType;

  protected:
    // three values that determine the linear polynom in 2D
    DomainType a_0_;
    RangeType p_a_0_; // p(a_0) = p_a_0

    DomainType a_1_;
    RangeType p_a_1_;

    DomainType a_2_;
    RangeType p_a_2_;

  private:

   EntityPointerType* it_;

  public:
    // Constructor for LinearLagrangeFunction2D
    inline explicit LinearLagrangeFunction2D ( DomainType a_0,
                                               RangeType p_a_0,
                                               DomainType a_1,
                                               RangeType p_a_1,
                                               DomainType a_2,
                                               RangeType p_a_2 )
    : a_0_( a_0 ),
      p_a_0_( p_a_0 ),
      a_1_( a_1 ),
      p_a_1_( p_a_1 ),
      a_2_( a_2 ),
      p_a_2_( p_a_2 ),
      it_ ( NULL )
    {
    }

    // Constructor for LinearLagrangeFunction2D
    inline explicit LinearLagrangeFunction2D ( EntityPointerType& it )
    : a_0_( it->geometry().corner(0) ),
      p_a_0_( 0.0 ),
      a_1_( it->geometry().corner(1) ),
      p_a_1_( 0.0 ),
      a_2_( it->geometry().corner(2) ),
      p_a_2_( 0.0 ),
      it_ ( &it )
    {
    }

    inline void evaluate ( const DomainType &x,
                           RangeType &y ) const
    {

      if ( p_a_0_ != p_a_0_ )
        {std :: cout << "p_a_0 is nan" << std :: endl;
         }//!abort();}

      if ( p_a_1_ != p_a_1_ )
        {std :: cout << "p_a_1 is nan" << std :: endl;
         }//!abort();}

      if ( p_a_2_ != p_a_2_ )
        {std :: cout << "p_a_2 is nan" << std :: endl;
         }//!abort();}

      if ( x[0] != x[0] )
        {std :: cout << "x[0] is nan" << std :: endl;
         }//!abort();}

      if ( x[1] != x[1] )
        {std :: cout << "x[1] is nan" << std :: endl;
         }//!abort();}

      RangeType lambda_1;
      RangeType lambda_0;

      if (a_0_[0] == a_2_[0])
       {
         lambda_1 = (x[0] - a_2_[0]) / (a_1_[0] - a_2_[0]);

         lambda_0 = ( (x[1] - a_2_[1]) - ( lambda_1 * (a_1_[1] - a_2_[1]) ) ) / (a_0_[1] - a_2_[1]);
       }
      else
       {
         lambda_1 = 
                ( (x[1]    - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (x[0]    - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) )
              / ( (a_1_[1] - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) );

         lambda_0 =
                ( ( x[0] - a_2_[0] ) / ( a_0_[0] - a_2_[0] ) ) - ( ( ( a_1_[0] - a_2_[0] ) / (a_0_[0] - a_2_[0] ) ) * lambda_1 );

       }



      y = (p_a_0_ * lambda_0) + (p_a_1_ * lambda_1) + (p_a_2_ * (1.0 - lambda_0 - lambda_1) );

      if ( lambda_1 != lambda_1 )
        {
         std :: cout << "lambda_1 is nan! Details:" << std :: endl << std :: endl;
         std :: cout << "a_0_ = " << a_0_ << std :: endl;
         std :: cout << "a_1_ = " << a_1_ << std :: endl;
         std :: cout << "a_2_ = " << a_2_ << std :: endl << std :: endl;
         std :: cout << "p_a_0_ = " << p_a_0_ << std :: endl;
         std :: cout << "p_a_1_ = " << p_a_1_ << std :: endl;
         std :: cout << "p_a_2_ = " << p_a_2_ << std :: endl << std :: endl;
         std :: cout << "a_0_[0] - a_2_[0] = " << a_0_[0] - a_2_[0] << std :: endl;
         std :: cout << "ganzer Nenner = " << ( (a_1_[1] - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) ) << std :: endl;
         abort();}

      if ( lambda_0 != lambda_0 )
        {std :: cout << "lambda_0 is nan" << std :: endl;
         abort();}

    }

    template< typename DiscreteFunctionType >
    inline void set_corners( const DiscreteFunctionType& disc_func )
    {
std :: cout << "HIER Line 167" << std :: endl;
      typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
      typedef typename EntityType :: template Codim< 2 > :: EntityPointer NodePointerType;

      LocalFunctionType loc_func = disc_func.localFunction( *(*it_) );
std :: cout << "HIER Line 172" << std :: endl;

      int number_of_nodes = (*(*it_)).template count<2>();
std :: cout << "HIER Line 175" << std :: endl;
      if (!( number_of_nodes == loc_func.baseFunctionSet().numBaseFunctions() ))
       { std :: cout << "Error! Inconsistency in 'linear-lagrange-interpolation.hh'." << std :: endl; }

std :: cout << "HIER Line 179" << std :: endl;
      for ( int i = 0; i < number_of_nodes; i += 1 )
       {
std :: cout << "HIER Line 182" << std :: endl;
          const NodePointerType node = (*(*it_)).template subEntity<2>(i);
std :: cout << "HIER Line 184" << std :: endl;

#if 1
          if ( !(node->geometry().corner(0) == (*it_)->geometry().corner(i)) )
           { std :: cout << "Error! Inconsistency in 'linear-lagrange-interpolation.hh'." << std :: endl; }
#endif

          if ( i == 0 )
           { p_a_0_ = loc_func[0]; }

          if ( i == 1 )
           { p_a_1_ = loc_func[1]; }

          if ( i == 2 )
           { p_a_2_ = loc_func[2]; }

      }

    }

  };


} // end namespace 
#endif //end LINEARLAGRANGEINTERPOLATION
