// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_LINEARLAGRANGEINTERPOLATION_HH
#define DUNE_LINEARLAGRANGEINTERPOLATION_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/memory.hh>
#include <boost/math/special_functions/fpclassify.hpp>
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>

namespace Dune {
//! 2D:
//! easy way to obtain a linear function by either steting 3 point and associated values
//! or giving a picewiese linear discrete function on an entity.
//! this class provides a general evaluate method
//! NOTE if you want to use the method with the 'entity-Version' FunctionSpaceImp needs to be a DiscreteFunctionSpaceImp
//! \TODO this needs to be based on a localfunction implementation
//! \attention 2D Simplex only
template< class FunctionSpaceImp >
class LinearLagrangeFunction2D
  : public Dune::Fem::Function< FunctionSpaceImp, LinearLagrangeFunction2D< FunctionSpaceImp > >
{
private:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef LinearLagrangeFunction2D< FunctionSpaceType >      ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef FunctionSpaceType                                DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity                    EntityType;
  typedef typename EntityType::EntityPointer               EntityPointerType;

  typedef DomainFieldType TimeType;
  typedef std::unique_ptr<Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType>> LocalFunctionType;
  // three values that determine the linear polynomial in 2D
  DomainType a_0_;
  RangeType p_a_0_; // p(a_0) = p_a_0

  DomainType a_1_;
  RangeType p_a_1_;

  DomainType a_2_;
  RangeType p_a_2_;

  EntityPointerType* entity_;

  inline void check_dofs() const {
    if (boost::math::isnan(p_a_0_))
    {
      DSC_LOG_ERROR << "p_a_0 is nan" << std::endl;
    }
    if (boost::math::isnan(p_a_1_))
    {
      DSC_LOG_ERROR << "p_a_1 is nan" << std::endl;
    }
    if (boost::math::isnan(p_a_2_))
    {
      DSC_LOG_ERROR << "p_a_2 is nan" << std::endl;
    }
  }

public:
  // Constructor for LinearLagrangeFunction2D
  inline explicit LinearLagrangeFunction2D(DomainType a_0,
                                           RangeType p_a_0,
                                           DomainType a_1,
                                           RangeType p_a_1,
                                           DomainType a_2,
                                           RangeType p_a_2)
    : a_0_(a_0)
      , p_a_0_(p_a_0)
      , a_1_(a_1)
      , p_a_1_(p_a_1)
      , a_2_(a_2)
      , p_a_2_(p_a_2)
      , entity_(nullptr)
      , localFunc_(nullptr)
  {
    DUNE_THROW(NotImplemented, "This constructor is not implemented, yet!");
    check_dofs();
  }

  // Constructor for LinearLagrangeFunction2D
  inline explicit LinearLagrangeFunction2D(DiscreteFunctionSpaceType& dfSpace, EntityPointerType& entity)
    : localFunc_(DSC::make_unique(dfSpace, *entity))
  {}
  /*: a_0_( entity->geometry().corner(0) )
      , p_a_0_(0.0)
      , a_1_( entity->geometry().corner(1) )
      , p_a_1_(0.0)
      , a_2_( entity->geometry().corner(2) )
      , p_a_2_(0.0)
      , entity_(&entity)
  {}
*/
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    assert(localFunc_);
    localFunc_->evaluate(x, y);

    /*RangeType lambda_1;
    RangeType lambda_0;
    if (a_0_[0] == a_2_[0])
    {
      lambda_1 = (x[0] - a_2_[0]) / (a_1_[0] - a_2_[0]);

      lambda_0 = ( (x[1] - a_2_[1]) - ( lambda_1 * (a_1_[1] - a_2_[1]) ) ) / (a_0_[1] - a_2_[1]);
    } else {
      lambda_1
        = ( (x[1] - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (x[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) )
          / ( (a_1_[1] - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) );

      lambda_0
        = ( (x[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) - ( ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) * lambda_1 );
    }

    y = (p_a_0_ * lambda_0) + (p_a_1_ * lambda_1) + ( p_a_2_ * (1.0 - lambda_0 - lambda_1) );

    if (boost::math::isnan(lambda_1))
    {
      DSC_LOG_ERROR << "lambda_1 is nan! Details:" << std::endl << std::endl
                    << "a_0_ = " << a_0_ << std::endl
                    << "a_1_ = " << a_1_ << std::endl
                    << "a_2_ = " << a_2_ << std::endl << std::endl
                    << "p_a_0_ = " << p_a_0_ << std::endl
                    << "p_a_1_ = " << p_a_1_ << std::endl
                    << "p_a_2_ = " << p_a_2_ << std::endl << std::endl
                    << "a_0_[0] - a_2_[0] = " << a_0_[0] - a_2_[0] << std::endl
                    << "Whole denominator = "
                    << ( (a_1_[1]
                              - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) ) << std::endl;
      DUNE_THROW(Dune::InvalidStateException,"NaN");
    }

    if (boost::math::isnan(lambda_0))
    {
      DUNE_THROW(Dune::InvalidStateException,"lambda_0 is nan");
    }*/
  } // evaluate

  template< typename DiscreteFunctionType >
  inline void set_corners(const DiscreteFunctionType& disc_func) {
    auto loc_func = disc_func.localFunction( *(*entity_) );
    const int baseSetSize = localFunc_->baseFunctionSet().size();
    assert( baseSetSize == loc_func.baseFunctionSet().size());
    for (int i=0; i<baseSetSize; ++i) {
      (*localFunc_)[i] = loc_func[i];
    }

    /*const int number_of_nodes = ( *(*entity_) ).template count< 2 >();

    if ( number_of_nodes != int( loc_func.baseFunctionSet().size() ) ) {
      DSC_LOG_ERROR << "Error! Inconsistency in 'linear-lagrange-interpolation.hh'." << std::endl;
    }

    for (int i = 0; i < number_of_nodes; ++i) {
      const NodePointerType node = ( *(*entity_) ).template subEntity< 2 >(i);

      if ( node->geometry().corner(0) != (*entity_)->geometry().corner(i) ) {
        DSC_LOG_ERROR << "Error! Inconsistency in 'linear-lagrange-interpolation.hh'." << std::endl;
      }

      if (i == 0)
      { p_a_0_ = loc_func[0]; }

      if (i == 1)
      { p_a_1_ = loc_func[1]; }

      if (i == 2)
      { p_a_2_ = loc_func[2]; }
    }
    check_dofs();*/
  } // set_corners

  private:
    LocalFunctionType localFunc_;
};

//    void DUNE_DEPRECATED_MSG( "This Dune::Geometry is still a reference to its implementation." )
//    deprecationWarning ( integral_constant< bool, true > ) {}


//! DEPRECATED should be eventually replaced by the LinearLagrangeFunction2D implementation
template< class FunctionSpaceImp >
class LinearLagrangeInterpolation2D
  : public Dune::Fem::Function< FunctionSpaceImp, LinearLagrangeInterpolation2D< FunctionSpaceImp > >
{
private:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef LinearLagrangeInterpolation2D< FunctionSpaceType >      ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef FunctionSpaceType                                DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity                    EntityType;
  typedef typename EntityType::EntityPointer               EntityPointerType;

  typedef DomainFieldType TimeType;
  typedef std::unique_ptr<Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType>> LocalFunctionType;
  // three values that determine the linear polynomial in 2D
  DomainType a_0_;
  RangeType p_a_0_; // p(a_0) = p_a_0

  DomainType a_1_;
  RangeType p_a_1_;

  DomainType a_2_;
  RangeType p_a_2_;

  inline void check_dofs() const {
    if (boost::math::isnan(p_a_0_))
    {
      DSC_LOG_ERROR << "p_a_0 is nan" << std::endl;
    }
    if (boost::math::isnan(p_a_1_))
    {
      DSC_LOG_ERROR << "p_a_1 is nan" << std::endl;
    }
    if (boost::math::isnan(p_a_2_))
    {
      DSC_LOG_ERROR << "p_a_2 is nan" << std::endl;
    }
  }

  void DUNE_DEPRECATED_MSG( "This Dune::LinearLagrangeInterpolation2D is deprecated and should be replaced by Dune::LinearLagrangeFunction2D." )
    deprecationWarning ( ) {}
  
public:
  // Constructor for LinearLagrangeInterpolation2D
  inline explicit LinearLagrangeInterpolation2D(DomainType a_0,
                                           RangeType p_a_0,
                                           DomainType a_1,
                                           RangeType p_a_1,
                                           DomainType a_2,
                                           RangeType p_a_2)
    : a_0_(a_0)
      , p_a_0_(p_a_0)
      , a_1_(a_1)
      , p_a_1_(p_a_1)
      , a_2_(a_2)
      , p_a_2_(p_a_2)
  {
    deprecationWarning();
  }

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    RangeType lambda_1;
    RangeType lambda_0;
    if (a_0_[0] == a_2_[0])
    {
      lambda_1 = (x[0] - a_2_[0]) / (a_1_[0] - a_2_[0]);

      lambda_0 = ( (x[1] - a_2_[1]) - ( lambda_1 * (a_1_[1] - a_2_[1]) ) ) / (a_0_[1] - a_2_[1]);
    } else {
      lambda_1
        = ( (x[1] - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (x[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) )
          / ( (a_1_[1] - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) );

      lambda_0
        = ( (x[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) - ( ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) * lambda_1 );
    }

    y = (p_a_0_ * lambda_0) + (p_a_1_ * lambda_1) + ( p_a_2_ * (1.0 - lambda_0 - lambda_1) );

    if (boost::math::isnan(lambda_1))
    {
      DSC_LOG_ERROR << "lambda_1 is nan! Details:" << std::endl << std::endl
                    << "a_0_ = " << a_0_ << std::endl
                    << "a_1_ = " << a_1_ << std::endl
                    << "a_2_ = " << a_2_ << std::endl << std::endl
                    << "p_a_0_ = " << p_a_0_ << std::endl
                    << "p_a_1_ = " << p_a_1_ << std::endl
                    << "p_a_2_ = " << p_a_2_ << std::endl << std::endl
                    << "a_0_[0] - a_2_[0] = " << a_0_[0] - a_2_[0] << std::endl
                    << "Whole denominator = "
                    << ( (a_1_[1]
                              - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) ) << std::endl;
      DUNE_THROW(Dune::InvalidStateException,"NaN");
    }

    if (boost::math::isnan(lambda_0))
    {
      DUNE_THROW(Dune::InvalidStateException,"lambda_0 is nan");
    }
  } // evaluate


};



} // end namespace
#endif // end LINEARLAGRANGEINTERPOLATION
