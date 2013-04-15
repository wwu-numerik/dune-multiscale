#ifndef DUNE_LINEARLAGRANGEINTERPOLATION_HH
#define DUNE_LINEARLAGRANGEINTERPOLATION_HH

#include <dune/fem/function/common/function.hh>
#include <dune/stuff/common/logging.hh>

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

  // three values that determine the linear polynom in 2D
  DomainType a_0_;
  RangeType p_a_0_; // p(a_0) = p_a_0

  DomainType a_1_;
  RangeType p_a_1_;

  DomainType a_2_;
  RangeType p_a_2_;

  EntityPointerType* entity_;

  inline void check_dofs() const {
    if (std::isnan(p_a_0_))
    {
      DSC_LOG_ERROR << "p_a_0 is nan" << std::endl;
    }
    if (std::isnan(p_a_1_))
    {
      DSC_LOG_ERROR << "p_a_1 is nan" << std::endl;
    }
    if (std::isnan(p_a_2_))
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
  {
    check_dofs();
  }

  // Constructor for LinearLagrangeFunction2D
  inline explicit LinearLagrangeFunction2D(EntityPointerType& entity)
    : a_0_( entity->geometry().corner(0) )
      , p_a_0_(0.0)
      , a_1_( entity->geometry().corner(1) )
      , p_a_1_(0.0)
      , a_2_( entity->geometry().corner(2) )
      , p_a_2_(0.0)
      , entity_(&entity)
  {}

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

    if (std::isnan(lambda_1))
    {
      DSC_LOG_ERROR << "lambda_1 is nan! Details:" << std::endl << std::endl
                    << "a_0_ = " << a_0_ << std::endl
                    << "a_1_ = " << a_1_ << std::endl
                    << "a_2_ = " << a_2_ << std::endl << std::endl
                    << "p_a_0_ = " << p_a_0_ << std::endl
                    << "p_a_1_ = " << p_a_1_ << std::endl
                    << "p_a_2_ = " << p_a_2_ << std::endl << std::endl
                    << "a_0_[0] - a_2_[0] = " << a_0_[0] - a_2_[0] << std::endl
                    << "ganzer Nenner = "
                    << ( (a_1_[1]
                              - a_2_[1]) - ( (a_0_[1] - a_2_[1]) * ( (a_1_[0] - a_2_[0]) / (a_0_[0] - a_2_[0]) ) ) ) << std::endl;
      DUNE_THROW(Dune::InvalidStateException,"NaN");
    }

    if (std::isnan(lambda_0))
    {
      DUNE_THROW(Dune::InvalidStateException,"lambda_0 is nan");
    }
  } // evaluate

  template< typename DiscreteFunctionType >
  inline void set_corners(const DiscreteFunctionType& disc_func) {
    typedef typename DiscreteFunctionType::LocalFunctionType        LocalFunctionType;
    typedef typename EntityType::template Codim< 2 >::EntityPointer NodePointerType;

    LocalFunctionType loc_func = disc_func.localFunction( *(*entity_) );

    const int number_of_nodes = ( *(*entity_) ).template count< 2 >();

    if ( !( number_of_nodes == int( loc_func.baseFunctionSet().size() ) ) )
    { DSC_LOG_ERROR << "Error! Inconsistency in 'linear-lagrange-interpolation.hh'." << std::endl; }

    for (int i = 0; i < number_of_nodes; i += 1)
    {
      const NodePointerType node = ( *(*entity_) ).template subEntity< 2 >(i);

      if ( !( node->geometry().corner(0) == (*entity_)->geometry().corner(i) ) )
      { DSC_LOG_ERROR << "Error! Inconsistency in 'linear-lagrange-interpolation.hh'." << std::endl; }

      if (i == 0)
      { p_a_0_ = loc_func[0]; }

      if (i == 1)
      { p_a_1_ = loc_func[1]; }

      if (i == 2)
      { p_a_2_ = loc_func[2]; }
    }
    check_dofs();
  } // set_corners
};
} // end namespace
#endif // end LINEARLAGRANGEINTERPOLATION
