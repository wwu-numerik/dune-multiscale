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

#ifndef DUNE_MS_MEANVALUE_HH
#define DUNE_MS_MEANVALUE_HH

// where the quadratures are defined
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/l2error.hh>

#include "misc/linear-lagrange-interpolation.hh"

namespace Dune {
/**  \class Meanvalue
   *  \brief The Meanvalue class provides a method to calculate the meanvalue of a discrete function
   *
   *  Actually only the meanvalue of discrete functions on the unit-cube can be calculated.
   *  If you want more, divide the return value of getMeanvalue() by the size of the domain.
   *  Since it's the unit cube for our purpose, the domain-size is 1 and therefore unimportant
   *  for us.
   *
 * Usage of the meanvalue class:
  * Example:
   *
   *  Meanvalue< DiscreteFunctionType > mymean;
   *  DiscreteFunctionSpaceType :: RangeType theMeanValue;
   *  theMeanValue = mymean.getMeanvalue( a_discreteFunction );
   *  std :: cout << "Meanvalue of the numerical solution: " << theMeanValue << std :: endl;
   *
   * Shift the discrete function to meanvalue zero:
   *
   * mymean.adapt<DiscreteFunctionType>( a_discreteFunction , theMeanValue );
   **/
template< class DiscreteFunctionType >
class Meanvalue
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionType::RangeType
  RangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType
  IteratorType;

  typedef typename DiscreteFunctionSpaceType::GridType
  GridType;

  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef DomainFieldType                                     TimeType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
  GridPartType;

  typedef typename GridType::template Codim< 0 >::Entity
  EntityType;

  typedef typename GridType::template Codim< 0 >::Geometry
  EnGeometryType;

  typedef typename EntityType::ctype
  coordType;

  typedef typename DiscreteFunctionType::LocalFunctionType
  LocalFunctionType;

  enum { dimension = GridType::dimension };
  enum { spacePolOrd = DiscreteFunctionSpaceType::polynomialOrder };

public:
  RangeType getMeanvalue(const DiscreteFunctionType& discFunc) const {
    int polOrd = (2 * spacePolOrd + 2);

    // get function space
    const DiscreteFunctionSpaceType& space = discFunc.space();
    const GridPartType& gridPart = space.gridPart();
    const auto& comm = gridPart.grid().comm();

    RangeType y(0.0);    // return value
    RangeType theMeanValue(0.0);

    for (const auto& entity : space)
    {
      // create quadrature for given geometry type
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get local function
      LocalFunctionType localfunc = discFunc.localFunction(entity);

      // get geoemetry of entity
      const EnGeometryType& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geo.integrationElement( quadrature.point(quadraturePoint) );

        localfunc.evaluate(quadrature, quadraturePoint, y);

        theMeanValue += det * y;
      }
    }

    theMeanValue = comm.sum(theMeanValue);

    return theMeanValue;
  }   // end of method

  template< class FunctionType >
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space,
                         const FunctionType& function) const {
    int polOrd = (2 * spacePolOrd + 2);

    RangeType y(0.0);    // return value
    RangeType theMeanValue(0.0);

    for (const auto& entity : space)
    {
      // create quadrature for given geometry type
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geo.integrationElement( quadrature.point(quadraturePoint) );

        function.evaluate(geo.global( quadrature.point(quadraturePoint) ), y);

        theMeanValue += det * y;
      }
    }

    return theMeanValue;
  }   // end of method

  // the case the function is a vector (for instance advection)
  template< class FunctionType >
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space,
                         const FunctionType& function,
                         const int& i /*in case there are several components*/) const {
    int polOrd = (2 * spacePolOrd + 2);

    RangeType y(0.0);    // return value
    RangeType theMeanValue(0.0);

    for (const auto& entity : space)
    {
      // create quadrature for given geometry type
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geo.integrationElement( quadrature.point(quadraturePoint) );

        function.evaluate(i, geo.global( quadrature.point(quadraturePoint) ), y);

        theMeanValue += det * y;
      }
    }

    return theMeanValue;
  }   // end of method

  // the case the function is a time-dependent vector (for instance advection). The Time t is fixed.
  template< class FunctionType >
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space,
                         const FunctionType& function,
                         const TimeType& t,
                         const int& i /*in case there are several components*/) const {
    int polOrd = (2 * spacePolOrd + 2);

    RangeType y(0.0);    // return value
    RangeType theMeanValue(0.0);

    for (const auto& entity : space)
    {
      // create quadrature for given geometry type
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geo.integrationElement( quadrature.point(quadraturePoint) );

        function.evaluate(i, geo.global( quadrature.point(quadraturePoint) ), t, y);

        theMeanValue += det * y;
      }
    }

    return theMeanValue;
  }   // end of method

  // the case the function is a matrix (for instance diffusion)
  template< class FunctionType >
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space,
                         const FunctionType& function,
                         const int& i,
                         const int& j) const {
    int polOrd = (2 * spacePolOrd + 2);

    RangeType y(0.0);    // return value
    RangeType theMeanValue(0.0);

    for (const auto& entity : space)
    {
      // create quadrature for given geometry type
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geo.integrationElement( quadrature.point(quadraturePoint) );

        function.evaluate(i, j, geo.global( quadrature.point(quadraturePoint) ), y);

        theMeanValue += det * y;
      }
    }

    return theMeanValue;
  }   // end of method

  // the case the function is a matrix (for instance diffusion) with Time t
  template< class FunctionType >
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space,
                         const FunctionType& function,
                         const TimeType& t,
                         const int& i,
                         const int& j) const {
    int polOrd = (2 * spacePolOrd + 2);

    RangeType y(0.0);    // return value
    RangeType theMeanValue(0.0);

    for (const auto& entity : space)
    {
      // create quadrature for given geometry type
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geo.integrationElement( quadrature.point(quadraturePoint) );

        function.evaluate(i, j, geo.global( quadrature.point(quadraturePoint) ), t, y);

        theMeanValue += det * y;
      }
    }

    return theMeanValue;
  }   // end of method

  //! Subdraktion des Mittelwertes von der DiscreteFunction
  template< class FunctionType >
  static void adapt(FunctionType& discreteFunction, RangeType& meanvalue) {
    typedef typename FunctionType::DofIteratorType DofIteratorType;

    const DofIteratorType end = discreteFunction.dend();
    for (DofIteratorType it = discreteFunction.dbegin(); it != end; ++it)
      *it -= meanvalue[0];
    // Dof Iterator verwenden um Verschiebung um Mittelwert

  } // adapt
}; // end of class Meanvalue

} // end namespace DUNE

#endif // ifndef DUNE_MS_MEANVALUE_HH
