#ifndef H1_ERROR_HH
#define H1_ERROR_HH

// where the quadratures are defined
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/l2error.hh>

namespace Dune {
template< class DiscreteFunctionType, int n = 0 >
class H1Error
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::RangeType                 RangeType;
  typedef typename DiscreteFunctionType::DomainType                DomainType;
  typedef typename DiscreteFunctionType::RangeFieldType            RangeFieldType;
  typedef typename DiscreteFunctionType::JacobianRangeType         JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::IteratorType         IteratorType;

  typedef typename DiscreteFunctionSpaceType::GridType     GridType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef typename GridType::template Codim< 0 >::Entity        EntityType;
  typedef typename GridType::template Codim< 0 >::EntityPointer EntityPointerType;
  typedef typename GridType::template Codim< 0 >::Geometry      EnGeometryType;
  typedef typename GridType::Traits::GlobalIdSet                IdSet;

  typedef typename IdSet::IdType IdType;

  typedef typename EntityType::ctype coordType;

  // hierarchic iterator over all elements of codim 0
  typedef typename EntityType::HierarchicIterator               HierarchicIteratorType;
  typedef typename GridType::template Codim< 0 >::LevelIterator LevelIteratorType;

  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  enum { dim = GridType::dimension };
  enum { spacePolOrd = DiscreteFunctionSpaceType::polynomialOrder };
  enum { dimRange = RangeType::dimension };

public:
  // H1 error:
  template< class FunctionType >
  RangeType semi_norm(const FunctionType& function,
                      const DiscreteFunctionType& disc_function,
                      const int polOrd = (2 * spacePolOrd + 2) ) const {

    RangeType y(0.0); // return value

    const IteratorType endit = disc_function.space().end();
    for (IteratorType it = disc_function.space().begin(); it != endit; ++it)
    {
      // entity
      const EntityType& entity = *it;

      // create quadrature for given geometry type
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geo = entity.geometry();

      const LocalFunctionType local_disc_func = disc_function.localFunction(entity);

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        JacobianRangeType gradient_of_function;
        JacobianRangeType gradient_of_disc_function;

        local_disc_func.jacobian(quadrature[quadraturePoint], gradient_of_disc_function);

        function.evaluateJacobian(geo.global( quadrature.point(quadraturePoint) ), gradient_of_function);

        const double det = quadrature.weight(quadraturePoint)
                           * geo.integrationElement( quadrature.point(quadraturePoint) );

        RangeType value(0.0);
        for (int i = 0; i < dim; ++i)
        { value += pow(gradient_of_function[0][i] - gradient_of_disc_function[0][i], 2.0); }

        y += det * value;
      }
    }

    return sqrt(y);
  } // end of method
};
}
#endif // ifndef H1_ERROR_HH
