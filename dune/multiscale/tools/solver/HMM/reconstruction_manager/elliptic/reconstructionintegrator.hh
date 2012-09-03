/**************************************************************************
  **       Title: RECONSTRUCTIONINTEGRATER
  **    $RCSfile$
  **   $Revision: 1723 $$Name$
  **       $Date: 2007-06-20 15:20:54 +0000 (Wed, 20 Jun 2007) $
  **   Copyright: GPL $Author: dune community $
  ** Description: Neue Beschreibung machen!!!!!!!!
  **
  **
  **
  **************************************************************************/

#ifndef DUNE_RECONSTRUCTIONINTEGRATER_HH
#define DUNE_RECONSTRUCTIONINTEGRATER_HH

// where the quadratures are defined
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {

/**
 *  \class RecInt
 *  \brief The RecInt class provides methods to calculate the meanvalues of local reconstructions of base functions
 * \attention method adapt doesn't work with ALBERTAGRID, because leakpointer produceses errors!!!
 */
template< class PeriodicDiscreteFunctionImp, class TensorType >
class RecInt // Reconstruction Integrator
{
  //! type of discrete functions
  typedef PeriodicDiscreteFunctionImp PeriodicDiscreteFunctionType;

  typedef typename PeriodicDiscreteFunctionImp::LocalFunctionType PeriodicLocalFunctionType;

  //! type of discrete function space
  typedef typename PeriodicDiscreteFunctionImp::DiscreteFunctionSpaceType
  PeriodicDiscreteFunctionSpaceType;

  //! type of grid partition
  typedef typename PeriodicDiscreteFunctionSpaceType::GridPartType PeriodicGridPartType;

  //! type of grid
  typedef typename PeriodicDiscreteFunctionSpaceType::GridType PeriodicGridType;

  typedef typename PeriodicDiscreteFunctionType::RangeType
  RangeType;

  typedef typename PeriodicDiscreteFunctionType::DomainType
  DomainType;

  typedef typename PeriodicDiscreteFunctionSpaceType::JacobianRangeType
  PeriodicJacobianRangeType;

  typedef typename PeriodicDiscreteFunctionSpaceType::IteratorType
  PeriodicIteratorType;

  typedef typename PeriodicGridPartType::IntersectionIteratorType PeriodicIntersectionIteratorType;

  typedef typename PeriodicGridType::template Codim< 0 >::Entity
  PeriodicEntityType;

  typedef typename PeriodicGridType::template Codim< 0 >::EntityPointer
  PeriodicEntityPointerType;

  typedef typename PeriodicGridType::template Codim< 0 >::Geometry
  PeriodicEntityGeometryType;

  typedef typename PeriodicGridType::template Codim< 1 >::Geometry
  PeriodicFaceGeometryType;

  typedef CachingQuadrature< PeriodicGridPartType, 0 > PeriodicEntityQuadratureType;

  typedef CachingQuadrature< PeriodicGridPartType, 1 > PeriodicFaceQuadratureType;

  enum { dimension = PeriodicGridType::dimension };
  enum { spacePolOrd = PeriodicDiscreteFunctionSpaceType::polynomialOrder };

public:
  template< class JacobianRangeImp >
  RangeType integrateGlobalBaseFunctions( const TensorType& tensor,
                                          const DomainType& globalPoint,
                                          const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                                          const JacobianRangeImp& grad_PHI_i,
                                          const JacobianRangeImp& grad_PHI_j,
                                          const int polOrd = (2 * spacePolOrd + 2) ) const {
    // Note that this method does not work for perforated structures!
    // (but the old version of the code also works for the case:
    // 1. A^{eps}(x) = A(x,x/eps) with A(x,.) Y-periodic AND
    // 2. Y* is non-perorated or periodically perforated )
    // (if required, ask Patrick Henning for the old version of the code)

    const double epsilon_est = DSC_CONFIG_GET("problem.epsilon_guess", 1.0f);
    const PeriodicGridPartType& gridPart = periodicDiscreteFunctionSpace.gridPart();
    typedef typename PeriodicGridPartType::GridType::Traits::CollectiveCommunication
    CommunicatorType;

    const CommunicatorType& comm = gridPart.grid().comm();
    RangeType result(0.0);

    const PeriodicIteratorType endit = periodicDiscreteFunctionSpace.end();
    for (PeriodicIteratorType it = periodicDiscreteFunctionSpace.begin(); it != endit; ++it)
    {
      // entity
      const PeriodicEntityType& entity = *it;
      // create quadrature for given geometry type
      const PeriodicEntityQuadratureType quadrature(entity, polOrd);
      // get geoemetry of entity
      const PeriodicEntityGeometryType& geometry = entity.geometry();
      // integrate
      const int quadratureNop = quadrature.nop();
      for (int localQuadPoint = 0; localQuadPoint < quadratureNop; ++localQuadPoint)
      {
        RangeType localIntegral = 0;
        RangeType a[dimension][dimension];

        DomainType y_eps; // x_j + \eps_{est} * y
        for (int k = 0; k < dimension; ++k)
          y_eps[k] = globalPoint[k] + epsilon_est* geometry.global( quadrature.point(localQuadPoint) )[k];

        for (int k = 0; k < dimension; ++k)
          for (int l = 0; l < dimension; ++l)
            tensor.evaluate(k, l, y_eps, a[k][l]);

        RangeType w[dimension];
        for (int k = 0; k < dimension; ++k)
          w[k] = 0;

        for (int k = 0; k < dimension; ++k)
          for (int l = 0; l < dimension; ++l)
            w[k] += a[k][l] * grad_PHI_i[0][l];

        for (int k = 0; k < dimension; ++k)
          localIntegral += w[k] * grad_PHI_j[0][k];

        const double entityVolume = quadrature.weight(localQuadPoint)
                                    * geometry.integrationElement( quadrature.point(localQuadPoint) );

        result += entityVolume * localIntegral;
      }
    }

    RangeType cell_volume = 1;

    for (int k = 0; k < dimension; ++k)
      cell_volume *= epsilon_est;

    result = cell_volume * comm.sum(result);

    return result;
  } // end of method

  template< class JacobianRangeImp >
  RangeType integrateCorrectorBaseFunctions
    ( const TensorType& tensor,
    const DomainType& globalPoint,
    const PeriodicDiscreteFunctionSpaceType
    & periodicDiscreteFunctionSpace,
    const PeriodicDiscreteFunctionType& corrector_PHI_i,
    const JacobianRangeImp& grad_PHI_j,
    const int polOrd = (2 * spacePolOrd + 2) ) const {
    // Note that this method does not work for perforated structures!
    // (but the old version of the code also works for the case:
    // 1. A^{eps}(x) = A(x,x/eps) with A(x,.) Y-periodic AND
    // 2. Y* is non-perorated or periodically perforated )
    // (if required, ask Patrick Henning for the old version of the code)

    const double epsilon_est = DSC_CONFIG_GET("problem.epsilon_guess", 1.0f);
    const PeriodicGridPartType& gridPart = periodicDiscreteFunctionSpace.gridPart();

    typedef typename PeriodicGridPartType::GridType::Traits::
      CollectiveCommunication
    CommunicatorType;

    const CommunicatorType& comm = gridPart.grid().comm();
    RangeType result(0.0);

    const PeriodicIteratorType endit = periodicDiscreteFunctionSpace.end();
    for (PeriodicIteratorType it = periodicDiscreteFunctionSpace.begin(); it != endit; ++it)
    {
      // entity
      const PeriodicEntityType& entity = *it;

      const PeriodicLocalFunctionType localfunc = corrector_PHI_i.localFunction(entity);

      // create quadrature for given geometry type
      PeriodicEntityQuadratureType quadrature(entity, polOrd);

      // get geoemetry of entity
      const PeriodicEntityGeometryType& geometry = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int localQuadPoint = 0; localQuadPoint < quadratureNop; ++localQuadPoint)
      {
        RangeType localIntegral = 0;

        PeriodicJacobianRangeType gradLocCor;
        localfunc.jacobian(quadrature[localQuadPoint], gradLocCor);
        // In comparison to the jacobian method for base function we do not need an additional transformation of the
        // gradient. This is due to the fact that the discrete functions are global functions whereas base functions
        // life
        // on the reference element!

        RangeType a[dimension][dimension];

        DomainType y_eps; // x_j + \eps_{est} * y
        for (int k = 0; k < dimension; ++k)
          y_eps[k] = globalPoint[k] + epsilon_est* geometry.global( quadrature.point(localQuadPoint) )[k];

        for (int k = 0; k < dimension; ++k)
          for (int i = 0; i < dimension; ++i)
            tensor.evaluate(i, k, y_eps, a[i][k]);

        RangeType w[dimension];
        for (int k = 0; k < dimension; ++k) w[k] = 0;

        for (int k = 0; k < dimension; ++k)
          for (int l = 0; l < dimension; ++l)
            w[k] += a[k][l] * gradLocCor[0][l];

        for (int k = 0; k < dimension; ++k)
          localIntegral += w[k] * grad_PHI_j[0][k];

        const double entityVolume = quadrature.weight(localQuadPoint)
                                    * geometry.integrationElement( quadrature.point(localQuadPoint) );

        result += entityVolume * localIntegral;
      }
    }

    RangeType epsweight = 1; // to calculate \epsilon^(dimension - 1):
    for (int k = 0; k < (dimension - 1); ++k)
      epsweight *= epsilon_est;

    result = epsweight * comm.sum(result);

    return result;
  } // end of method
}; // end of class RecInt
} // end namespace
#endif // ifndef DUNE_RECONSTRUCTIONINTEGRATER_HH
