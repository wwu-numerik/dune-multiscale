// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ANALYTICALHOMOGENIZER_HH
#define DUNE_ANALYTICALHOMOGENIZER_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#else
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

// where the quadratures are defined
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/stuff/common/logging.hh>

namespace Dune {

//! das alles klappt (mathematisch) nur in 2-D!!! Für Tensoren, die:
//! 1. irgendwelche Elliptizitätsbedingungen erfüllen
//! 2. A(x,y) = A(y)
//! 3. a_i_j(y) = a_i_j(y_1,y_2) = a_i_j(y_1)
//! 4. symmetrisch
template< class GridImp, class TensorImp >
class AnalyticalHomogenizer
{
private:
  typedef GridImp GridType;

  typedef LeafGridPart< GridType > GridPartType;

  typedef FunctionSpace< double, double, 2, 1 > FunctionSpaceType;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >
  DiscreteFunctionSpaceType;

  typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef TensorImp TensorType;

  typedef typename DiscreteFunctionSpaceType::IteratorType
  IteratorType;

  typedef typename GridType::template Codim< 0 >::Entity
  EntityType;

  typedef typename GridType::template Codim< 0 >::Geometry
  EnGeometryType;

  typedef typename EntityType::ctype
  coordType;

  enum { dimension = GridType::dimension };
  enum { spacePolOrd = DiscreteFunctionSpaceType::polynomialOrder };

  typedef FieldMatrix< RangeType, dimension, dimension > TensorMatrixType;


  // dgf file that describes the perforated domain
  std::string& filename_;

public:
  AnalyticalHomogenizer(std::string& filename)
    : filename_(filename)
  {}

  FieldMatrix< RangeType, dimension, dimension > getHomTensor
    ( const TensorType& tensor, int polOrd = (2 * spacePolOrd + 2) ) const {
    std::cout
    << "WARNING! Use of deprecated homogenizer, which requires deprecated use of 'evaluate' in Diffusion-class!"
    << std::endl;

    GridPtr< GridType > gridptr(filename_);
    gridptr->globalRefine(12);
    GridPartType gridPart(*gridptr);

    DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);

    FieldMatrix< RangeType, dimension, dimension > tensorHom(0.0);

    TensorMatrixType a;
    for (const EntityType& entity : discreteFunctionSpace)
    {
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geometry = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geometry.integrationElement( quadrature.point(quadraturePoint) );

        tensor.evaluate(0, 0,
                        geometry.global( quadrature.point(quadraturePoint) ),
                        a[0][0]);

        tensorHom[0][0] += det * (1 / a[0][0]);
      }
    }

    tensorHom[0][0] = 1 / tensorHom[0][0];

    for (const EntityType& entity : discreteFunctionSpace)
    {
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geometry = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geometry.integrationElement( quadrature.point(quadraturePoint) );

        tensor.evaluate(0, 0,
                        geometry.global( quadrature.point(quadraturePoint) ),
                        a[0][0]);

        tensor.evaluate(0, 1,
                        geometry.global( quadrature.point(quadraturePoint) ),
                        a[0][1]);

        tensorHom[0][1] += det * (a[0][1] / a[0][0]);
      }
    }

    tensorHom[0][1] *= tensorHom[0][0];
    tensorHom[1][0] = tensorHom[0][1];

    for (const EntityType& entity : discreteFunctionSpace)
    {
      const CachingQuadrature< GridPartType, 0 > quadrature(entity, polOrd);

      // get geoemetry of entity
      const EnGeometryType& geometry = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * geometry.integrationElement( quadrature.point(quadraturePoint) );

        tensor.evaluate(0, 0,
                        geometry.global( quadrature.point(quadraturePoint) ),
                        a[0][0]);

        tensor.evaluate(0, 1,
                        geometry.global( quadrature.point(quadraturePoint) ),
                        a[0][1]);

        tensor.evaluate(1, 0,
                        geometry.global( quadrature.point(quadraturePoint) ),
                        a[1][0]);

        tensor.evaluate(1, 1,
                        geometry.global( quadrature.point(quadraturePoint) ),
                        a[1][1]);

        tensorHom[1][1] += det
                           * ( a[1][1] - ( (a[0][1] * a[1][0]) / a[0][0] ) );
      }
    }

    tensorHom[1][1] += (tensorHom[1][0] * tensorHom[0][1])
                       / tensorHom[0][0];

    DSC_LOG_DEBUG << "analytical: A_homogenized[0][0] = " << tensorHom[0][0] << std::endl
                  << "analytical: A_homogenized[0][1] = " << tensorHom[0][1] << std::endl
                  << "analytical: A_homogenized[1][0] = " << tensorHom[1][0] << std::endl
                  << "analytical: A_homogenized[1][1] = " << tensorHom[1][1] << std::endl;

    return tensorHom;
  } // end of method
}; // end of class
} // end namespace
#endif // ifndef DUNE_ANALYTICALHOMOGENIZER_HH
