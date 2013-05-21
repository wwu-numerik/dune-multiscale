// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_COMMON_TRAITS_HH
#define DUNE_MULTISCALE_COMMON_TRAITS_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/tuples.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {

template <class T>
class GridPtr;
template <class T, class R>
class GridFunctionAdapter;
template <class T>
class LagrangeParallelMatrixAdapter;
template <class T>
class ParallelScalarProduct;
template <class T, class R>
class DataOutput;
template <bool T>
class LagrangeMatrixSetup;
template <class T, class R>
class AdaptationManager;
template <class T, class R, class S>
class SparseRowMatrixOperator;

namespace Multiscale {

namespace Problem {
namespace PROBLEM_NAME {
class ModelProblemData;
class FirstSource;
class SecondSource;
class Diffusion;
class MassTerm;
class DefaultDummyFunction;
class ExactSolution;
} //namespace PROBLEM_NAME
} //namespace Problem

//! type construction for the HMM algorithm
struct CommonTraits {
  //! --------- typedefs for the macro grid and the corresponding discrete space -------------
  typedef Dune::GridSelector::GridType
      GridType;
  // Dune::InteriorBorder_Partition or Dune::All_Partition >?
  // see:
  // http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
  typedef Dune::AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;
  typedef Dune::GridPtr< GridType >                                       GridPointerType;
  typedef Dune::FunctionSpace< double, double, WORLDDIM, 1 >              FunctionSpaceType;
  //!-----------------------------------------------------------------------------------------

  //! --------- typedefs for the coefficient and data functions ------------------------------
  typedef Problem::PROBLEM_NAME::ModelProblemData ModelProblemDataType;
  // type of first source term (right hand side of differential equation or type of 'f')
  typedef Problem::PROBLEM_NAME::FirstSource FirstSourceType;
  // type of second source term 'G' (second right hand side of differential equation 'div G')
  typedef Problem::PROBLEM_NAME::SecondSource SecondSourceType;
  // type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
  typedef Problem::PROBLEM_NAME::Diffusion DiffusionType;
  // type of mass (or reaction) term (i.e. 'm' or 'c')
  typedef Problem::PROBLEM_NAME::MassTerm MassTermType;
  // default type for any missing coefficient function (e.g. advection,...)
  typedef Problem::PROBLEM_NAME::DefaultDummyFunction DefaultDummyFunctionType;
  //!-----------------------------------------------------------------------------------------

  //! ---------  typedefs for the standard discrete function space (macroscopic) -------------

  // type of exact solution (in general unknown)
  typedef Problem::PROBLEM_NAME::ExactSolution ExactSolutionType;

  typedef FunctionSpaceType::DomainType DomainType;
  //! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
  typedef FunctionSpaceType::RangeType RangeType;
  //! defines the function space to which the numerical solution belongs to
  //! see dune/fem/lagrangebase.hh
  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >  // 1=POLORDER
      DiscreteFunctionSpaceType;
  typedef DiscreteFunctionSpaceType::DomainFieldType            TimeType;
  typedef DiscreteFunctionSpaceType::JacobianRangeType          JacobianRangeType;
  typedef GridType::Codim< 0 >::Entity                          EntityType;
  typedef GridType::Codim< 0 >::EntityPointer                   EntityPointerType;
  typedef GridType::Codim< 0 >::Geometry                        EntityGeometryType;
  typedef GridType::Codim< 1 >::Geometry                        FaceGeometryType;
  typedef DiscreteFunctionSpaceType::BaseFunctionSetType        BaseFunctionSetType;
  typedef Dune::CachingQuadrature< GridPartType, 0 >                  EntityQuadratureType;
  typedef Dune::CachingQuadrature< GridPartType, 1 >                  FaceQuadratureType;
  typedef DiscreteFunctionSpaceType::RangeFieldType             RangeFieldType;
  typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef DiscreteFunctionType::LocalFunctionType               LocalFunctionType;
  typedef DiscreteFunctionType::DofIteratorType                 DofIteratorType;

  //! --------------------- the standard matrix traits -------------------------------------
  struct MatrixTraits
  {
    typedef DiscreteFunctionSpaceType                          RowSpaceType;
    typedef DiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef Dune::LagrangeMatrixSetup< false >                       StencilType;
    typedef Dune::ParallelScalarProduct< DiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef Dune::LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };




  //!------------------------- for adaptive grid refinement ---------------------------------
  //! type of restrict-prolong operator
  typedef Dune::RestrictProlongDefault< DiscreteFunctionType > RestrictProlongOperatorType;
  //! type of the adaption manager
  typedef Dune::AdaptationManager< GridType, RestrictProlongOperatorType > AdaptationManagerType;
  //!---------------------------------------------------------------------------------------

  typedef std::vector< RangeType > RangeVector;
  typedef std::vector< RangeVector > RangeVectorVector;

  static const int assembler_order = 2* DiscreteFunctionSpaceType::polynomialOrder + 2;

  typedef Dune::SparseRowMatrixOperator< DiscreteFunctionType, DiscreteFunctionType, MatrixTraits > FEMMatrix;

};

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_COMMON_TRAITS_HH
