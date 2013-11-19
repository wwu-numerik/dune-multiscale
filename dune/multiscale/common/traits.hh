// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_COMMON_TRAITS_HH
#define DUNE_MULTISCALE_COMMON_TRAITS_HH


#include <dune/multiscale/common/la_backend.hh>
#include <dune/common/tuples.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/stuff/common/memory.hh>

namespace Dune {

template <class T>
struct GridPtr;

namespace Fem {
template <class T, class R>
class GridFunctionAdapter;
template <class T, class R>
class DataOutput;
template <class T, class R>
class DataWriter;
template <class T, class R>
class AdaptationManager;
} // namespace Fem

namespace Stuff {
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class GlobalFunction;
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class GlobalConstantFunction;
}
namespace Multiscale {
namespace Problem {

struct DiffusionBase;
struct LowerOrderBase;
class NeumannDataBase;
class DirichletDataBase;
class IModelProblemData;

} // namespace Problem

//! type construction for the HMM algorithm
struct CommonTraits {
  //! --------- typedefs for the macro grid and the corresponding discrete space -------------
  typedef Dune::GridSelector::GridType GridType;
  typedef GridType::Codim<0>::Entity EntityType;
  // Dune::InteriorBorder_Partition or Dune::All_Partition >?
  // see:
  // http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
  typedef Dune::Fem::AdaptiveLeafGridPart<GridType /*,Dune::All_Partition*/> GridPartType;
  typedef Dune::GridPtr<GridType> GridPointerType;
  typedef Dune::Fem::FunctionSpace<double, double, GridType::dimension, 1> FunctionSpaceType;
  //!-----------------------------------------------------------------------------------------

  typedef Dune::Stuff::GlobalFunction<EntityType, FunctionSpaceType::DomainFieldType, FunctionSpaceType::dimDomain,
                                         FunctionSpaceType::RangeFieldType, FunctionSpaceType::dimRange> FunctionBaseType;

  typedef Dune::Stuff::GlobalConstantFunction<EntityType, FunctionSpaceType::DomainFieldType, FunctionSpaceType::dimDomain,
                                        FunctionSpaceType::RangeFieldType,
                                        FunctionSpaceType::dimRange> ConstantFunctionBaseType;
  //! --------- typedefs for the coefficient and data functions ------------------------------
  typedef Problem::IModelProblemData ModelProblemDataType;
  // type of first source term (right hand side of differential equation or type of 'f')
  typedef FunctionBaseType FirstSourceType;
  // type of second source term 'G' (second right hand side of differential equation 'div G')
  typedef FunctionBaseType SecondSourceType;
  // type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
  typedef Problem::DiffusionBase DiffusionType;
  // type of (possibly non-linear) lower order term F( x, u(x), grad u(x) )
  typedef Problem::LowerOrderBase LowerOrderTermType;
  // type of inhomogeneous Dirichlet boundary condition
  typedef FunctionBaseType DirichletBCType;
  // type of inhomogeneous Neumann boundary condition
  typedef FunctionBaseType NeumannBCType;
  // type of dirichlet data
  typedef Problem::DirichletDataBase DirichletDataType;
  // type of neumann data
  typedef Problem::NeumannDataBase NeumannDataType;
  // type of mass (or reaction) term (i.e. 'm' or 'c')
  typedef FunctionBaseType MassTermType;
  // default type for any missing coefficient function (e.g. advection,...)
  typedef FunctionBaseType DefaultDummyFunctionType;
  //!-----------------------------------------------------------------------------------------

  //! ---------  typedefs for the standard discrete function space (macroscopic) -------------

  // type of exact solution (in general unknown)
  typedef FunctionBaseType ExactSolutionType;

  typedef FunctionSpaceType::DomainType DomainType;
  //! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
  typedef FunctionSpaceType::RangeType RangeType;
  //! defines the function space to which the numerical solution belongs to
  //! see dune/fem/lagrangebase.hh
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType, st_lagrangespace_order>
      DiscreteFunctionSpaceType;
  typedef DiscreteFunctionSpaceType::DomainFieldType TimeType;
  typedef DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef GridPartType::IntersectionType FaceType;
  typedef GridType::Codim<0>::EntityPointer EntityPointerType;
  typedef GridType::Codim<0>::Geometry EntityGeometryType;
  typedef GridType::Codim<1>::Geometry FaceGeometryType;
  //!TODO carry the rename over to the type def'ed name
  typedef DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  typedef BackendChooser<DiscreteFunctionSpaceType>::DiscreteFunctionType DiscreteFunctionType;
  typedef std::shared_ptr<DiscreteFunctionType> DiscreteFunction_ptr;
  typedef BackendChooser<DiscreteFunctionSpaceType>::LinearOperatorType LinearOperatorType;

  //!------------------------- for adaptive grid refinement ---------------------------------
  //! type of restrict-prolong operator
  typedef Dune::Fem::RestrictProlongDefault<DiscreteFunctionType> RestrictProlongOperatorType;
  //! type of the adaption manager
  typedef Dune::Fem::AdaptationManager<GridType, RestrictProlongOperatorType> AdaptationManagerType;
  //!---------------------------------------------------------------------------------------

  typedef std::vector<RangeType> RangeVector;
  typedef std::vector<RangeVector> RangeVectorVector;

  static const int assembler_order = 2 * DiscreteFunctionSpaceType::polynomialOrder + 2;
};

template< class T = CommonTraits::DiscreteFunctionType>
std::shared_ptr<T> make_df_ptr(const std::string name, const typename T::DiscreteFunctionSpaceType& space)
{
  return std::make_shared<T>(name, space);
}

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_COMMON_TRAITS_HH
