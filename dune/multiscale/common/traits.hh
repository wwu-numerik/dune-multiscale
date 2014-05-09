// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_COMMON_TRAITS_HH
#define DUNE_MULTISCALE_COMMON_TRAITS_HH

#include <dune/multiscale/common/la_backend.hh>
#include <dune/common/tuples.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/functions/global.hh>

namespace Dune {

template <class T>
struct GridPtr;

namespace PDELab {

class NoConstraints;
//class ConformingDirichletConstraints;

} //namespace PDELab

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

namespace GDT {
template <class T, class R>
class DiscreteFunction;
template <class T, class R>
class ConstDiscreteFunction;
}
namespace Multiscale {
namespace Problem {

struct DiffusionBase;
struct LowerOrderBase;
class NeumannDataBase;
class DirichletDataBase;
class IModelProblemData;

} // namespace Problem

//! Common Type, duh
struct CommonTraits {

  typedef Dune::GridSelector::GridType GridType;
  static constexpr int dimRange = 1;
  static constexpr int dimDomain = GridType::dimension;
  static constexpr int world_dim = dimDomain
                                     ;
//  typedef Dune::SGrid<world_dim, world_dim> GridType;
//  typedef Dune::YaspGrid<world_dim> GridType;
  typedef GridType::Codim<0>::Entity EntityType;
  typedef Dune::Fem::AdaptiveLeafGridPart<GridType> GridPartType;
  typedef Dune::GridPtr<GridType> GridPointerType;
  typedef double FieldType;
  typedef Dune::Fem::FunctionSpace<FieldType, FieldType, dimDomain, dimRange> FunctionSpaceType;

  typedef Stuff::GlobalFunctionInterface<EntityType, FieldType, dimDomain, FieldType, dimRange> FunctionBaseType;
  typedef Stuff::GlobalFunctionInterface<EntityType, FieldType, dimDomain, FieldType, dimDomain, dimDomain> DiffusionFunctionBaseType;
  typedef Stuff::Functions::Constant< EntityType, FieldType, dimDomain, FieldType, dimRange > ConstantFunctionBaseType;

  typedef Problem::IModelProblemData ModelProblemDataType;
  //! type of first source term (right hand side of differential equation or type of 'f')
  typedef FunctionBaseType SourceType;

  //! type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
  typedef Problem::DiffusionBase DiffusionType;
  //! type of (possibly non-linear) lower order term F( x, u(x), grad u(x) )
  typedef Problem::LowerOrderBase LowerOrderTermType;
  //! type of inhomogeneous Dirichlet boundary condition
  typedef FunctionBaseType DirichletBCType;
  //! type of inhomogeneous Neumann boundary condition
  typedef FunctionBaseType NeumannBCType;
  //! type of dirichlet data
  typedef Problem::DirichletDataBase DirichletDataType;
  //! type of neumann data
  typedef Problem::NeumannDataBase NeumannDataType;

  //! type of exact solution (in general unknown)
  typedef FunctionBaseType ExactSolutionType;
  static constexpr unsigned int exact_solution_space_order = 3 * st_lagrangespace_order;

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
  typedef FieldType RangeFieldType;
  typedef FieldType DomainFieldType;

  typedef BackendChooser<DiscreteFunctionSpaceType>::DiscreteFunctionType DiscreteFunctionType;
  typedef std::shared_ptr<DiscreteFunctionType> DiscreteFunction_ptr;
  typedef BackendChooser<DiscreteFunctionSpaceType>::LinearOperatorType LinearOperatorType;

  typedef std::vector<RangeType> RangeVector;
  typedef std::vector<RangeVector> RangeVectorVector;

  static constexpr int polynomial_order = DiscreteFunctionSpaceType::polynomialOrder;
  static constexpr int quadrature_order = 2 * polynomial_order + 2;

  //GDT::ChooseSpaceBackend::pdelab UsePdelab;

  typedef GDT::Spaces::ContinuousLagrangeProvider< GridType, DSG::ChooseLayer::leaf,
                                                   GDT::ChooseSpaceBackend::fem,
                                                   st_lagrangespace_order, FieldType, dimRange > GdtSpaceProviderType;

  static constexpr auto st_gdt_grid_level = 0;
  typedef GdtSpaceProviderType::Type GdtSpaceType;
  typedef GdtSpaceType::GridViewType GridViewType;
  typedef BackendChooser<DiscreteFunctionSpaceType>::GdtMatrixType GdtMatrixType;
  typedef BackendChooser<DiscreteFunctionSpaceType>::GdtVectorType GdtVectorType;
  typedef GDT::DiscreteFunction< GdtSpaceType, GdtVectorType >      GdtDiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction< GdtSpaceType, GdtVectorType > GdtConstDiscreteFunctionType;
  typedef ConstantFunctionBaseType GdtConstantFunctionType;
};

template <class T = CommonTraits::DiscreteFunctionType>
std::shared_ptr<T> make_df_ptr(const std::string name, const typename T::DiscreteFunctionSpaceType& space) {
  return std::make_shared<T>(name, space);
  //  return DSC::make_unique<T>(name, space);
}

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_COMMON_TRAITS_HH
