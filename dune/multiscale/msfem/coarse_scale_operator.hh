// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MSFEM_ELLIPTIC_DiscreteEllipticMSFEMOperator_HH
#define MSFEM_ELLIPTIC_DiscreteEllipticMSFEMOperator_HH


#include <assert.h>
#include <boost/noncopyable.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/multiscale/common/dirichletconstraints.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/fem/matrix_object.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <ostream>
#include <type_traits>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {
/**
 * \todo docme
 */
class CoarseScaleOperator : boost::noncopyable {
private:
  typedef CommonTraits::DiscreteFunctionType CoarseDiscreteFunction;
  typedef CommonTraits::DiscreteFunctionType FineDiscreteFunction;
  typedef MacroMicroGridSpecifier MacroMicroGridSpecifierType;

  typedef CommonTraits::DiffusionType DiffusionModel;

  typedef typename CoarseDiscreteFunction::DiscreteFunctionSpaceType CoarseDiscreteFunctionSpace;
  typedef typename FineDiscreteFunction::DiscreteFunctionSpaceType FineDiscreteFunctionSpace;

  typedef typename FineDiscreteFunctionSpace::FunctionSpaceType FunctionSpace;

  typedef typename FineDiscreteFunctionSpace::GridPartType FineGridPart;
  typedef typename FineDiscreteFunctionSpace::GridType FineGrid;

  typedef typename FineDiscreteFunctionSpace::RangeFieldType RangeFieldType;
  typedef typename FineDiscreteFunctionSpace::DomainType DomainType;
  typedef typename FineDiscreteFunctionSpace::RangeType RangeType;
  typedef typename FineDiscreteFunctionSpace::JacobianRangeType JacobianRangeType;

  typedef LocalProblemSolver LocalProblemSolverType;

  static const int dimension = FineGridPart::GridType::dimension;

  typedef typename FineDiscreteFunctionSpace::BasisFunctionSetType FineBaseFunctionSet;
  typedef typename FineDiscreteFunctionSpace::EntityType FineEntity;
  typedef typename FineEntity::EntityPointer FineEntityPointer;

  typedef typename CoarseDiscreteFunctionSpace::BasisFunctionSetType CoarseBaseFunctionSet;
  typedef typename CommonTraits::EntityType CoarseEntity;

  typedef typename CommonTraits::LinearOperatorType MatrixType;

public:
  CoarseScaleOperator(const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
                                LocalGridList& subgrid_list, const DiffusionModel& diffusion_op);

  void apply_inverse(const CoarseDiscreteFunction& b, CoarseDiscreteFunction& x);

private:
  const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace_;
  LocalGridList& subgrid_list_;
  const DiffusionModel& diffusion_operator_;
  const bool petrovGalerkin_;
  MatrixType global_matrix_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef MSFEM_ELLIPTIC_DiscreteElliptic_HH
