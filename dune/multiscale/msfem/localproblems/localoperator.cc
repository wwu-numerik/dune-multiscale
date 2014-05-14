#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "localoperator.hh"

#include <assert.h>
#include <boost/assert.hpp>
#include <dune/common/exceptions.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/matrix_object.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/common/fmatrix.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {

  
  std::vector<MsFEMTraits::CoarseBaseFunctionSetType::JacobianRangeType> LocalProblemOperator::coarseBaseJacs_;
  std::vector<CommonTraits::BaseFunctionSetType::JacobianRangeType> LocalProblemOperator::dirichletJacs_;
  bool LocalProblemOperator::cached_;
  
LocalProblemOperator::LocalProblemOperator(const CoarseSpaceType& coarse_space,
                                           const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace,
                                           const DiffusionOperatorType& diffusion_op)
  : subDiscreteFunctionSpace_(subDiscreteFunctionSpace)
  , diffusion_operator_(diffusion_op)
  , coarse_space_(coarse_space)
  , system_matrix_("Local Problem System Matrix", subDiscreteFunctionSpace, subDiscreteFunctionSpace)
{
  assemble_matrix();
}

void LocalProblemOperator::assemble_matrix()
    // x_T is the barycenter of the macro grid element T
{
  system_matrix_.reserve(DSFe::diagonalAndNeighborStencil(system_matrix_));
  system_matrix_.clear();

  // local grid basis functions:
  std::vector<RangeType> phi(subDiscreteFunctionSpace_.blockMapper().maxNumDofs());

  // gradient of micro scale base function:
  std::vector<typename BaseFunctionSetType::JacobianRangeType> gradient_phi(
      subDiscreteFunctionSpace_.blockMapper().maxNumDofs());
  typename BaseFunctionSetType::JacobianRangeType diffusion_in_gradient_phi;

  for (const auto& sub_grid_entity : subDiscreteFunctionSpace_) {
    const auto& sub_grid_geometry = sub_grid_entity.geometry();

    DSFe::LocalMatrixProxy<LocalProblemSolver::LinearOperatorType> local_matrix(system_matrix_, sub_grid_entity,
                                                                                sub_grid_entity);

    const auto& baseSet = local_matrix.domainBasisFunctionSet();
    const auto numBaseFunctions = baseSet.size();

    // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to
    // use a higher order quadrature:
    const auto quadrature = DSFe::make_quadrature(sub_grid_entity, subDiscreteFunctionSpace_);

    const auto numQuadraturePoints = quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      // local (barycentric) coordinates (with respect to local grid entity)
      const auto& local_point = quadrature.point(quadraturePoint);
      const auto global_point = sub_grid_geometry.global(local_point);

      const double weight = quadrature.weight(quadraturePoint) * sub_grid_geometry.integrationElement(local_point);

      baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi);

      for (unsigned int i = 0; i < numBaseFunctions; ++i) {
        // A( x, \nabla \phi(x) )
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi[i], diffusion_in_gradient_phi);
        for (unsigned int j = 0; j < numBaseFunctions; ++j) {
          // stiffness contribution
          local_matrix.add(j, i, weight * (diffusion_in_gradient_phi[0] * gradient_phi[j][0]));
          // mass contribution (just for stabilization!)
          // local_matrix.add( j, i, 0.00000001 * weight * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
        }
      }
    }
  }
  DSG::BoundaryInfos::AllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> boundaryInfo;
  DirichletConstraints<MsFEMTraits::LocalGridDiscreteFunctionType> constraints(boundaryInfo, subDiscreteFunctionSpace_);
  constraints.applyToOperator(system_matrix_);
  system_matrix_.communicate();
} // assemble_matrix

long LocalProblemOperator::getNumQuadPoints(const MsFEMTraits::LocalGridDiscreteFunctionSpaceType& discreteFunctionSpace) const {
  const auto quadOrder = (2*discreteFunctionSpace.order()+2);
  int m = 0;
  while ((2*m-1 < quadOrder) && m < quadOrder)
    ++m;
  return std::pow(m, CommonTraits::GridType::dimension);
}
  
void LocalProblemOperator::assemble_all_local_rhs(const CoarseEntityType& coarseEntity,
                                                  MsFEMTraits::LocalSolutionVectorType& allLocalRHS) const {
  BOOST_ASSERT_MSG(allLocalRHS.size() > 0, "You need to preallocate the necessary space outside this function!");

  //! @todo correct the error message below (+1 for simplecial, +2 for arbitrary), as there's no finespace any longer
  //  BOOST_ASSERT_MSG(
  //      (DSG::is_simplex_grid(coarse_space_) && allLocalRHS.size() == GridType::dimension + 1) ||
  //          (!(DSG::is_simplex_grid(coarse_space_)) &&
  //           static_cast<long long>(allLocalRHS.size()) ==
  //               static_cast<long long>(specifier.fineSpace().mapper().maxNumDofs() + 2)),
  //      "You need to allocate storage space for the correctors for all unit vector/all coarse basis functions"
  //      " and the dirichlet- and neuman corrector");

  // build unit vectors (needed for cases where rhs is assembled for unit vectors instead of coarse
  // base functions)
  JacobianRangeType unitVectors[dimension];
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      if (i == j) {
        unitVectors[i][0][j] = 1.0;
      } else {
        unitVectors[i][0][j] = 0.0;
      }
    }
  }

  // get dirichlet and neumann data
  const auto& neumannData = *Dune::Multiscale::Problem::getNeumannData();
  const auto& discreteFunctionSpace = allLocalRHS[0]->space();

  //! @todo we should use the dirichlet constraints here somehow
  LocalGridDiscreteFunctionType dirichletExtension("dirichletExtension", discreteFunctionSpace);
  dirichletExtension.clear();
  CommonTraits::DiscreteFunctionType dirichletExtensionCoarse("Dirichlet Extension Coarse", coarse_space_);
  dirichletExtensionCoarse.clear();
  //! @todo is this needed or could it be replaced by a method from dirichletconstraints.hh?
  project_dirichlet_values(dirichletExtensionCoarse);
  Dune::Stuff::HeterogenousProjection<> projection;
  projection.project(dirichletExtensionCoarse, dirichletExtension);

  // set entries to zero:
  for (auto& rhs : allLocalRHS)
    rhs->clear();

  // get the base function set of the coarse space for the given coarse entity
  const auto& coarseBaseSet = coarse_space_.basisFunctionSet(coarseEntity);
  std::vector<CoarseBaseFunctionSetType::JacobianRangeType> coarseBaseFuncJacs(coarseBaseSet.size());

  auto numBasisFunctions = discreteFunctionSpace.blockMapper().maxNumDofs();
  std::vector<RangeType> phi(numBasisFunctions);

  const bool is_simplex_grid = DSG::is_simplex_grid(coarse_space_);
  const auto numBoundaryCorrectors = is_simplex_grid ? 1u : 2u;
  const auto numInnerCorrectors = allLocalRHS.size() - numBoundaryCorrectors;

  int coarseJacCacheCounter = 0;
  int dirichletJacCacheCounter = 0;
  JacobianRangeType dirichletJac(0.0);
  CoarseBaseFunctionSetType::JacobianRangeType coarseBaseJac;
  // reserve memory for the jacobian caches
  if (!cached_) {
    // compute the number of quadrature points resulting from a standard quadrature on the current space
    const auto& numQuadPoints = getNumQuadPoints(discreteFunctionSpace);
    const auto& gridSize = discreteFunctionSpace.grid_view().grid().size(0);
    coarseBaseJacs_.reserve(gridSize * numQuadPoints * numInnerCorrectors);
    dirichletJacs_.reserve(gridSize * numQuadPoints);
  }
  
  for (auto& localGridCell : discreteFunctionSpace) {
    const auto& geometry = localGridCell.geometry();
    const bool hasBoundaryIntersection = localGridCell.hasBoundaryIntersections();
    auto dirichletLF = dirichletExtension.localFunction(localGridCell);

    const auto quadrature = DSFe::make_quadrature(localGridCell, discreteFunctionSpace);
    const auto numQuadraturePoints = quadrature.nop();
    std::vector<std::vector<JacobianRangeType>> gradient_phi(numQuadraturePoints,
                                                             std::vector<JacobianRangeType>(numBasisFunctions));

    const auto& baseSet = discreteFunctionSpace.basisFunctionSet(localGridCell);
    const auto numBaseFunctions = baseSet.size();
    for (const auto& qp : DSC::valueRange(numQuadraturePoints))
      baseSet.jacobianAll(quadrature[qp], gradient_phi[qp]);
    
    for (std::size_t coarseBaseFunc = 0; coarseBaseFunc < allLocalRHS.size(); ++coarseBaseFunc) {
      auto rhsLocalFunction = allLocalRHS[coarseBaseFunc]->localFunction(localGridCell);

      // correctors with index < numInnerCorrectors are for the basis functions, corrector at
      // position numInnerCorrectors is for the neumann values, corrector at position numInnerCorrectors+1
      // for the dirichlet values.
      if (coarseBaseFunc < numInnerCorrectors || coarseBaseFunc == numInnerCorrectors + 1) {
        for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
          const auto& local_point = quadrature.point(quadraturePoint);
          // global point in the subgrid
          const auto global_point = geometry.global(local_point);

          const double weight = quadrature.weight(quadraturePoint) * geometry.integrationElement(local_point);

          JacobianRangeType diffusion(0.0);
          if (coarseBaseFunc < numInnerCorrectors) {
            if (is_simplex_grid)
              diffusion_operator_.diffusiveFlux(global_point, unitVectors[coarseBaseFunc], diffusion);
            else {
              const DomainType quadInCoarseLocal = coarseEntity.geometry().local(global_point);
              if (!cached_) {
                coarseBaseSet.jacobianAll(quadInCoarseLocal, coarseBaseFuncJacs);
                coarseBaseJac = coarseBaseFuncJacs[coarseBaseFunc];
                coarseBaseJacs_.push_back(coarseBaseJac);
              } else
                coarseBaseJac = coarseBaseJacs_.at(coarseJacCacheCounter++);
              diffusion_operator_.diffusiveFlux(global_point, coarseBaseJac, diffusion);
            }
          } else {
            if (!cached_) {
              dirichletLF.jacobian(local_point, dirichletJac);
              dirichletJacs_.push_back(dirichletJac);
            }
            else
              dirichletJac = dirichletJacs_.at(dirichletJacCacheCounter++);
            diffusion_operator_.diffusiveFlux(global_point, dirichletJac, diffusion);
          }
          
          for (unsigned int i = 0; i < numBaseFunctions; ++i) {
            rhsLocalFunction[i] -= weight * (diffusion[0] * gradient_phi[quadraturePoint][i][0]);
          }
        }
      }

      // boundary integrals
      if (coarseBaseFunc == numInnerCorrectors && hasBoundaryIntersection) {
        const auto intEnd = discreteFunctionSpace.gridPart().iend(localGridCell);
        for (auto iIt = discreteFunctionSpace.gridPart().ibegin(localGridCell); iIt != intEnd; ++iIt) {
          const auto& intersection = *iIt;
          if (DMP::is_neumann(intersection)) {
            const auto orderOfIntegrand =
                (CommonTraits::polynomial_order - 1) + 2 * (CommonTraits::polynomial_order + 1);
            const auto quadOrder = std::ceil((orderOfIntegrand + 1) / 2);
            const auto faceQuad = DSFe::make_quadrature(intersection, discreteFunctionSpace, quadOrder);
            RangeType neumannValue(0.0);
            const auto numQuadPoints = faceQuad.nop();
            // loop over all quadrature points
            for (unsigned int iqP = 0; iqP < numQuadPoints; ++iqP) {
              // get local coordinate of quadrature point
              const auto& xLocal = faceQuad.localPoint(iqP);
              const auto& faceGeometry = intersection.geometry();

              // the following does not work because subgrid does not implement geometryInInside()
              // const auto& insideGeometry    = intersection.geometryInInside();
              // const typename FaceQuadratureType::CoordinateType& xInInside = insideGeometry.global(xLocal);
              // therefore, we have to do stupid things:
              const auto& xGlobal = faceGeometry.global(xLocal);
              auto insidePtr = intersection.inside();
              const auto& insideEntity = *insidePtr;
              const auto& xInInside = insideEntity.geometry().local(xGlobal);
              const double factor = faceGeometry.integrationElement(xLocal) * faceQuad.weight(iqP);

              neumannData.evaluate(xGlobal, neumannValue);
              baseSet.evaluateAll(xInInside, phi);
              for (unsigned int i = 0; i < numBaseFunctions; ++i) {
                rhsLocalFunction[i] -= factor * (neumannValue * phi[i]);
              }
            }
          }
        }
      }
    }
  }
  
  // dirichlet jacobians and coarse base func jacobians were cached
  cached_ = true;

  DSG::BoundaryInfos::AllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> boundaryInfo;
  DirichletConstraints<MsFEMTraits::LocalGridDiscreteFunctionType> constraints(boundaryInfo, discreteFunctionSpace);
  for (auto& rhsIt : allLocalRHS) {
    constraints.setValue(0.0, *rhsIt);
  }
  return;
}

void LocalProblemOperator::project_dirichlet_values(CommonTraits::DiscreteFunctionType& function) const {
  /*  // make sure, we are on a hexahedral element
    BOOST_ASSERT_MSG(function.space().grid_view()->grid().leafIndexSet().geomTypes(0).size()==1 &&
           function.space().grid_view()->grid().leafIndexSet().geomTypes(0)[0].isCube(),
           "This method only works for hexahedral elements at the moment!");*/

  const auto& gridPart = function.space().gridPart();
  const auto& dirichletData = *Multiscale::Problem::getDirichletData();
  //  Fem::GridFunctionAdapter<Multiscale::Problem::DirichletDataBase,
  //  typename LocalGridDiscreteFunctionType::LocalGridDiscreteFunctionSpaceType::GridPartType> gf("dirichlet",
  // dirichletData ,
  // gridPart);
  RangeType dirichletVal(0.0);
  for (const auto& localCell : function.space()) {
    if (localCell.hasBoundaryIntersections())
      for (const auto& intersection : DSC::intersectionRange(gridPart, localCell)) {
        if (DMP::is_dirichlet(intersection)) {
          auto funcLocal = function.localFunction(localCell);
          const auto& lagrangePointSet = function.space().lagrangePointSet(localCell);
          const auto faceNumber = intersection.indexInInside();
          for (auto lp : DSC::lagrangePointSetRange<1>(function.space(), localCell, faceNumber)) {
            auto lagrangePoint = lagrangePointSet.point(lp);
            auto lagrangePointGlobal = localCell.geometry().global(lagrangePoint);
            dirichletData.evaluate(lagrangePointGlobal, dirichletVal);
            funcLocal[lp] = dirichletVal;
          }
        }
      }
  }
  return;
}

void LocalProblemOperator::apply_inverse(const MsFEMTraits::LocalGridDiscreteFunctionType &current_rhs, MsFEMTraits::LocalGridDiscreteFunctionType &current_solution)
{
  if (!current_rhs.dofsValid())
    DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");

  const auto solver =
      Dune::Multiscale::Problem::getModelData()->symmetricDiffusion() ? std::string("cg") : std::string("bcgs");
  typedef BackendChooser<LocalGridDiscreteFunctionSpaceType>::InverseOperatorType LocalInverseOperatorType;
  const auto localProblemSolver = DSC::make_unique<LocalInverseOperatorType>(system_matrix_, 1e-8, 1e-8, 20000,
                                               DSC_CONFIG_GET("msfem.localproblemsolver_verbose", false), solver,
                                               DSC_CONFIG_GET("preconditioner_type", std::string("sor")), 1);
  localProblemSolver->apply(current_rhs, current_solution);

  if (!current_solution.dofsValid())
    DUNE_THROW(Dune::InvalidStateException, "Current solution of the local msfem problem invalid!");
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
