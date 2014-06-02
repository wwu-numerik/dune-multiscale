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
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/functionals/l2.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {
  
LocalProblemOperator::LocalProblemOperator(const CoarseSpaceType& coarse_space,
                                           const LocalGridDiscreteFunctionSpaceType& space,
                                           const DiffusionOperatorType& diffusion_op)
  : localSpace_(space)
  , diffusion_operator_(diffusion_op)
  , coarse_space_(coarse_space)
  , system_matrix_(space.mapper().size(), space.mapper().size(),
                   EllipticOperatorType::pattern(space))
  , system_assembler_(localSpace_)
  , elliptic_operator_(diffusion_operator_, system_matrix_, localSpace_)
  , constraints_(Problem::getModelData()->subBoundaryInfo(), space.mapper().maxNumDofs(), space.mapper().maxNumDofs())

{
  assemble_matrix();
}

void LocalProblemOperator::assemble_matrix()
    // x_T is the barycenter of the macro grid element T
{
  system_assembler_.add(elliptic_operator_);


} // assemble_matrix

void LocalProblemOperator::assemble_all_local_rhs(const CoarseEntityType& coarseEntity,
                                                  MsFEMTraits::LocalSolutionVectorType& allLocalRHS) {
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
  constexpr auto dimension = CommonTraits::GridType::dimension;
  CommonTraits::JacobianRangeType unitVectors[dimension];
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      if (i == j) {
        unitVectors[i][0][j] = 1.0;
      } else {
        unitVectors[i][0][j] = 0.0;
      }
    }
  }

  LocalGridDiscreteFunctionType dirichletExtension(localSpace_, "dirichletExtension");
  CommonTraits::DiscreteFunctionType dirichletExtensionCoarse(coarse_space_, "Dirichlet Extension Coarse");

  GDT::Operators::DirichletProjectionLocalizable< CommonTraits::GridViewType, CommonTraits::DirichletDataType,
      CommonTraits::DiscreteFunctionType >
      coarse_dirichlet_projection_operator(*(coarse_space_.grid_view()),
                                    DMP::getModelData().boundaryInfo(),
                                    *DMP::getDirichletData(),
                                    dirichletExtensionCoarse);
  system_assembler_.add(coarse_dirichlet_projection_operator);
  system_assembler_.assemble();
  Dune::Stuff::HeterogenousProjection<> projection;
  projection.project(dirichletExtensionCoarse, dirichletExtension);

  const bool is_simplex_grid = DSG::is_simplex_grid(coarse_space_);
  const auto numBoundaryCorrectors = is_simplex_grid ? 1u : 2u;
  const auto numInnerCorrectors = allLocalRHS.size() - numBoundaryCorrectors;
  //!*********** anfang neu gdt

  std::size_t coarseBaseFunc = 0;
  for (; coarseBaseFunc < numInnerCorrectors; ++coarseBaseFunc)
  {
    if (is_simplex_grid) {
//      diffusionsauswertung in unitVectors[coarseBaseFunc]
    } else {
//      diffusionsauswertung in
//      const DomainType quadInCoarseLocal = coarseEntity.geometry().local(QuadraturPunkt);
//      coarseBaseSet.jacobianAll(quadInCoarseLocal, coarseBaseFuncJacs);
    }
//    baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
//    for (unsigned int i = 0; i < numBaseFunctions; ++i) {
//      rhsLocalFunction[i] -= weight * (diffusions_auswertung * gradient_phi[i][0]);
//    }
  }

  coarseBaseFunc++; // coarseBaseFunc == numInnerCorrectors
  //neumann correktor
  GDT::Functionals::L2Face< CommonTraits::NeumannDataType, CommonTraits::GdtVectorType, MsFEMTraits::LocalSpaceType >
      neumann_functional(*Dune::Multiscale::Problem::getNeumannData(),
                         allLocalRHS[coarseBaseFunc], localSpace_);
  system_assembler_.add(neumann_functional);


  coarseBaseFunc++;// coarseBaseFunc == 1 + numInnerCorrectors
  //dirichlet correktor
  {
//    const auto dirichletLF = dirichletExtension.local_function(entity);
//    dirichletLF.jacobian(local_point, dirichletJac);
//    diffusion_operator_.diffusiveFlux(global_point, dirichletJac, diffusion);
//    for (unsigned int i = 0; i < numBaseFunctions; ++i) {
//      rhsLocalFunction[i] -= weight * (diffusion[0] * gradient_phi[i][0]);
//    }
  }

  //dirichlet-0 for all rhs
  
  LocalGridDiscreteFunctionType dirichlet_projection(localSpace_);
  GDT::Operators::DirichletProjectionLocalizable< CommonTraits::GridViewType, CommonTraits::DirichletDataType, CommonTraits::DiscreteFunctionType >
      dirichlet_projection_operator(*(space.grid_view()),
                                    allLocalDirichletInfo_,
                                    *dirichletZero_,
                                    dirichlet_projection);
  system_assembler_.add(dirichlet_projection_operator,
                       new GDT::ApplyOn::BoundaryEntities< GridViewType >());


  system_assembler_.add(constraints_, system_matrix_, new GDT::ApplyOn::BoundaryEntities< MsFEMTraits::LocalGridViewType>());
  for (auto& rhs : allLocalRHS )
    system_assembler_.add(constraints_, rhs, new GDT::ApplyOn::BoundaryEntities< MsFEMTraits::GridViewType >());
  system_assembler_.assemble();

  //!*********** ende neu gdt

#if 0 // alter dune-fem code
  // get the base function set of the coarse space for the given coarse entity
  const auto& coarseBaseSet = coarse_space_.basisFunctionSet(coarseEntity);
  std::vector<CoarseBaseFunctionSetType::JacobianRangeType> coarseBaseFuncJacs(coarseBaseSet.size());

  // gradient of micro scale base function:
  std::vector<JacobianRangeType> gradient_phi(localSpace.blockMapper().maxNumDofs());
  std::vector<RangeType> phi(localSpace.blockMapper().maxNumDofs());


  for (auto& localGridCell : localSpace) {
    const auto& geometry = localGridCell.geometry();
    const bool hasBoundaryIntersection = localGridCell.hasBoundaryIntersections();
    auto dirichletLF = dirichletExtension.localFunction(localGridCell);
    JacobianRangeType dirichletJac(0.0);

    for (std::size_t coarseBaseFunc = 0; coarseBaseFunc < allLocalRHS.size(); ++coarseBaseFunc) {
      auto rhsLocalFunction = allLocalRHS[coarseBaseFunc]->localFunction(localGridCell);

      const auto& baseSet = rhsLocalFunction.basisFunctionSet();
      const auto numBaseFunctions = baseSet.size();

      // correctors with index < numInnerCorrectors are for the basis functions, corrector at
      // position numInnerCorrectors is for the neumann values, corrector at position numInnerCorrectors+1
      // for the dirichlet values.
      if (coarseBaseFunc < numInnerCorrectors || coarseBaseFunc == numInnerCorrectors + 1) {
        const auto quadrature = DSFe::make_quadrature(localGridCell, localSpace);
        const auto numQuadraturePoints = quadrature.nop();
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
              coarseBaseSet.jacobianAll(quadInCoarseLocal, coarseBaseFuncJacs);
              diffusion_operator_.diffusiveFlux(global_point, coarseBaseFuncJacs[coarseBaseFunc], diffusion);
            }
          } else {
            dirichletLF.jacobian(local_point, dirichletJac);
            diffusion_operator_.diffusiveFlux(global_point, dirichletJac, diffusion);
          }
          baseSet.jacobianAll(quadrature[quadraturePoint], gradient_phi);
          for (unsigned int i = 0; i < numBaseFunctions; ++i) {
            rhsLocalFunction[i] -= weight * (diffusion[0] * gradient_phi[i][0]);
          }
        }
      }

      // boundary integrals
      if (coarseBaseFunc == numInnerCorrectors && hasBoundaryIntersection) {
        const auto intEnd = localSpace.gridPart().iend(localGridCell);
        for (auto iIt = localSpace.gridPart().ibegin(localGridCell); iIt != intEnd; ++iIt) {
          const auto& intersection = *iIt;
          if (DMP::is_neumann(intersection)) {
            const auto orderOfIntegrand =
                (CommonTraits::polynomial_order - 1) + 2 * (CommonTraits::polynomial_order + 1);
            const auto quadOrder = std::ceil((orderOfIntegrand + 1) / 2);
            const auto faceQuad = DSFe::make_quadrature(intersection, localSpace, quadOrder);
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
#endif

}


void LocalProblemOperator::apply_inverse(const MsFEMTraits::LocalGridDiscreteFunctionType &current_rhs, MsFEMTraits::LocalGridDiscreteFunctionType &current_solution)
{
  if (!current_rhs.dofs_valid())
    DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");

  const auto solver =
      Dune::Multiscale::Problem::getModelData()->symmetricDiffusion() ? std::string("cg") : std::string("bcgs");
  typedef BackendChooser<LocalGridDiscreteFunctionSpaceType>::InverseOperatorType LocalInverseOperatorType;
  const auto localProblemSolver = DSC::make_unique<LocalInverseOperatorType>(system_matrix_, 1e-8, 1e-8, 20000,
                                               DSC_CONFIG_GET("msfem.localproblemsolver_verbose", false), solver,
                                               DSC_CONFIG_GET("preconditioner_type", std::string("sor")), 1);
  localProblemSolver->apply(current_rhs, current_solution);

  if (!current_solution.dofs_valid())
    DUNE_THROW(Dune::InvalidStateException, "Current solution of the local msfem problem invalid!");
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
