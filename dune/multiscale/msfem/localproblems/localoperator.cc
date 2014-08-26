#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "localoperator.hh"

#include <assert.h>
#include <boost/assert.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
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
  , local_diffusion_operator_(diffusion_operator_)
  , coarse_space_(coarse_space)
  , system_matrix_(space.mapper().size(), space.mapper().size(),
                   EllipticOperatorType::pattern(space))
  , system_assembler_(localSpace_)
  , elliptic_operator_(local_diffusion_operator_, system_matrix_, localSpace_)
  , dirichletConstraints_(Problem::getModelData()->subBoundaryInfo(), space.mapper().maxNumDofs(), space.mapper().maxNumDofs())
{
  system_assembler_.add(elliptic_operator_);
}

void LocalProblemOperator::assemble_all_local_rhs(const CoarseEntityType& coarseEntity,
                                                  MsFEMTraits::LocalSolutionVectorType& allLocalRHS) {
  BOOST_ASSERT_MSG(allLocalRHS.size() > 0, "You need to preallocate the necessary space outside this function!");

  //! @todo correct the error message below (+1 for simplicial, +2 for arbitrary), as there's no finespace any longer
  //  BOOST_ASSERT_MSG(
  //      (DSG::is_simplex_grid(coarse_space_) && allLocalRHS.size() == GridType::dimension + 1) ||
  //          (!(DSG::is_simplex_grid(coarse_space_)) &&
  //           static_cast<long long>(allLocalRHS.size()) ==
  //               static_cast<long long>(specifier.fineSpace().mapper().maxNumDofs() + 2)),
  //      "You need to allocate storage space for the correctors for all unit vector/all coarse basis functions"
  //      " and the dirichlet- and neuman corrector");

  LocalGridDiscreteFunctionType dirichletExtensionLocal(localSpace_, "dirichletExtension");
  CommonTraits::DiscreteFunctionType dirichletExtensionCoarse(coarse_space_, "Dirichlet Extension Coarse");

  GDT::SystemAssembler<CommonTraits::DiscreteFunctionSpaceType> global_system_assembler_(coarse_space_);
  GDT::Operators::DirichletProjectionLocalizable< CommonTraits::GridViewType, Problem::DirichletDataBase,
      CommonTraits::DiscreteFunctionType >
      coarse_dirichlet_projection_operator(*(coarse_space_.grid_view()),
                                    DMP::getModelData()->boundaryInfo(),
                                    *DMP::getDirichletData(),
                                    dirichletExtensionCoarse);
  global_system_assembler_.add(coarse_dirichlet_projection_operator, new GDT::ApplyOn::BoundaryEntities<CommonTraits::GridViewType>());
  global_system_assembler_.assemble();
  GDT::Operators::LagrangeProlongation<MsFEMTraits::LocalGridViewType> projection(*localSpace_.grid_view());
  projection.apply(dirichletExtensionCoarse, dirichletExtensionLocal);

  const bool is_simplex_grid = DSG::is_simplex_grid(coarse_space_);
  const auto numBoundaryCorrectors = is_simplex_grid ? 1u : 2u;
  const auto numInnerCorrectors = allLocalRHS.size() - numBoundaryCorrectors;
  //!*********** anfang neu gdt

  if (is_simplex_grid)
      DUNE_THROW(NotImplemented, "special treatment for simplicial grids missing");

  typedef GDT::Functionals::L2Volume< Problem::LocalDiffusionType, CommonTraits::GdtVectorType,
      MsFEMTraits::LocalSpaceType, MsFEMTraits::LocalGridViewType,
      CoarseBasisProduct > RhsFunctionalType;
  std::vector<std::unique_ptr<RhsFunctionalType>> rhs_functionals(numInnerCorrectors);
  std::size_t coarseBaseFunc = 0;
  const auto coarseBaseFunctionSet = coarse_space_.base_function_set(coarseEntity);
  for (; coarseBaseFunc < numInnerCorrectors; ++coarseBaseFunc)
  {
    assert(allLocalRHS[coarseBaseFunc]);
    CoarseBasisProduct eval(coarseBaseFunctionSet, local_diffusion_operator_, coarseBaseFunc);
    GDT::LocalFunctional::Codim0Integral<CoarseBasisProduct> local_rhs_functional(eval);
    auto& rhs_vector = allLocalRHS[coarseBaseFunc]->vector();
    rhs_functionals[coarseBaseFunc] = DSC::make_unique<RhsFunctionalType>(local_diffusion_operator_, rhs_vector,
                                                                          localSpace_, local_rhs_functional);
    system_assembler_.add(*rhs_functionals[coarseBaseFunc]);
  }

  // coarseBaseFunc == numInnerCorrectors
  //neumann corrector
  const auto local_neumann = DMP::getNeumannData()->transfer<MsFEMTraits::LocalEntityType>();
  GDT::Functionals::L2Face< decltype(local_neumann), CommonTraits::GdtVectorType, MsFEMTraits::LocalSpaceType >
      neumann_functional(local_neumann,
                         allLocalRHS[coarseBaseFunc]->vector(), localSpace_);
  system_assembler_.add(neumann_functional);

  coarseBaseFunc++;// coarseBaseFunc == 1 + numInnerCorrectors
  //dirichlet corrector
  typedef DirichletProduct DirichletEvaluationType;
  typedef GDT::Functionals::L2Volume<Problem::LocalDiffusionType, CommonTraits::GdtVectorType,
      MsFEMTraits::LocalSpaceType, MsFEMTraits::LocalGridViewType,
      DirichletEvaluationType > DirichletCorrectorFunctionalType;
  DirichletEvaluationType eval(dirichletExtensionLocal, local_diffusion_operator_);
  GDT::LocalFunctional::Codim0Integral<DirichletEvaluationType> dl_corrector_functional(eval);
  auto& dl_vector = allLocalRHS[coarseBaseFunc]->vector();
  DirichletCorrectorFunctionalType dl_func(local_diffusion_operator_, dl_vector,
                                           localSpace_, dl_corrector_functional);
  system_assembler_.add(dl_func);

  //dirichlet-0 for all rhs
  typedef GDT::ApplyOn::BoundaryEntities< MsFEMTraits::LocalGridViewType > OnLocalBoundaryEntities;
  system_assembler_.add(dirichletConstraints_, system_matrix_, new OnLocalBoundaryEntities());
  for (auto& rhs : allLocalRHS )
    system_assembler_.add(dirichletConstraints_, rhs->vector(), new OnLocalBoundaryEntities());

  system_assembler_.assemble();
}


void LocalProblemOperator::apply_inverse(const MsFEMTraits::LocalGridDiscreteFunctionType &current_rhs,
                                         MsFEMTraits::LocalGridDiscreteFunctionType &current_solution)
{
  if (!current_rhs.dofs_valid())
    DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");

  typedef BackendChooser<LocalGridDiscreteFunctionSpaceType>::InverseOperatorType LocalInverseOperatorType;
  const auto localProblemSolver = DSC::make_unique<LocalInverseOperatorType>(system_matrix_,
                                                                             current_rhs.space().communicator());
                                                                             /*1e-8, 1e-8, 20000,
                                               DSC_CONFIG_GET("msfem.localproblemsolver_verbose", false), solver,
                                               DSC_CONFIG_GET("preconditioner_type", std::string("sor")), 1);*/
  localProblemSolver->apply(current_rhs.vector(), current_solution.vector());

  if (!current_solution.dofs_valid())
    DUNE_THROW(Dune::InvalidStateException, "Current solution of the local msfem problem invalid!");
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {


#if 0 // alter dune-fem code
LocalProblemOperator::assemble_all_local_rhs
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
