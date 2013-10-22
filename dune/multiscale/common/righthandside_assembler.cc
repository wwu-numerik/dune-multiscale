#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/ranges.hh>
#include <memory>

#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/common/dirichletconstraints.hh"
#include "dune/multiscale/msfem/msfem_grid_specifier.hh"
#include "dune/multiscale/msfem/localproblems/subgrid-list.hh"
#include "dune/multiscale/msfem/localproblems/localsolutionmanager.hh"
#include "righthandside_assembler.hh"


void Dune::Multiscale::RightHandSideAssembler::assemble(const Dune::Multiscale::CommonTraits::FirstSourceType &f, Dune::Multiscale::RightHandSideAssembler::DiscreteFunctionType &rhsVector) {
  rhsVector.clear();
  for (const auto& entity : rhsVector.space()) {
    const auto& geometry = entity.geometry();
    auto elementOfRHS = rhsVector.localFunction(entity);
    const auto baseSet = rhsVector.space().basisFunctionSet(entity);
    const auto quadrature = make_quadrature(entity, rhsVector.space(), quadratureOrder);
    const auto numDofs = elementOfRHS.numDofs();
    for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
      // the return values:
      RangeType f_x;
      std::vector<RangeType> phi_x(numDofs);
      // to save: A \nabla PHI_H * \nabla phi_h;
      const auto det = geometry.integrationElement(quadrature.point(quadraturePoint));

      f.evaluate(geometry.global(quadrature.point(quadraturePoint)), f_x);
      baseSet.evaluateAll(quadrature[quadraturePoint], phi_x);

      for (int i = 0; i < numDofs; ++i) {
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x[i]);
      }
    }
  }
}


void Dune::Multiscale::RightHandSideAssembler::assemble(
    const Dune::Multiscale::CommonTraits::FirstSourceType &f,
    const Dune::Multiscale::CommonTraits::DiffusionType &A,
    const Dune::Multiscale::RightHandSideAssembler::DiscreteFunctionType &
        dirichlet_extension,
    const Dune::Multiscale::CommonTraits::NeumannBCType &neumann_bc,
    Dune::Multiscale::RightHandSideAssembler::DiscreteFunctionType &rhsVector) {
  rhsVector.clear();

  for (const auto& entity : rhsVector.space()) {

    const auto& geometry = entity.geometry();
    auto elementOfRHS = rhsVector.localFunction(entity);
    const auto baseSet = rhsVector.space().basisFunctionSet(entity);

    const int numDofs = elementOfRHS.numDofs();

    std::vector<RangeType> phi_x(numDofs);
    // gradient of base function and gradient of old_u_H
    std::vector<JacobianRangeType> grad_phi_x(numDofs);

    const auto loc_dirichlet_extension = dirichlet_extension.localFunction(entity);
    const auto quadrature = make_quadrature(entity, rhsVector.space(), quadratureOrder);

    const auto& lagrangePointSet = rhsVector.space().lagrangePointSet(entity);

    for (const auto& intersection : DSC::intersectionRange(rhsVector.space().gridPart(), entity)) {
      if (Problem::isNeumannBoundary(intersection)) {
        const auto face = intersection.indexInInside();
        const auto faceQuadrature = make_quadrature(intersection, rhsVector.space(), quadratureOrder);
        const auto numFaceQuadraturePoints = faceQuadrature.nop();

        for (auto faceQuadraturePoint : DSC::valueRange(numFaceQuadraturePoints)) {
          baseSet.evaluateAll(faceQuadrature[faceQuadraturePoint], phi_x);
          baseSet.jacobianAll(faceQuadrature[faceQuadraturePoint], grad_phi_x);

          const auto local_point_entity = faceQuadrature.point(faceQuadraturePoint);
          const auto global_point = geometry.global(local_point_entity);
          const auto local_point_face = intersection.geometry().local(global_point);

          RangeType neumann_value(0.0);
          neumann_bc.evaluate(global_point, neumann_value);

          const double face_weight = intersection.geometry().integrationElement(local_point_face) *
                                     faceQuadrature.weight(faceQuadraturePoint);

          for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
            elementOfRHS[lp] += neumann_value * face_weight * phi_x[lp];
          }
        }
      }
    }

    // the return values:
    RangeType f_x;

    JacobianRangeType gradient_dirichlet_extension;
    JacobianRangeType diffusive_flux_in_gradient_dirichlet_extension;

    for (auto quadraturePoint : DSC::valueRange(quadrature.nop())) {
      // local (barycentric) coordinates (with respect to entity)
      const auto& local_point = quadrature.point(quadraturePoint);
      const auto global_point = geometry.global(local_point);

      const double weight = geometry.integrationElement(local_point) * quadrature.weight(quadraturePoint);
      // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
      f.evaluate(global_point, f_x);
      // evaluate the current base function at the current quadrature point and save its value in 'z':
      baseSet.evaluateAll(quadrature[quadraturePoint], phi_x); // i = i'te Basisfunktion;
      // evaluate the gradient of the current base function at the current quadrature point and save its value in
      // 'returnGradient':
      baseSet.jacobianAll(quadrature[quadraturePoint], grad_phi_x);
      // get gradient of dirichlet extension:
      loc_dirichlet_extension.jacobian(quadrature[quadraturePoint], gradient_dirichlet_extension);
      A.diffusiveFlux(global_point, gradient_dirichlet_extension, diffusive_flux_in_gradient_dirichlet_extension);

      for (int i = 0; i < numDofs; ++i) {
        elementOfRHS[i] += weight * (f_x * phi_x[i]);
        elementOfRHS[i] -= weight * (diffusive_flux_in_gradient_dirichlet_extension[0] * grad_phi_x[i][0]);
      }
    }
  }
}

void Dune::Multiscale::RightHandSideAssembler::assemble_for_MsFEM_symmetric(const Dune::Multiscale::CommonTraits::FirstSourceType &f, Dune::Multiscale::MsFEM::MacroMicroGridSpecifier &specifier, Dune::Multiscale::MsFEM::SubGridList &subgrid_list, Dune::Multiscale::RightHandSideAssembler::DiscreteFunctionType &rhsVector) {

  // gather some problem data
  auto diffusionPtr = Problem::getDiffusion();
  const auto& diffusion = *diffusionPtr;
  auto neumannDataPtr = Problem::getNeumannData();
  const auto& neumannData = *neumannDataPtr;

  DiscreteFunctionType dirichletExtension("Dirichlet Extension", specifier.fineSpace());
  dirichletExtension.clear();
  Dune::Multiscale::copyDirichletValues(rhsVector.space(), dirichletExtension);

  // set rhsVector to zero:
  rhsVector.clear();
  const auto& coarseGridLeafIndexSet = specifier.coarseSpace().gridPart().grid().leafIndexSet();
  RangeType f_x;
  for (const auto& coarse_grid_entity : rhsVector.space()) {
    const auto coarseEntityIndex = coarseGridLeafIndexSet.index(coarse_grid_entity);

    const auto& coarseGeometry = coarse_grid_entity.geometry();
    auto rhsLocalFunction = rhsVector.localFunction(coarse_grid_entity);
    const auto numLocalBaseFunctions = rhsLocalFunction.numDofs();
    const auto& coarse_grid_baseSet = specifier.coarseSpace().basisFunctionSet(coarse_grid_entity);

    // --------- add standard contribution of right hand side -------------------------
    {
      const auto quadrature = make_quadrature(coarse_grid_entity, rhsVector.space(), quadratureOrder);
      std::vector<RangeType> phi_x_vec(numLocalBaseFunctions);
      const auto numQuadraturePoints = quadrature.nop();
      for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
        const double det = coarseGeometry.integrationElement(quadrature.point(quadraturePoint));
        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_y':
        f.evaluate(coarseGeometry.global(quadrature.point(quadraturePoint)), f_x);
        coarse_grid_baseSet.evaluateAll(quadrature[quadraturePoint], phi_x_vec);
        for (int i = 0; i < numLocalBaseFunctions; ++i) {
          rhsLocalFunction[i] += det * quadrature.weight(quadraturePoint) * (f_x * phi_x_vec[i]);
        }
      }
    }
    // ----------------------------------------------------------------------------------

    // --------- add corrector contribution of right hand side --------------------------
    // Load local solutions
    MsFEM::LocalSolutionManager localSolutionManager(coarse_grid_entity, subgrid_list, specifier);
    localSolutionManager.loadLocalSolutions();
    auto& localSolutions = localSolutionManager.getLocalSolutions();
    assert(localSolutions.size() > 0);

    // iterator for the micro grid ( grid for the reference element T_0 )
    const auto& subGrid = subgrid_list.getSubGrid(coarse_grid_entity);
    for (const auto& localEntity : DSC::viewRange(subGrid.leafView())) {
      const auto& hostCell = subGrid.getHostEntity<0>(localEntity);
      const auto enclosingCoarseCellIndex = subgrid_list.getEnclosingMacroCellIndex(hostCell);
      auto dirichletExtensionLF = dirichletExtension.localFunction(*hostCell);
      if (enclosingCoarseCellIndex == coarseEntityIndex) {
        // higher order quadrature, since A^{\epsilon} is highly variable
        const auto localQuadrature =
            make_quadrature(localEntity, localSolutionManager.getLocalDiscreteFunctionSpace());

        // evaluate all local solutions and their jacobians in all quadrature points
        std::vector<std::vector<RangeType>> allLocalSolutionEvaluations(
                                              localSolutions.size(), std::vector<RangeType>(localQuadrature.nop(), 0.0));
        std::vector<std::vector<JacobianRangeType>> allLocalSolutionJacobians(
                                                      localSolutions.size(), std::vector<JacobianRangeType>(localQuadrature.nop(), JacobianRangeType(0.0)));
        for (auto lsNum : DSC::valueRange(localSolutions.size())) {
          auto localFunction = localSolutions[lsNum]->localFunction(localEntity);
          // this evaluates the local solutions in all quadrature points...
          localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
          // while this automatically evaluates their jacobians.
          localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionJacobians[lsNum]);

          const auto& subGridPart = localSolutionManager.getSubGridPart();
          for (const auto& intersection : DSC::intersectionRange(subGridPart.grid().leafView(), localEntity)) {
            if (Problem::isNeumannBoundary(intersection)) {
              const int orderOfIntegrand = (polynomialOrder - 1) + 2 * (polynomialOrder + 1);
              const int quadOrder = std::ceil((orderOfIntegrand + 1) / 2);
              // get type of face quadrature. Is done in this scope because Patricks methods use another type.
              const auto faceQuad = make_quadrature(intersection, localSolutions[lsNum]->space(), quadOrder);
              RangeType neumannValue(0.0);
              const auto numQuadPoints = faceQuad.nop();
              // loop over all quadrature points

              std::vector<RangeType> phi_x_vec(numLocalBaseFunctions);
              std::vector<RangeType> localSolutionOnFace(numLocalBaseFunctions);
              localFunction.evaluateQuadrature(faceQuad, localSolutionOnFace);

              for (unsigned int iqP = 0; iqP < numQuadPoints; ++iqP) {
                // get local coordinate of quadrature point
                const auto& xLocal = faceQuad.localPoint(iqP);
                const auto& faceGeometry = intersection.geometry();

                // the following does not work because subgrid does not implement geometryInInside()
                // const auto& insideGeometry    = intersection.geometryInInside();
                // const typename FaceQuadratureType::CoordinateType& xInInside = insideGeometry.global(xLocal);
                // therefore, we have to do stupid things:
                const auto& xGlobal = faceGeometry.global(xLocal);
                const auto& xInCoarseLocal = coarse_grid_entity.geometry().local(xGlobal);
                const double factor = faceGeometry.integrationElement(xLocal) * faceQuad.weight(iqP);

                neumannData.evaluate(xGlobal, neumannValue);
                coarse_grid_baseSet.evaluateAll(xInCoarseLocal, phi_x_vec);
                for (auto i : DSC::valueRange(numLocalBaseFunctions)) {
                  assert((long long)i < (long long)phi_x_vec.size());
                  assert(iqP < localSolutionOnFace.size());
                  rhsLocalFunction[i] += factor * (neumannValue * (phi_x_vec[i] + localSolutionOnFace[iqP]));
                }
              }
            }
          }
        }

        const auto& localGeometry = localEntity.geometry();
        RangeType corrector_phi_x;
        for (size_t qP = 0; qP < localQuadrature.nop(); ++qP) {
          // local (barycentric) coordinates (with respect to entity)
          const auto& quadPoint = localQuadrature.point(qP);
          const auto quadPointGlobal = localGeometry.global(quadPoint);

          const double quadWeight = localQuadrature.weight(qP) * localGeometry.integrationElement(quadPoint);

          // evaluate gradient of basis function
          const auto quadPointLocalInCoarse = coarseGeometry.local(quadPointGlobal);
          std::vector<JacobianRangeType> gradient_Phi_vec(numLocalBaseFunctions);
          coarse_grid_baseSet.jacobianAll(quadPointLocalInCoarse, gradient_Phi_vec);

          for (int coarseBF = 0; coarseBF < numLocalBaseFunctions; ++coarseBF) {
            JacobianRangeType diffusive_flux(0.0);

            JacobianRangeType reconstructionGradPhi(gradient_Phi_vec[coarseBF]);

            if (specifier.simplexCoarseGrid()) {
              assert(localSolutions.size() == GridSelector::dimgrid + localSolutionManager.numBoundaryCorrectors());
              DUNE_THROW(NotImplemented, "Boundary values are not implemented for simplex grids yet!");
            } else {
              assert(localSolutions.size() == numLocalBaseFunctions + localSolutionManager.numBoundaryCorrectors());
              // local corrector for coarse base func
              corrector_phi_x = allLocalSolutionEvaluations[coarseBF][qP];
              // element part of boundary conditions
              JacobianRangeType directionOfFlux(0.0);
              //! @attention At this point we assume, that the quadrature points on the subgrid and hostgrid
              //! are the same (dirichletExtensionLF is a localfunction on the hostgrid, quadPoint stems from
              //! a quadrature on the subgrid)!!
              dirichletExtensionLF.jacobian(quadPoint, directionOfFlux);
              // add dirichlet-corrector
              directionOfFlux += allLocalSolutionJacobians[numLocalBaseFunctions + 1][qP];
              // subtract neumann-corrector
              directionOfFlux -= allLocalSolutionJacobians[numLocalBaseFunctions][qP];

              diffusion.diffusiveFlux(quadPointGlobal, directionOfFlux, diffusive_flux);
              reconstructionGradPhi += allLocalSolutionJacobians[coarseBF][qP];
            }
            f.evaluate(quadPointGlobal, f_x);
            double val = quadWeight * (f_x * corrector_phi_x);
            rhsLocalFunction[coarseBF] += val;
            val = quadWeight * (diffusive_flux[0] * reconstructionGradPhi[0]);
            rhsLocalFunction[coarseBF] -= val;
          }
        }
      }
    }
  }

  // set dirichlet dofs to zero
  Dune::Multiscale::getConstraintsCoarse(rhsVector.space()).setValue(0.0, rhsVector);
}
