#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/fem/misc/threads/domainthreaditerator.hh>
#include <memory>

#include <dune/stuff/fem/functions/integrals.hh>
#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/common/dirichletconstraints.hh"
#include "dune/multiscale/msfem/msfem_grid_specifier.hh"
#include "dune/multiscale/msfem/localproblems/subgrid-list.hh"
#include "dune/multiscale/msfem/localproblems/localsolutionmanager.hh"
#include "righthandside_assembler.hh"


void Dune::Multiscale::RightHandSideAssembler::assemble(const Dune::Multiscale::CommonTraits::FirstSourceType &f,
                                                        Dune::Multiscale::RightHandSideAssembler::DiscreteFunctionType &rhsVector) {
  rhsVector.clear();
  for (const auto& entity : rhsVector.space()) {
    const auto& geometry = entity.geometry();
    auto elementOfRHS = rhsVector.localFunction(entity);
    const auto baseSet = rhsVector.space().basisFunctionSet(entity);
    const auto quadrature = DSFe::make_quadrature(entity, rhsVector.space(), quadratureOrder);
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
  assert(false);//boundary treatment missing
//  const auto boundary = Problem::getModelData()->boundaryInfo();
//  const auto dirichlet_data = Problem::getDirichletData();
//  //! \TODO use the static thingies
//  DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType> constraints(*boundary, rhsVector.space());
//  constraints(*dirichlet_data, rhsVector);
}


void Dune::Multiscale::RightHandSideAssembler::assemble_hmm_lod(
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

    const auto& lagrangePointSet = rhsVector.space().lagrangePointSet(entity);

    for (const auto& intersection : DSC::intersectionRange(rhsVector.space().gridPart(), entity)) {
      if (Problem::isNeumannBoundary(intersection)) {
        const auto face = intersection.indexInInside();
        const auto faceQuadrature = DSFe::make_quadrature(intersection, rhsVector.space(), quadratureOrder);
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

    const auto loc_dirichlet_extension = dirichlet_extension.localFunction(entity);
    const auto quadrature = DSFe::make_quadrature(entity, rhsVector.space(), quadratureOrder);

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

void Dune::Multiscale::RightHandSideAssembler::assemble_for_MsFEM_symmetric(
    const Dune::Multiscale::CommonTraits::FirstSourceType &f,
    DMM::MacroMicroGridSpecifier &specifier, DMM::LocalGridList &subgrid_list,
    Dune::Multiscale::RightHandSideAssembler::DiscreteFunctionType &rhsVector) {
  auto diffusionPtr = Problem::getDiffusion();
  const auto& diffusion = *diffusionPtr;
  auto neumannDataPtr = Problem::getNeumannData();
  const auto& neumannData = *neumannDataPtr;

  rhsVector.clear();
  RangeType f_x;
  Fem::DomainDecomposedIteratorStorage< CommonTraits::GridPartType > threadIterators(rhsVector.space().gridPart());

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
  for (const auto& coarse_grid_entity : threadIterators) {
    const auto& coarseGeometry = coarse_grid_entity.geometry();
    auto rhsLocalFunction = rhsVector.localFunction(coarse_grid_entity);
    const auto numLocalBaseFunctions = rhsLocalFunction.numDofs();
    const auto& coarseBaseSet = specifier.coarseSpace().basisFunctionSet(coarse_grid_entity);

    // --------- add corrector contribution of right hand side --------------------------
    // Load local solutions
    MsFEM::LocalSolutionManager localSolutionManager(coarse_grid_entity, subgrid_list, specifier);
    localSolutionManager.loadLocalSolutions();
    auto& localSolutions = localSolutionManager.getLocalSolutions();
    assert(localSolutions.size() > 0);
    MsFEM::MsFEMTraits::LocalGridDiscreteFunctionType dirichletExtension("Dirichlet Extension",
                                                                         localSolutionManager.space());
    dirichletExtension.clear();
    Dune::Multiscale::copyDirichletValues(rhsVector.space(), dirichletExtension);

    // iterator for the micro grid ( grid for the reference element T_0 )
    const auto& subGrid = subgrid_list.getSubGrid(coarse_grid_entity);
    auto view = subGrid.leafView();
    for (const auto& localEntity : DSC::viewRange(view)) {    
      if (subgrid_list.covers(coarse_grid_entity, localEntity)) {
        // higher order quadrature, since A^{\epsilon} is highly variable
        const auto localQuadrature =
            DSFe::make_quadrature(localEntity, localSolutionManager.space());
        auto dirichletExtensionLF = dirichletExtension.localFunction(localEntity);

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

          // assemble intersection-part
          const auto& subGridPart = localSolutionManager.grid_part();
          for (const auto& intersection : DSC::intersectionRange(subGridPart.grid().leafView(), localEntity)) {
            if (Problem::isNeumannBoundary(intersection)) {
              const int orderOfIntegrand = (polynomialOrder - 1) + 2 * (polynomialOrder + 1);
              const int quadOrder = std::ceil((orderOfIntegrand + 1) / 2);
              const auto faceQuad = DSFe::make_quadrature(intersection, localSolutions[lsNum]->space(), quadOrder);
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
                coarseBaseSet.evaluateAll(xInCoarseLocal, phi_x_vec);
                for (auto i : DSC::valueRange(numLocalBaseFunctions)) {
                  assert((long long)i < (long long)phi_x_vec.size());
                  assert(iqP < localSolutionOnFace.size());
                  rhsLocalFunction[i] += factor * (neumannValue * (phi_x_vec[i] + localSolutionOnFace[iqP]));
                }
              }
            }
          }
        }

        // assemble element-part
        const auto& localGeometry = localEntity.geometry();
        for (size_t qP = 0; qP < localQuadrature.nop(); ++qP) {
          // local (barycentric) coordinates (with respect to entity)
          const auto& quadPoint = localQuadrature.point(qP);
          const auto quadPointGlobal = localGeometry.global(quadPoint);

          const double quadWeight = localQuadrature.weight(qP) * localGeometry.integrationElement(quadPoint);

          // evaluate gradient of basis function
          const auto quadPointLocalInCoarse = coarseGeometry.local(quadPointGlobal);
          std::vector<RangeType> coarseBaseEvals(numLocalBaseFunctions);
          std::vector<JacobianRangeType> coarseBaseJacs(numLocalBaseFunctions);
          coarseBaseSet.evaluateAll(quadPointLocalInCoarse, coarseBaseEvals);
          coarseBaseSet.jacobianAll(quadPointLocalInCoarse, coarseBaseJacs);

          for (int coarseBF = 0; coarseBF < numLocalBaseFunctions; ++coarseBF) {
            JacobianRangeType diffusive_flux(0.0);

            JacobianRangeType reconstructionGradPhi(coarseBaseJacs[coarseBF]);
            RangeType reconstructionPhi(coarseBaseEvals[coarseBF]);

            if (specifier.simplexCoarseGrid()) {
              assert(localSolutions.size() == dimension + localSolutionManager.numBoundaryCorrectors());
              for (const auto& i : DSC::valueRange(dimension))
                reconstructionPhi +=
                    coarseBaseJacs[coarseBF][0][i] * allLocalSolutionEvaluations[i][qP];
              //! @todo add the dirichlet and neumann-correctors!
            } else {
              assert(localSolutions.size() == numLocalBaseFunctions + localSolutionManager.numBoundaryCorrectors());
              // local corrector for coarse base func
              reconstructionPhi += allLocalSolutionEvaluations[coarseBF][qP];
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

            rhsLocalFunction[coarseBF] += quadWeight * (f_x * reconstructionPhi);
            rhsLocalFunction[coarseBF] -= quadWeight * (diffusive_flux[0] * reconstructionGradPhi[0]);
          }
        }
      }
    }
  } // for
  } // omp region

  // set dirichlet dofs to zero
  Dune::Multiscale::getConstraintsCoarse(rhsVector.space()).setValue(0.0, rhsVector);
}
