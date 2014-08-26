#include <config.h>

#include "coarse_rhs_functional.hh"

#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/functions/femadapter.hh>
#include <dune/stuff/common/profiler.hh>
#include <memory>

#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/stuff/functions/femadapter.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

size_t RhsCodim0Integral::numTmpObjectsRequired() const
{
  return numTmpObjectsRequired_;
}

void RhsCodim0Integral::apply(MsFEM::MsFEMTraits::LocalGridDiscreteFunctionType& dirichletExtension,
                              LocalSolutionManager &localSolutionManager,
                              const MsFEMTraits::LocalEntityType &localGridEntity,
                              const RhsCodim0Integral::TestLocalfunctionSetInterfaceType &testBase,
                              Dune::DynamicVector<CommonTraits::RangeFieldType> &ret,
                              std::vector<Dune::DynamicVector<CommonTraits::RangeFieldType> > &tmpLocalVectors) const
{
  auto& diffusion_operator = DMP::getDiffusion();
  const auto& coarse_scale_entity = testBase.entity();
  const auto& f = DMP::getSource();
  const auto& diffusion = DMP::getDiffusion();

  // quadrature
  typedef Dune::QuadratureRules<CommonTraits::DomainFieldType, CommonTraits::dimDomain> VolumeQuadratureRules;
  typedef Dune::QuadratureRule<CommonTraits::DomainFieldType, CommonTraits::dimDomain> VolumeQuadratureType;
  const size_t integrand_order = diffusion_operator->order() + testBase.order() + over_integrate_;
  assert(integrand_order < std::numeric_limits< int >::max());
  const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(localGridEntity.type(), int(integrand_order));
  // check matrix and tmp storage
  const size_t numLocalBaseFunctions = testBase.size();

  ret *= 0.0;
  assert(tmpLocalVectors.size() >= numTmpObjectsRequired_);

  const auto numQuadraturePoints = volumeQuadrature.size();
  const auto& localSolutions = localSolutionManager.getLocalSolutions();

  const auto numLocalSolutions = localSolutions.size();
  typedef CommonTraits::DiscreteFunctionSpaceType::BaseFunctionSetType::RangeType RangeType;
  typedef CommonTraits::DiscreteFunctionSpaceType::BaseFunctionSetType::JacobianRangeType JacobianRangeType;
  // evaluate the jacobians of all local solutions in all quadrature points
  std::vector<std::vector<JacobianRangeType>> allLocalSolutionJacobians(
        numLocalSolutions, std::vector<JacobianRangeType>(numQuadraturePoints, JacobianRangeType(0.0)));
  std::vector<std::vector<RangeType>> allLocalSolutionEvaluations(
        numLocalSolutions, std::vector<RangeType>(numQuadraturePoints, RangeType(0.0)));
  for (auto lsNum : DSC::valueRange(numLocalSolutions)) {
    const auto localFunction = localSolutions[lsNum]->local_function(localGridEntity);
    //      assert(localSolutionManager.space().indexSet().contains(localGridEntity));
    localFunction->jacobian(volumeQuadrature, allLocalSolutionJacobians[lsNum]);
    localFunction->evaluate(volumeQuadrature, allLocalSolutionEvaluations[lsNum]);
  }

  RangeType f_x;
  auto dirichletExtensionLF = dirichletExtension.local_function(localGridEntity);
  // loop over all quadrature points
  const auto quadPointEndIt = volumeQuadrature.end();
  std::size_t localQuadraturePoint = 0;
  for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt, ++localQuadraturePoint) {
    const auto x = quadPointIt->position();
    const auto coarseBaseJacs = testBase.jacobian(x);
    const auto coarseBaseEvals = testBase.evaluate(x);
    // integration factors
    const double integrationFactor = localGridEntity.geometry().integrationElement(x);
    const double quadratureWeight = quadPointIt->weight();
    const auto quadPointGlobal = localGridEntity.geometry().global(x);

    // compute integral
    for (size_t ii = 0; ii < numLocalBaseFunctions; ++ii) {
      auto& retRow = ret[ii];
      JacobianRangeType diffusive_flux(0.0);

      JacobianRangeType reconstructionGradPhi(coarseBaseJacs[ii]);
      RangeType reconstructionPhi(coarseBaseEvals[ii]);

      assert(localSolutions.size() == numLocalBaseFunctions + localSolutionManager.numBoundaryCorrectors());
      // local corrector for coarse base func
      reconstructionPhi += allLocalSolutionEvaluations[ii][localQuadraturePoint];
      // element part of boundary conditions
      JacobianRangeType directionOfFlux(0.0);
      //! @attention At this point we assume, that the quadrature points on the subgrid and hostgrid
      //! are the same (dirichletExtensionLF is a localfunction on the hostgrid, quadPoint stems from
      //! a quadrature on the subgrid)!!
      dirichletExtensionLF->jacobian(x, directionOfFlux);
      // add dirichlet-corrector
      directionOfFlux += allLocalSolutionJacobians[numLocalBaseFunctions + 1][localQuadraturePoint];
      // subtract neumann-corrector
      directionOfFlux -= allLocalSolutionJacobians[numLocalBaseFunctions][localQuadraturePoint];

      diffusion->diffusiveFlux(quadPointGlobal, directionOfFlux, diffusive_flux);
      reconstructionGradPhi += allLocalSolutionJacobians[ii][localQuadraturePoint];

      f->evaluate(quadPointGlobal, f_x);

      retRow += integrationFactor * quadratureWeight * (f_x * reconstructionPhi);
      retRow -= integrationFactor * quadratureWeight * (diffusive_flux[0] * reconstructionGradPhi[0]);
    } // compute integral
  } // loop over all quadrature points
}

std::vector<size_t> RhsCodim0Vector::numTmpObjectsRequired() const
{
  return { numTmpObjectsRequired_, localFunctional_.numTmpObjectsRequired() };
}

void RhsCodim0Vector::assembleLocal(const CommonTraits::GdtSpaceType &testSpace,
                                    const CommonTraits::EntityType &coarse_grid_entity,
                                    CommonTraits::GdtVectorType &systemVector,
                                    std::vector<std::vector<Dune::DynamicVector<CommonTraits::RangeFieldType> > > &tmpLocalVectorContainer,
                                    Dune::DynamicVector<size_t> &tmpIndices) const
{
  // check
  assert(tmpLocalVectorContainer.size() >= 2);
  assert(tmpLocalVectorContainer[0].size() >= numTmpObjectsRequired_);
  assert(tmpLocalVectorContainer[1].size() >= localFunctional_.numTmpObjectsRequired());

  Multiscale::MsFEM::LocalSolutionManager localSolutionManager(testSpace,
                                                               coarse_grid_entity,
                                                               localGridList_);
  localSolutionManager.load();
  const auto& localSolutions = localSolutionManager.getLocalSolutions();
  assert(localSolutions.size() > 0);

  MsFEM::MsFEMTraits::LocalGridDiscreteFunctionType dirichletExtension(localSolutionManager.space(), "Dirichlet Extension");
  //! \todo fill with actual values

  for (const auto& localGridEntity : DSC::viewRange(*localSolutionManager.space().grid_view())) {
    // ignore overlay elements
    if (!localGridList_.covers(coarse_grid_entity, localGridEntity))
      continue;

    // get and clear vector
    auto& localVector = tmpLocalVectorContainer[0][0];
    localVector *= 0.0;
    auto& tmpFunctionalVectors = tmpLocalVectorContainer[1];
    // apply local functional (result is in localVector)
    localFunctional_.apply(dirichletExtension, localSolutionManager, localGridEntity,
                           testSpace.base_function_set(coarse_grid_entity),
                           localVector,
                           tmpFunctionalVectors);
    // write local vector to global
    const size_t size = testSpace.mapper().numDofs(coarse_grid_entity);
    assert(tmpIndices.size() >= size);
    testSpace.mapper().globalIndices(coarse_grid_entity, tmpIndices);
    for (size_t ii = 0; ii < size; ++ii) {
      systemVector.add_to_entry(tmpIndices[ii], localVector[ii]);
    } // write local matrix to global
  }
}

CoarseRhsFunctional::CoarseRhsFunctional(const CoarseRhsFunctionalTraits::FunctionType &, CoarseRhsFunctional::VectorType &vec, const CoarseRhsFunctional::SpaceType &spc, LocalGridList &localGridList)
  : FunctionalBaseType(vec, spc)
  , AssemblerBaseType(spc)
  , local_assembler_(local_functional_, localGridList)
{
  this->add_codim0_assembler(local_assembler_, this->vector());
}


void CoarseRhsFunctional::assemble()
{
  AssemblerBaseType::assemble();
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#if 0 //alter referenzcde
rhsVector.clear();
RangeType f_x;

for (const auto& coarse_grid_entity : threadIterators) {
  const auto& coarseGeometry = coarse_grid_entity.geometry();
  auto rhsLocalFunction = rhsVector.localFunction(coarse_grid_entity);
  const auto numLocalBaseFunctions = rhsLocalFunction.numDofs();
  const auto& coarseBaseSet = coarse_space.basisFunctionSet(coarse_grid_entity);

  // --------- add corrector contribution of right hand side --------------------------
  // Load local solutions
  MsFEM::LocalSolutionManager localSolutionManager(coarse_space, coarse_grid_entity, subgrid_list);
  localSolutionManager.load();
  auto& localSolutions = localSolutionManager.getLocalSolutions();
  assert(localSolutions.size() > 0);
  MsFEM::MsFEMTraits::LocalGridDiscreteFunctionType dirichletExtension("Dirichlet Extension",
                                                                       localSolutionManager.space());
  dirichletExtension.clear();
  Dune::Multiscale::copyDirichletValues(rhsVector.space(), dirichletExtension);

  // iterator for the micro grid ( grid for the reference element T_0 )
  const auto& subGrid = subgrid_list.getSubGrid(coarse_grid_entity);
  auto view = subGrid.leafGridView();
  for (const auto& localEntity : DSC::viewRange(view)) {
    if (!subgrid_list.covers(coarse_grid_entity, localEntity))
      continue;

    // higher order quadrature, since A^{\epsilon} is highly variable
    const auto localQuadrature = DSFe::make_quadrature(localEntity, localSolutionManager.space());
    auto dirichletExtensionLF = dirichletExtension.localFunction(localEntity);

    // evaluate all local solutions and their jacobians in all quadrature points
    std::vector<std::vector<RangeType>> allLocalSolutionEvaluations(
          localSolutions.size(), std::vector<RangeType>(localQuadrature.nop(), 0.0));
      std::vector<std::vector<JacobianRangeType>> allLocalSolutionJacobians(
            localSolutions.size(), std::vector<JacobianRangeType>(localQuadrature.nop(), JacobianRangeType(0.0)));
      // add contributions for all inner correctors
      for (auto lsNum : DSC::valueRange(numLocalBaseFunctions))
      {
        auto localFunction = localSolutions[lsNum]->localFunction(localEntity);
        // this evaluates the local solutions in all quadrature points...
        localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
        // while this automatically evaluates their jacobians.
        localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionJacobians[lsNum]);

        // assemble intersection-part
        const auto& subGridPart = localSolutionManager.grid_part();
        for (const auto& intersection : DSC::intersectionRange(subGridPart.grid().leafGridView(), localEntity)) {
          if (DMP::is_neumann(intersection)) {
            const int orderOfIntegrand =
                (CommonTraits::polynomial_order - 1) + 2 * (CommonTraits::polynomial_order + 1);
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
              const auto& xGlobal = faceGeometry.global(xLocal);
              const auto& xInCoarseLocal = coarse_grid_entity.geometry().local(xGlobal);
              const double factor = faceGeometry.integrationElement(xLocal) * faceQuad.weight(iqP);

              neumannData.evaluate(xGlobal, neumannValue);
              coarseBaseSet.evaluateAll(xInCoarseLocal, phi_x_vec);
              assert((long long)lsNum < (long long)phi_x_vec.size());
              assert(iqP < localSolutionOnFace.size());
              rhsLocalFunction[lsNum] -= factor * (neumannValue * (phi_x_vec[lsNum] + localSolutionOnFace[iqP]));
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

          if (isSimplexGrid) {
            assert(localSolutions.size() == dimension + localSolutionManager.numBoundaryCorrectors());
            for (const auto& i : DSC::valueRange(dimension))
              reconstructionPhi += coarseBaseJacs[coarseBF][0][i] * allLocalSolutionEvaluations[i][qP];
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
  } // for

  // set dirichlet dofs to zero
  Dune::Multiscale::getConstraintsCoarse(rhsVector.space()).setValue(0.0, rhsVector);
#endif //0 //alter referenzcde
