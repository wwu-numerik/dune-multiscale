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
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>

namespace Dune {
namespace Multiscale {

size_t RhsCodim0Integral::numTmpObjectsRequired() const { return numTmpObjectsRequired_; }

void RhsCodim0Integral::apply(MsFEMTraits::LocalGridDiscreteFunctionType& dirichletExtension,
                              LocalSolutionManager& localSolutionManager,
                              const MsFEMTraits::LocalEntityType& localGridEntity,
                              const RhsCodim0Integral::TestLocalfunctionSetInterfaceType& testBase,
                              Dune::DynamicVector<CommonTraits::RangeFieldType>& ret,
                              std::vector<Dune::DynamicVector<CommonTraits::RangeFieldType>>& /*tmpLocalVectors*/) const {
  const auto& f = DMP::getSource();
  const auto& diffusion = DMP::getDiffusion();

  // quadrature
  typedef Dune::QuadratureRules<CommonTraits::DomainFieldType, CommonTraits::dimDomain> VolumeQuadratureRules;
  typedef Dune::QuadratureRule<CommonTraits::DomainFieldType, CommonTraits::dimDomain> VolumeQuadratureType;
  const size_t integrand_order = diffusion->order() + testBase.order() + over_integrate_;
  assert(integrand_order < std::numeric_limits<int>::max());
  const VolumeQuadratureType& volumeQuadrature =
      VolumeQuadratureRules::rule(localGridEntity.type(), int(integrand_order));
  // check matrix and tmp storage
  const size_t numLocalBaseFunctions = testBase.size();

  ret *= 0.0;

  const auto numQuadraturePoints = volumeQuadrature.size();
  const auto& localSolutions = localSolutionManager.getLocalSolutions();

  const auto numLocalSolutions = localSolutions.size();
  typedef CommonTraits::SpaceType::BaseFunctionSetType::RangeType RangeType;
  typedef CommonTraits::SpaceType::BaseFunctionSetType::JacobianRangeType JacobianRangeType;
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
  const auto dirichletExtensionLF = dirichletExtension.local_function(localGridEntity);
  // loop over all quadrature points
  const auto quadPointEndIt = volumeQuadrature.end();
  std::size_t localQuadraturePoint = 0;
  for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt;
       ++quadPointIt, ++localQuadraturePoint) {
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
  }   // loop over all quadrature points
}

std::vector<size_t> RhsCodim0Vector::numTmpObjectsRequired() const {
  return {numTmpObjectsRequired_, localFunctional_.numTmpObjectsRequired()};
}

void RhsCodim0Vector::assembleLocal(
    const CommonTraits::SpaceType& testSpace, const CommonTraits::EntityType& coarse_grid_entity,
    CommonTraits::GdtVectorType& systemVector,
    std::vector<std::vector<Dune::DynamicVector<CommonTraits::RangeFieldType>>>& tmpLocalVectorContainer,
    Dune::DynamicVector<size_t>& tmpIndices) const {
  // check
  assert(tmpLocalVectorContainer.size() >= 2);
  assert(tmpLocalVectorContainer[0].size() >= numTmpObjectsRequired_);
  assert(tmpLocalVectorContainer[1].size() >= localFunctional_.numTmpObjectsRequired());

  Multiscale::LocalSolutionManager localSolutionManager(testSpace, coarse_grid_entity, localGridList_);
  localSolutionManager.load();
  const auto& localSolutions = localSolutionManager.getLocalSolutions();
  assert(localSolutions.size() > 0);

  MsFEMTraits::LocalGridDiscreteFunctionType dirichletExtension(localSolutionManager.space(), "Dirichlet Extension");
  //! \todo fill with actual values

  for (const auto& localGridEntity : DSC::entityRange(localSolutionManager.space().grid_view())) {
    // ignore overlay elements
    if (!localGridList_.covers(coarse_grid_entity, localGridEntity))
      continue;

    // get and clear vector
    auto& localVector = tmpLocalVectorContainer[0][0];
    localVector *= 0.0;
    auto& tmpFunctionalVectors = tmpLocalVectorContainer[1];
    // apply local functional (result is in localVector)
    localFunctional_.apply(dirichletExtension, localSolutionManager, localGridEntity,
                           testSpace.base_function_set(coarse_grid_entity), localVector, tmpFunctionalVectors);
    // write local vector to global
    const size_t size = testSpace.mapper().numDofs(coarse_grid_entity);
    assert(tmpIndices.size() >= size);
    testSpace.mapper().globalIndices(coarse_grid_entity, tmpIndices);
    for (size_t ii = 0; ii < size; ++ii) {
      systemVector.add_to_entry(tmpIndices[ii], localVector[ii]);
    } // write local matrix to global
  }
}

CoarseRhsFunctional::CoarseRhsFunctional(CoarseRhsFunctional::VectorType& vec,
                                         const CoarseRhsFunctional::SpaceType& spc, LocalGridList& localGridList, const CommonTraits::InteriorGridViewType& interior)
  : FunctionalBaseType(vec, spc, interior)
  , AssemblerBaseType(spc, interior)
  , local_assembler_(local_functional_, localGridList) {
  this->add_codim0_assembler(local_assembler_, this->vector());
}

void CoarseRhsFunctional::assemble() { AssemblerBaseType::assemble(); }

} // namespace Multiscale {
} // namespace Dune {
