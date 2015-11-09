#include <config.h>

#include "coarse_scale_assembler.hh"

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>
#include <sstream>

namespace Dune {
namespace Multiscale {

MsFEMCodim0Integral::MsFEMCodim0Integral(const Problem::DiffusionBase& diffusion, const size_t over_integrate)
  : over_integrate_(over_integrate)
  , diffusion_(diffusion) {}

size_t MsFEMCodim0Integral::numTmpObjectsRequired() const { return numTmpObjectsRequired_; }

void MsFEMCodim0Integral::apply(
    LocalSolutionManager& localSolutionManager, const MsFEMTraits::LocalEntityType& localGridEntity,
    const MsFEMCodim0Integral::TestLocalfunctionSetInterfaceType& testBase,
    const MsFEMCodim0Integral::AnsatzLocalfunctionSetInterfaceType& ansatzBase,
    Dune::DynamicMatrix<CommonTraits::RangeFieldType>& ret,
    std::vector<Dune::DynamicMatrix<CommonTraits::RangeFieldType>>& /*tmpLocalMatrices*/) const {
  const auto& diffusion_operator = diffusion_;

  // quadrature
  typedef Dune::QuadratureRules<CommonTraits::DomainFieldType, CommonTraits::dimDomain> VolumeQuadratureRules;
  typedef Dune::QuadratureRule<CommonTraits::DomainFieldType, CommonTraits::dimDomain> VolumeQuadratureType;
  const size_t integrand_order = diffusion_operator.order() + ansatzBase.order() + testBase.order() + over_integrate_;
  assert(integrand_order < std::numeric_limits<int>::max());
  const VolumeQuadratureType& volumeQuadrature =
      VolumeQuadratureRules::rule(localGridEntity.type(), int(integrand_order));
  // check matrix and tmp storage
  const size_t rows = testBase.size();
  const size_t cols = ansatzBase.size();
  ret *= 0.0;
  assert(ret.rows() >= rows);
  assert(ret.cols() >= cols);

  const auto numQuadraturePoints = volumeQuadrature.size();
  const auto& localSolutions = localSolutionManager.getLocalSolutions();
  // number of local solutions without the boundary correctors. Those are only needed for the right hand side
  const auto numLocalSolutions = localSolutions.size() - localSolutionManager.numBoundaryCorrectors();
  typedef CommonTraits::SpaceType::BaseFunctionSetType::RangeType RangeType;
  typedef CommonTraits::SpaceType::BaseFunctionSetType::JacobianRangeType JacobianRangeType;
  // evaluate the jacobians of all local solutions in all quadrature points
  std::vector<std::vector<JacobianRangeType>> allLocalSolutionEvaluations(
      numLocalSolutions, std::vector<JacobianRangeType>(numQuadraturePoints, RangeType(0.0)));
  for (auto lsNum : DSC::valueRange(numLocalSolutions)) {
    const auto localFunction = localSolutions[lsNum]->local_function(localGridEntity);
    //      assert(localSolutionManager.space().indexSet().contains(localGridEntity));
    localFunction->jacobian(volumeQuadrature, allLocalSolutionEvaluations[lsNum]);
  }

  // loop over all quadrature points
  const auto quadPointEndIt = volumeQuadrature.end();
  std::size_t localQuadraturePoint = 0;
  for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt;
       ++quadPointIt, ++localQuadraturePoint) {
    const auto x = quadPointIt->position();
    const auto coarseBaseJacs = testBase.jacobian(x);
    // integration factors
    const double integrationFactor = localGridEntity.geometry().integrationElement(x);
    const double quadratureWeight = quadPointIt->weight();
    const auto global_point_in_U_T = localGridEntity.geometry().global(x);
    CommonTraits::DiffusionFunctionBaseType::RangeType diffusion_eval;
    diffusion_operator.evaluate(global_point_in_U_T, diffusion_eval);
    // compute integral
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        // Compute the gradients of the i'th and j'th local problem solutions
        assert(allLocalSolutionEvaluations.size() == rows /*numMacroBaseFunctions*/);
        const auto& gradLocProbSoli = allLocalSolutionEvaluations[ii][localQuadraturePoint];
        const auto& gradLocProbSolj = allLocalSolutionEvaluations[jj][localQuadraturePoint];

        auto reconstructionGradPhii = coarseBaseJacs[ii];
        reconstructionGradPhii += gradLocProbSoli;
        auto reconstructionGradPhij = coarseBaseJacs[jj];
        reconstructionGradPhij += gradLocProbSolj;

        const auto& arg = reconstructionGradPhii[0];
        CommonTraits::DiffusionFunctionBaseType::RangeType::row_type diffusive_flux;
        diffusion_eval.mv(arg, diffusive_flux);
        const RangeType local_integral = diffusive_flux * reconstructionGradPhij[0];
        //! TODO check indexing. Correct wrt pre-gdt, but still
        ret[jj][ii] += local_integral * integrationFactor * quadratureWeight;
      }
    } // compute integral
  }   // loop over all quadrature points
}

MsFemCodim0Matrix::MsFemCodim0Matrix(const MsFemCodim0Matrix::LocalOperatorType& op, LocalGridList& localGridList)
  : localOperator_(op)
  , localGridList_(localGridList) {}

const MsFemCodim0Matrix::LocalOperatorType& MsFemCodim0Matrix::localOperator() const { return localOperator_; }

std::vector<size_t> MsFemCodim0Matrix::numTmpObjectsRequired() const {
  return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
}

void MsFemCodim0Matrix::assembleLocal(
    const CommonTraits::SpaceType& testSpace, const CommonTraits::SpaceType& ansatzSpace,
    const CommonTraits::EntityType& coarse_grid_entity, CommonTraits::LinearOperatorType& systemMatrix,
    std::vector<std::vector<Dune::DynamicMatrix<CommonTraits::RangeFieldType>>>& tmpLocalMatricesContainer,
    std::vector<Dune::DynamicVector<size_t>>& tmpIndicesContainer) const {
  // check
  assert(tmpLocalMatricesContainer.size() >= 1);
  assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
  assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
  assert(tmpIndicesContainer.size() >= 2);
  // get and clear matrix

  Multiscale::LocalSolutionManager localSolutionManager(testSpace, coarse_grid_entity, localGridList_);
  localSolutionManager.load();
  const auto& localSolutions = localSolutionManager.getLocalSolutions();
  assert(localSolutions.size() > 0);

  for (const auto& localGridEntity : DSC::entityRange(localSolutionManager.space().grid_view())) {
    // ignore overlay elements
    if (!localGridList_.covers(coarse_grid_entity, localGridEntity))
      continue;

    auto& localMatrix = tmpLocalMatricesContainer[0][0];
    localMatrix *= 0.0;
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // apply local operator (result is in localMatrix)
    localOperator_.apply(localSolutionManager, localGridEntity, testSpace.base_function_set(coarse_grid_entity),
                         ansatzSpace.base_function_set(coarse_grid_entity), localMatrix, tmpOperatorMatrices);
    // write local matrix to global
    auto& globalRows = tmpIndicesContainer[0];
    auto& globalCols = tmpIndicesContainer[1];
    const size_t rows = testSpace.mapper().numDofs(coarse_grid_entity);
    const size_t cols = ansatzSpace.mapper().numDofs(coarse_grid_entity);
    assert(globalRows.size() >= rows);
    assert(globalCols.size() >= cols);
    testSpace.mapper().globalIndices(coarse_grid_entity, globalRows);
    ansatzSpace.mapper().globalIndices(coarse_grid_entity, globalCols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& localRow = localMatrix[ii];
      const size_t globalII = globalRows[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const size_t globalJJ = globalCols[jj];
        systemMatrix.add_to_entry(globalII, globalJJ, localRow[jj]);
      }
    } // write local matrix to global
  }
}

// constructor

} // namespace Multiscale {
} // namespace Dune {
