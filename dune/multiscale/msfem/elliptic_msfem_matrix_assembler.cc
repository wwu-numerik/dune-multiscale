#include <config.h>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <sstream>

#include "dune/multiscale/msfem/msfem_traits.hh"
#include "elliptic_msfem_matrix_assembler.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

DiscreteEllipticMsFEMOperator::DiscreteEllipticMsFEMOperator(
    MacroMicroGridSpecifierType& specifier, const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
    // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
    // n(T)-layers:
    MsFEMTraits::SubGridListType& subgrid_list, const DiffusionModel& diffusion_op)
  : specifier_(specifier)
  , coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
  , subgrid_list_(subgrid_list)
  , diffusion_operator_(diffusion_op)
  , petrovGalerkin_(false) {}

void DiscreteEllipticMsFEMOperator::assemble_matrix(MatrixType& global_matrix) const {
  // the local problem:
  // Let 'T' denote a coarse grid element and
  // let 'U(T)' denote the environment of 'T' that corresponds with the subgrid.

  // if Petrov-Galerkin-MsFEM
  if (petrovGalerkin_)
    DSC_LOG_INFO << "Assembling Petrov-Galerkin-MsFEM Matrix." << std::endl;
  else // if classical (symmetric) MsFEM
    DSC_LOG_INFO << "Assembling MsFEM Matrix." << std::endl;

  //!TODO diagonal stencil reicht
  global_matrix.reserve(DSFe::diagonalAndNeighborStencil(global_matrix));
  global_matrix.clear();

  const auto& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace_.gridPart().grid().leafIndexSet();

  Fem::DomainDecomposedIteratorStorage< CommonTraits::GridPartType > threadIterators(coarseDiscreteFunctionSpace_.gridPart());
  threadIterators.update();

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
  for (const auto& coarse_grid_entity : threadIterators) {
    const auto& coarse_grid_geometry = coarse_grid_entity.geometry();
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    const auto global_index_entity = coarseGridLeafIndexSet.index(coarse_grid_entity);

    DSFe::LocalMatrixProxy<MatrixType> local_matrix(global_matrix, coarse_grid_entity, coarse_grid_entity);

    const auto& coarse_grid_baseSet = local_matrix.domainBasisFunctionSet();
    const auto numMacroBaseFunctions = coarse_grid_baseSet.size();

    // Load local solutions
    Multiscale::MsFEM::LocalSolutionManager localSolutionManager(coarse_grid_entity, subgrid_list_, specifier_);
    localSolutionManager.loadLocalSolutions();
    const auto& localSolutions = localSolutionManager.getLocalSolutions();
    assert(localSolutions.size() > 0);
    std::vector<typename CoarseBaseFunctionSet::JacobianRangeType> gradientPhi(numMacroBaseFunctions);

    // iterator for the micro grid ( grid for the reference element T_0 )
    for (const auto& localGridEntity : localSolutionManager.getLocalDiscreteFunctionSpace()) {
      // check if "localGridEntity" (which is an entity of U(T)) is in T:
      // -------------------------------------------------------------------
      const auto& hostEntity = localSolutionManager.getSubGridPart().grid().getHostEntity<0>(localGridEntity);
      // ignore overlay elements
      if (global_index_entity == subgrid_list_.getEnclosingMacroCellIndex(hostEntity)) {
        assert(hostEntity->partitionType() == InteriorEntity);

        const auto& local_grid_geometry = localGridEntity.geometry();

        // higher order quadrature, since A^{\epsilon} is highly variable
        const auto localQuadrature =
            DSFe::make_quadrature(localGridEntity, localSolutionManager.getLocalDiscreteFunctionSpace());
        const auto numQuadraturePoints = localQuadrature.nop();

        // number of local solutions without the boundary correctors. Those are only needed for the right hand side
        const auto numLocalSolutions = localSolutions.size() - localSolutionManager.numBoundaryCorrectors();
        // evaluate the jacobians of all local solutions in all quadrature points
        std::vector<std::vector<JacobianRangeType>> allLocalSolutionEvaluations(
            numLocalSolutions, std::vector<JacobianRangeType>(localQuadrature.nop(), JacobianRangeType(0.0)));
        for (auto lsNum : DSC::valueRange(numLocalSolutions)) {
          auto& sll = localSolutions[lsNum];
          assert(sll.get());
          assert(sll->dofsValid());
          auto localFunction = sll->localFunction(localGridEntity);
          localFunction.evaluateQuadrature(localQuadrature, allLocalSolutionEvaluations[lsNum]);
        }

        for (size_t localQuadraturePoint = 0; localQuadraturePoint < numQuadraturePoints; ++localQuadraturePoint) {
          // local (barycentric) coordinates (with respect to entity)
          const auto& local_subgrid_point = localQuadrature.point(localQuadraturePoint);

          auto global_point_in_U_T = local_grid_geometry.global(local_subgrid_point);
          const double weight_local_quadrature = localQuadrature.weight(localQuadraturePoint) *
                                                 local_grid_geometry.integrationElement(local_subgrid_point);

          // evaluate the jacobian of the coarse grid base set
          const auto& local_coarse_point = coarse_grid_geometry.local(global_point_in_U_T);
          coarse_grid_baseSet.jacobianAll(local_coarse_point, gradientPhi);

          for (unsigned int i = 0; i < numMacroBaseFunctions; ++i) {
            for (unsigned int j = 0; j < numMacroBaseFunctions; ++j) {
              RangeType local_integral(0.0);

              // Compute the gradients of the i'th and j'th local problem solutions
              JacobianRangeType gradLocProbSoli(0.0), gradLocProbSolj(0.0);
              if (specifier_.simplexCoarseGrid()) {
                assert(allLocalSolutionEvaluations.size() == dimension);
                // ∇ Phi_H + ∇ Q( Phi_H ) = ∇ Phi_H + ∂_x1 Phi_H ∇Q( e_1 ) + ∂_x2 Phi_H ∇Q( e_2 )
                for (int k = 0; k < dimension; ++k) {
                  gradLocProbSoli.axpy(gradientPhi[i][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
                  gradLocProbSolj.axpy(gradientPhi[j][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
                }
              } else {
                assert(allLocalSolutionEvaluations.size() == numMacroBaseFunctions);
                gradLocProbSoli = allLocalSolutionEvaluations[i][localQuadraturePoint];
                gradLocProbSolj = allLocalSolutionEvaluations[j][localQuadraturePoint];
              }

              JacobianRangeType reconstructionGradPhii(gradientPhi[i]);
              reconstructionGradPhii += gradLocProbSoli;
              JacobianRangeType reconstructionGradPhij(gradientPhi[j]);
              reconstructionGradPhij += gradLocProbSolj;
              JacobianRangeType diffusive_flux(0.0);
              diffusion_operator_.diffusiveFlux(global_point_in_U_T, reconstructionGradPhii, diffusive_flux);
              if (petrovGalerkin_)
                local_integral += weight_local_quadrature * (diffusive_flux[0] * gradientPhi[j][0]);
              else
                local_integral += weight_local_quadrature * (diffusive_flux[0] * reconstructionGradPhij[0]);

              // add entries
              local_matrix.add(j, i, local_integral);
            }
          }
        }
      }
    }
  } // for
  } // omp region

  // set unit rows for dirichlet dofs
  Dune::Multiscale::getConstraintsCoarse(coarseDiscreteFunctionSpace_).applyToOperator(global_matrix);
  global_matrix.communicate();
} // assemble_matrix

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
