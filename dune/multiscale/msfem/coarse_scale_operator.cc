#include <config.h>

#include "coarse_scale_operator.hh"

#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>
#include <sstream>



namespace Dune {
namespace Multiscale {
namespace MsFEM {


Stuff::LA::SparsityPatternDefault CoarseScaleOperator::pattern(const CoarseScaleOperator::RangeSpaceType &range_space, const CoarseScaleOperator::SourceSpaceType &source_space, const CoarseScaleOperator::GridViewType &grid_view)
{
  return range_space.compute_volume_pattern(grid_view, source_space);
}

CoarseScaleOperator::CoarseScaleOperator(const CoarseScaleOperator::SourceSpaceType &src_spc, LocalGridList &localGridList)
  : OperatorBaseType(global_matrix_, src_spc)
  , AssemblerBaseType(src_spc)
  , global_matrix_(src_spc.mapper().size(), src_spc.mapper().size(), CommonTraits::EllipticOperatorType::pattern(src_spc))
  , local_assembler_(local_operator_, localGridList)
{
  this->add_codim0_assembler(local_assembler_, this->matrix());
}

void CoarseScaleOperator::assemble()
{
  AssemblerBaseType::assemble();
}

void CoarseScaleOperator::apply_inverse(const CoarseScaleOperator::CoarseDiscreteFunction& rhs,
                                        CoarseScaleOperator::CoarseDiscreteFunction& solution) {
  BOOST_ASSERT_MSG(rhs.dofs_valid(), "Coarse scale RHS DOFs need to be valid!");
  DSC_PROFILER.startTiming("msfem.solveCoarse");
  const typename BackendChooser<CoarseDiscreteFunctionSpace>::InverseOperatorType inverse(
        global_matrix_, rhs.space().communicator());/*, 1e-8, 1e-8, DSC_CONFIG_GET("msfem.solver.iterations", rhs.size()),
      DSC_CONFIG_GET("msfem.solver.verbose", false), "bcgs",
      DSC_CONFIG_GET("msfem.solver.preconditioner_type", std::string("sor")));*/
  inverse.apply(rhs.vector(), solution.vector());
  if (!solution.dofs_valid())
    DUNE_THROW(InvalidStateException, "Degrees of freedom of coarse solution are not valid!");
  DSC_PROFILER.stopTiming("msfem.solveCoarse");
  DSC_LOG_DEBUG << "Time to solve coarse MsFEM problem: " << DSC_PROFILER.getTiming("msfem.solveCoarse") << "ms."
               << std::endl;
}

CoarseScaleOperator::MatrixType &CoarseScaleOperator::system_matrix() { return global_matrix_; }

const CoarseScaleOperator::MatrixType &CoarseScaleOperator::system_matrix() const { return global_matrix_; }

// constructor

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#if 0 //alter referenzcode
CoarseScaleOperator::CoarseScaleOperator(const CoarseDiscreteFunctionSpace& coarseDiscreteFunctionSpace,
                                         LocalGridList& subgrid_list, const CommonTraits::DiffusionType& diffusion_op)
  : coarseDiscreteFunctionSpace_(coarseDiscreteFunctionSpace)
  , subgrid_list_(subgrid_list)
  , diffusion_operator_(diffusion_op)
  , petrovGalerkin_(false)
  , global_matrix_(coarseDiscreteFunctionSpace.mapper().size(), coarseDiscreteFunctionSpace.mapper().size(),
                   CommonTraits::EllipticOperatorType::pattern(coarseDiscreteFunctionSpace))
{
    DSC_PROFILER.startTiming("msfem.assembleMatrix");

  const bool is_simplex_grid = DSG::is_simplex_grid(coarseDiscreteFunctionSpace_);
  assert(!is_simplex_grid); //special casing isn't ported yet.

  for (const auto& coarse_grid_entity : DSC::viewRange(*coarseDiscreteFunctionSpace_.grid_view())) {
    int cacheCounter = 0;
    const auto& coarse_grid_geometry = coarse_grid_entity.geometry();
    assert(coarse_grid_entity.partitionType() == InteriorEntity);

    //      DSFe::LocalMatrixProxy<MatrixType> local_matrix(global_matrix_, coarse_grid_entity, coarse_grid_entity);
    auto local_matrix = nullptr;

    const auto& coarse_grid_baseSet = local_matrix.domainBasisFunctionSet();
    const auto numMacroBaseFunctions = coarse_grid_baseSet.size();

    Multiscale::MsFEM::LocalSolutionManager localSolutionManager(coarseDiscreteFunctionSpace_,
                                                                 coarse_grid_entity,
                                                                 subgrid_list_);
    localSolutionManager.load();
    const auto& localSolutions = localSolutionManager.getLocalSolutions();
    assert(localSolutions.size() > 0);

    for (const auto& localGridEntity : DSC::viewRange(*localSolutionManager.space().grid_view())) {
      // ignore overlay elements
      if (!subgrid_list_.covers(coarse_grid_entity, localGridEntity))
        continue;
      const auto& local_grid_geometry = localGridEntity.geometry();

      // higher order quadrature, since A^{\epsilon} is highly variable
      const auto localQuadrature = DSFe::make_quadrature(localGridEntity, localSolutionManager.space());
      const auto numQuadraturePoints = localQuadrature.nop();

      // number of local solutions without the boundary correctors. Those are only needed for the right hand side
      const auto numLocalSolutions = localSolutions.size() - localSolutionManager.numBoundaryCorrectors();
      // evaluate the jacobians of all local solutions in all quadrature points
      std::vector<std::vector<JacobianRangeType>> allLocalSolutionEvaluations(
            numLocalSolutions, std::vector<JacobianRangeType>(localQuadrature.nop(), JacobianRangeType(0.0)));
      for (auto lsNum : DSC::valueRange(numLocalSolutions)) {
        auto& sll = localSolutions[lsNum];
        assert(sll.get());
        //! \todo re-enable
        //            assert(localSolutionManager.space().indexSet().contains(localGridEntity));
        auto localFunction = sll->local_function(localGridEntity);
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
            CoarseDiscreteFunctionSpace::RangeType local_integral(0.0);

            // Compute the gradients of the i'th and j'th local problem solutions
            JacobianRangeType gradLocProbSoli(0.0), gradLocProbSolj(0.0);
            if (is_simplex_grid) {
              assert(allLocalSolutionEvaluations.size() == CommonTraits::GridType::dimension);
              // ∇ Phi_H + ∇ Q( Phi_H ) = ∇ Phi_H + ∂_x1 Phi_H ∇Q( e_1 ) + ∂_x2 Phi_H ∇Q( e_2 )
              for (int k = 0; k < CommonTraits::GridType::dimension; ++k) {
                gradLocProbSoli.axpy(gradientPhi[i][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
                gradLocProbSolj.axpy(gradientPhi[j][0][k], allLocalSolutionEvaluations[k][localQuadraturePoint]);
              }
            } else {
              assert(allLocalSolutionEvaluations.size() == numMacroBaseFunctions);
              gradLocProbSoli = allLocalSolutionEvaluations[i][localQuadraturePoint];
              gradLocProbSolj = allLocalSolutionEvaluations[j][localQuadraturePoint];
            }

            JacobianRangeType reconstructionGradPhii(coarseBaseJacs_[cacheCounter][i]);
            reconstructionGradPhii += gradLocProbSoli;
            JacobianRangeType reconstructionGradPhij(coarseBaseJacs_[cacheCounter][j]);
            reconstructionGradPhij += gradLocProbSolj;
            JacobianRangeType diffusive_flux(0.0);
            diffusion_operator_.diffusiveFlux(global_point_in_U_T, reconstructionGradPhii, diffusive_flux);
            if (petrovGalerkin_)
              local_integral += weight_local_quadrature * (diffusive_flux[0] * coarseBaseJacs_[cacheCounter][j][0]);
            else
              local_integral += weight_local_quadrature * (diffusive_flux[0] * reconstructionGradPhij[0]);

            // add entries
            local_matrix.add(j, i, local_integral);
          }
        }
      }
    }
  } // for


  // set unit rows for dirichlet dofs
  Dune::Multiscale::getConstraintsCoarse(coarseDiscreteFunctionSpace_).applyToOperator(global_matrix_);

  DSC_PROFILER.stopTiming("msfem.assembleMatrix");
  DSC_LOG_DEBUG << "Time to assemble and communicate MsFEM matrix: " << DSC_PROFILER.getTiming("msfem.assembleMatrix") << "ms\n";
}
#endif
